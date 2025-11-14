!==============================================================================
!
!  Reversible Jump MCMC Sampling with parallel tempering for 
!  receiver function, surface (fundamental mode Rayeleigh) wave dispersion, Rayleigh wave ellipticity, and magnetotelluric data inversion (including Buckingham porous media model)
!  with iterative full covariance matrix estimation
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, May 9 2014
!  jand@uvic.ca                       (250) 472 4026
!  http://web.uvic.ca~/jand/
!  Last change: May 9 2014
!
!  Based on Green 1995, Malinverno 2002, Bodin Sambridge 2009, 
!  Agostinetti Malinverno 2010, Dettmer etal 2010, 2012, 1013
!  Reversible Jump MCMC Sampling with parallel tempering for receiver function, surface (fundamental mode Rayeleigh) wave dispersion
!  Created plotting functions for RF and SWD inversion results
!
!  Pejman Shahsavari, University of Calgary, Oct 1 2023
!  pejman.shahssvari@ucalgary.ca
!  Added Rayeleigh wave ellipticity, magnetotelluric, iterative covariance matrix estimation   
!  Made changes and created extra necessary plotting functions 
!
!==============================================================================

PROGRAM  RJMCMC_RF

!=======================================================================
USE RJMCMC_COM
USE NR
!USE SAC_I_O
USE ieee_arithmetic
IMPLICIT NONE
!INCLUDE 'raysum/params.h'
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: i,j,ipar,ipar2,ifr,ilay,ifreq,imcmc,ithin,isource,ikeep,ikeep2,iDres,is,ie,ik,iaz
INTEGER(KIND=IB)  :: idat,iang,it,ic,it2,ivo,isource1,isource2,io,ismp

TYPE (objstruc),DIMENSION(2)            :: objm     ! Objects on master node
TYPE (objstruc)                         :: obj      ! Object on slave node
TYPE (objstruc)                         :: objnew   ! Object on slave node
!TYPE (covstruc),ALLOCATABLE,DIMENSION(:):: cov      ! Structure for inverse data Covariance Matrices
REAL(KIND=RP)                           :: ran_uni  ! Uniform random number
REAL(KIND=RP)                           :: ran_nor  ! Normal random number
REAL(KIND=RP),DIMENSION(:),ALLOCATABLE  :: tmpvoro

!!---------------------------------------------------------------------
!!    Iterative covariance update
!!---------------------------------------------------------------------
INTEGER(KIND=IB) :: covIter_nsamples, CHAINTHIN_COVest_period                         
INTEGER(KIND=IB) :: iprocess, ios, stat, stat2, stat3 
CHARACTER(LEN=10) :: word

!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds

REAL(KIND=RP)               :: tstart, tend              ! Overall time 
REAL(KIND=RP)               :: tstart2, tend2            ! Time for one forward model computation
REAL(KIND=RP)               :: tstartsnd, tendsnd        ! Communication time
REAL(KIND=RP)               :: tstartcmp, tendcmp, tcmp  ! Forward computation time

icovIter = 0_IB

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )
OPEN(UNIT=20,FILE=filebasefile,STATUS='OLD',ACTION='READ')
READ(20,*) filebaselen
READ(20,*) filebase
CLOSE(20)

IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
IF(rank == src)WRITE(6,*) '~~~             Reversible Jump MCMC Sampling              ~~~  '
IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
IF(rank == src)WRITE(6,*) '...running on ',NTHREAD,' cores'

!!
!! Read parameter file and allocate global prior limits:
!!
CALL READPARFILE()
IF (rank==src) CALL PRINTPAR2()
CALL READCOVPARFile()

!! to prevent Bug
IF (I_RV==0) ICOV_iterUpdate_RV = 0
IF (I_SWD==0) ICOV_iterUpdate_SWD = 0
IF (I_ELL==0) ICOV_iterUpdate_ELL = 0
IF (I_MT==0) ICOV_iterUpdate_MT = 0

ICOViter_datasets = (/ICOV_iterUpdate_RV, ICOV_iterUpdate_SWD, ICOV_iterUpdate_ELL, ICOV_iterUpdate_MT /) 
cov_converged_datasets = 1
WHERE (ICOViter_datasets == 1) cov_converged_datasets = 0 

!IF (ICOV_iterUpdate_RV==1 .OR. ICOV_iterUpdate_SWD==1 .OR. ICOV_iterUpdate_ELL==1 .OR. ICOV_iterUpdate_MT==1)
!  ICOV_iterUpdate = 1
!END IF

IF ( ANY(ICOViter_datasets == 1) )  ICOV_iterUpdate = 1

IF (ICOV_iterUpdate==1) THEN
    cov_converged = .FALSE.
ELSE
    cov_converged = .TRUE.
    !!MAXcovIter = 0_IB
END IF
!! The original code was for fixed (non-iteratively updated) cov matrix. After update for iterative cov estimate, 
!! the original approach is retrived by setting (ICOV_iterUpdate=0 and/or cov_converged=.TRUE.) at the zeroth (and in this case the only) iteration.

!! Allocate objects for sampling
CALL ALLOC_OBJ(objm(1))
CALL ALLOC_OBJ(objm(2))
CALL ALLOC_OBJ(obj)
CALL ALLOC_OBJ(objnew)

CALL ALLOC_COVmat()

IF(rank == src)THEN
  !!
  !!  File to save posterior samples
  !!
  IF (IMAP==0) THEN
      OPEN(NEWUNIT(usample),FILE=samplefile,FORM='formatted',STATUS='REPLACE', &
      ACTION='WRITE')
      OPEN(NEWUNIT(ustep),FILE=stepsizefile,FORM='formatted',STATUS='REPLACE', &
      ACTION='WRITE')
      IF (ICOV_iterUpdate==1) OPEN(NEWUNIT(usample_covIter),FILE=samplefile_covIter,FORM='formatted',STATUS='REPLACE', &
      ACTION='WRITE')
  END IF

ENDIF
IF (rank==src) WRITE(*,*) 'NRF1, NMODE = NMODE_ELL', NRF1, NMODE, NMODE_ELL
ncount1 = NFPMX+NAP+3*NRF1+3*NRF1+2*NMODE+2*NMODE_ELL+1 !! 1 is added in the end for sdparMT (my comment)
ncount2 = NFPMX+1+3*NRF1+3*NRF1+2*NMODE+2*NMODE_ELL + 1
ncount3 = NRF2*NTIME + NMODE2*NDAT_SWD + NMODE_ELL2*NDAT_ELL + NMT2*2*NDAT_MT + 1

!ALLOCATE( sample(NKEEP,ncount1),tmpvoro(NFPMX) )
ALLOCATE( tmpvoro(NFPMX) )
ALLOCATE( tmpmap(ncount2) )
ALLOCATE( taper_dpred(NTIME) )
!sample       = 0._RP
ALLOCATE( icount(NTHREAD) )

!!
!!  Make MPI structure to match objstruc
!!
IF(IMAP == 0) THEN
  CALL MAKE_MPI_STRUC_SP(objm(1),objtype1)
  CALL MAKE_MPI_STRUC_SP(objm(2),objtype2)
  CALL MAKE_MPI_STRUC_SP(obj,objtype3)
END IF
!!------------------------------------------------------------------------
!!  Read in data
!!------------------------------------------------------------------------
CALL READDATA(obj)
!!------------------------------------------------------------------------
!!
!! Population / parallel tempering:
!!
!!------------------------------------------------------------------------
IF(NTHREAD > 1)THEN
  NPTCHAINS  = NTHREAD    ! Number of chains (equal to number of slaves)
  NT = NPTCHAINS-NPTCHAINS1+1
  ALLOCATE( NCHAINT(NT),beta_pt(NPTCHAINS) )
  NCHAINT = 1
  NCHAINT(1) = NPTCHAINS1

!!------------------------------------------------------------------------
!!  Set up tempering schedule
!!------------------------------------------------------------------------
  IF(rank == src)WRITE(6,*)'NT',NT
  IF(rank == src)WRITE(6,*)'NCHAINT',NCHAINT
  it2 = 1_IB
  DO it=1,NT
    DO ic=1,NCHAINT(it)
      beta_pt(it2) = 1._RP/dTlog**REAL(it-1_IB,RP)
      IF(rank == src)WRITE(6,218) 'Chain ',it2,':  beta = ',beta_pt(it2),'  T = ',1._RP/beta_pt(it2)
      it2 = it2 + 1
    ENDDO
  ENDDO
ENDIF
IF(rank == src)WRITE(6,*) ''
218 FORMAT(a6,i4,a10,F8.6,a6,F8.2)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!
!! Initialize random seeds on each core (Call RANDOM_SEED only once in the whole code. PARALLEL_SEED calls it)
!!
!CALL RANDOM_SEED


CALL PARALLEL_SEED()

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!CALL ellipticity_gpell_init(ELL_verbose)
!CALL dispersion_curve_init(ELL_verbose)

IF (ICOV_iterUpdate==1) THEN
    IF (ICOVest==1) THEN
        NsampleDres = 1_IB
        ALLOCATE( sampleDres(NsampleDres,ncount3) )
    ELSE IF (ICOVest==2) THEN
        usample_res_covIter = NEWUNIT()   !! required by all MPI processes
        !usample_res_covIter = 743         !! required by all MPI processes
       IF (rank==src) OPEN(UNIT=usample_res_covIter,FILE=samplefile_res_covIter,FORM='formatted',STATUS='REPLACE', &
                           ACTION='WRITE')
    END IF
END IF
stat  = 0_IB
stat2 = 0_IB
stat3 = 0_IB

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!
icovIter = 0_IB
DO      !! loop on the cov iterations 

IF ( (.NOT.cov_converged) .AND. (IMAP==0) ) THEN
    IF (rank==src) THEN
    WRITE(*,*) ''
    WRITE(*,*) 'Starting covariance iteration ', icovIter
    WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
    WRITE(*,*) ''
    END IF
END IF

!!
!! re-initialize necessary variables at the beginning of each cov iteration
!!
IF (icovIter>0_IB) THEN   
    !! variables iniatialized in RJMCMC_COM module
    raysumfail = 0_IB
    !kmin     = 0  !! changes inside READPARFILE()
    !kmax     = 0  !! changes inside READPARFILE()
    !NKEEP      = 1E1_IB   !! changes inside READPARFILE()
    !armxH      = 0.2_RP   !! changes inside READPARFILE()
    !armxV      = 0.2_RP   !! changes inside READPARFILE()
    !armxSWD    = 0.5_RP   !! changes inside READPARFILE()
    !armxELL    = 0.5_RP   !! changes inside READPARFILE()
    !modname = 'sample.geom'  !! doesn't change during the code
    ibirth  = 0
    ideath  = 0
    ibirths = 0
    ideaths = 0
    ioutside   = 0
    ireject    = 0
    iaccept = 0
    iaccept_delay = 0
    ireject_delay = 0
    i_sdpert = 0
    ishearfail = 0
    i_ref_nlay = 0
    iconv    = 0
    iconv2   = 0
    iconv3   = 0
    iarfail  = 0
    !ICHAINTHIN = 1E0_IB  !! changes inside READPARFILE()
    imcmc2 = 1
    !NFIELD = 62 !! doesn't change during the code
    obj%ireject_bd = 0
    obj%iaccept_bd = 0
    obj%ireject_bds = 0
    obj%iaccept_bds = 0
    obj%idxar = 0
    obj%gvoroidx = 0
    obj%S = 0
    objnew%ireject_bd = 0
    objnew%iaccept_bd = 0
    objnew%ireject_bds = 0
    objnew%iaccept_bds = 0  !! I have been coservative. Reinitializing objnew is probably not necessary
    objnew%idxar = 0
    objnew%gvoroidx = 0
    objnew%S = 0
    !! variables initialized by READPARFILE()
    isoflag = .FALSE.  !! inside ALLOC_RAYSUM()
    sdevm = 0._RP    
END IF
    !!------------------------------------- 


IF (ICOV_iterUpdate==0) THEN
    NKEEP2 = NKEEP
    IF (rank==src) ALLOCATE( sample(NKEEP2,ncount1) , STAT=stat )  
    units = usample
ELSEIF (icovIter==0_IB) THEN
    covIter_nsamples = covIter_zero_nsamples
    IF (rank==src) ALLOCATE( sample(covIter_nsamples,ncount1), STAT=stat )   
    NKEEP2 = NKEEP_covIter
    units = usample_covIter
    IF (ICOVest==2) THEN 
        NKEEP3 = NKEEP_covIter_res
        IF (rank==src) ALLOCATE( sample2(NKEEP3, ncount3), STAT=stat2 )   
        CHAINTHIN_COVest_period = CHAINTHIN_COVest_period_zeroIter
    END IF
    IF (iSAVEsample_only_zeroIter==1) iSAVEsample_covIter = 1    !! just to avoid a bug when (iSAVEsample_only_zeroIter==1) .AND. iSAVEsample_covIter = 0
ELSEIF (icovIter==1_IB) THEN
    IF (ALLOCATED(sample)) DEALLOCATE(sample)
    covIter_nsamples = covIter_period
    IF (rank==src) ALLOCATE( sample(covIter_nsamples,ncount1), STAT=stat )   
    IF (ICOVest==2) THEN
        CHAINTHIN_COVest_period = CHAINTHIN_COVest_period_nonzeroIter
        NsampleDres = covIter_zero_nsamples / CHAINTHIN_COVest_period_zeroIter   !!integer division 
        ALLOCATE( sampleDres( NsampleDres, ncount3  ), STAT=stat3  )
    END IF
    IF (iSAVEsample_only_zeroIter==1) iSAVEsample_covIter = 0
ELSEIF ( (ICOVest==2) .AND. (icovIter==2_IB) ) THEN
    IF (ALLOCATED(sampleDres)) DEALLOCATE(sampleDres)    
    NsampleDres = covIter_period / CHAINTHIN_COVest_period_nonzeroIter 
    ALLOCATE( sampleDres( NsampleDres, ncount3  ), STAT=stat3  )
END IF

IF (stat/=0) THEN   !! just process src can change stat and stat2. They have been initialized to 0 for all processes
    WRITE(*,*) 'Memory allocation for the array of posterior samples failed:'
    WRITE(*,*) 'at covariance iteration: ', icovIter 
    WRITE(*,*) 'while covariance did not converged.'  
    STOP
END IF

IF (stat2/=0) THEN   !! stat, stat2 and stat3 must just be used once here. If another variable for checking the STAT of ALLOCATE is necessary, declare a new one
    WRITE(*,*) 'Memory allocation for the saving array of residuals of posterior samples failed at cov iteration 0: '
    STOP
END IF

IF (stat3/=0) THEN
    WRITE(*,*) 'Memory allocation for the reading array of residuals of posterior samples failed in process: ', rank   
    WRITE(*,*) 'at covariance iteration: ', icovIter 
    STOP
END IF  

!!
!! Read mapfile to start (reads the starting parameter set):
!!
OPEN(UNIT=20,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) tmpmap
CLOSE(20)

obj%k       = INT(tmpmap(1),IB)
obj%NFP     = (obj%k * NPL)
obj%voro    = 0._RP
obj%voroidx = 0
obj%sdparR  = 0._RP
obj%sdparV  = 0._RP
obj%sdparT  = 0._RP
obj%arpar   = 0._RP
obj%sdparSWD= 0._RP
obj%arparSWD= 0._RP
obj%sdparELL= 0._RP
obj%arparELL= 0._RP
obj%sdparMT = 0._RP

obj%sdaveH  = 0._RP
obj%sdaveV  = 0._RP
obj%sdaveSWD= 0._RP
obj%sdaveELL= 0._RP
obj%sdaveMT = 0._RP

DO ivo = 1,obj%k
  ipar = (ivo-1)*NPL+2
  obj%voro(ivo,:) = tmpmap(ipar:ipar+NPL-1)
  obj%voroidx(ivo,:) = 1
  DO ipar2 = 1,NPL
    IF(obj%voro(ivo,ipar2) < -99._RP)THEN
      obj%voroidx(ivo,ipar2) = 0
    ELSE
      obj%voroidx(ivo,ipar2) = 1
    ENDIF
  ENDDO
ENDDO
IF(ICOV >= 1)THEN
  IF(I_RV == -1 .OR. I_RV == 1) obj%sdparR = tmpmap(NFPMX+1+1:NFPMX+1+NRF1)
  IF(I_RV == 1) obj%sdparV   = tmpmap(NFPMX+1+1+NRF1:NFPMX+1+2*NRF1)
  IF(I_T == 1) obj%sdparT   = tmpmap(NFPMX+1+1+2*NRF1:NFPMX+1+3*NRF1)
END IF
IF(ICOV_SWD >= 1) THEN
  IF(I_SWD == 1) obj%sdparSWD = tmpmap(NFPMX+1+1+3*NRF1:NFPMX+1+3*NRF1+NMODE)
END IF
IF(ICOV_ELL >= 1) THEN
  IF(I_ELL == 1) obj%sdparELL = tmpmap(NFPMX+1+3*NRF1+NMODE+1:NFPMX+1+3*NRF1+NMODE+NMODE_ELL)
END IF
IF(ICOV_MT >= 1) THEN
  IF(I_MT == 1) obj%sdparMT = tmpmap(NFPMX+1+3*NRF1+NMODE+NMODE_ELL+1)
END IF

IF(IAR == 1)THEN
  obj%arpar = tmpmap(NFPMX+1+3*NRF1+NMODE+NMODE_ELL+2:NFPMX+3*NRF1+NMODE+NMODE_ELL+3*NRF1+2)
  obj%arparSWD = tmpmap(NFPMX+1+3*NRF1+NMODE+NMODE_ELL+3*NRF1+2:NFPMX+3*NRF1+NMODE+NMODE_ELL+3*NRF1+NMODE+2)
  obj%arparELL = tmpmap(NFPMX+1+3*NRF1+NMODE+NMODE_ELL+3*NRF1+NMODE+2:NFPMX+3*NRF1+NMODE+NMODE_ELL+3*NRF1+NMODE+NMODE_ELL+2)
  obj%idxar = 1
  DO ipar = 1,NRF1
    IF(obj%arpar(ipar) < minlimar(ipar)) obj%idxar(ipar) = 0
  ENDDO
  obj%idxarSWD = 1
  DO ipar = 1,NMODE
    IF(obj%arparSWD(ipar) < minlimarSWD(ipar)) obj%idxarSWD(ipar) = 0
  ENDDO
  obj%idxarELL = 1
  DO ipar = 1,NMODE_ELL
    IF(obj%arparELL(ipar) < minlimarELL(ipar)) obj%idxarELL(ipar) = 0
  ENDDO
ENDIF

IF(I_VARPAR == 0)CALL INTERPLAYER_novar(obj)
IF(I_VARPAR == 1)CALL INTERPLAYER(obj)
IF(rank /= src)obj%beta = beta_pt(rank)

!ALLOCATE( icount(NTHREAD) )
icount = 0
tstart = MPI_WTIME()

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

IF(IMAP == 1)THEN
  IF(rank == src)THEN 
   WRITE(6,*) 'IMAP activated, exiting after computing replica for MAP.'
   CALL CHECKBOUNDS(obj)
   tstart2 = MPI_WTIME()
   !ALLOCATE(cur_modMT(NLMX,2)) !!!!!!!!!Pejman: for training dataset for the surrogate model
   !DO ic=1,100
     !CALL PRINTPAR(obj)
     !CALL PROPOSAL(obj,objnew,12,1,1.)
     !CALL PRINTPAR(objnew)
     !STOP
     IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
     IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)
     tend2 = MPI_WTIME()
     WRITE(6,*) 'time 1 = ',tend2-tstart2
   !ENDDO
   tend2 = MPI_WTIME()
   WRITE(6,*) 'time = ',tend2-tstart2
   WRITE(6,*) 'logL = ',obj%logL
   CALL SAVEREPLICA(obj)
   IF(ioutside == 1)WRITE(*,*)'FAILED in starting model'
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_FINALIZE( ierr )
  STOP
ELSE
  IF(rank == 1)THEN
    IF (icovIter==0_IB) THEN
        WRITE(6,*) 'Starting model:'
    ELSE
        WRITE(6,*) 'MAP model from previous cov iteration:'
    END IF
    CALL PRINTPAR(obj)
    CALL CHECKBOUNDS(obj)
    IF(ioutside == 1)WRITE(*,*)'FAILED in starting model'
    tstart2 = MPI_WTIME()
    IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
    IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)
    tend2 = MPI_WTIME()
    WRITE(6,*) 'logL = ',obj%logL
    WRITE(6,*) 'time = ',tend2-tstart2
    CALL CHECKBOUNDS(obj)
    IF(ioutside == 1)THEN 
      PRINT*,'Outside bounds.'
      STOP
    ENDIF
    tend2 = MPI_WTIME()
    WRITE(6,*) 'time for linrot = ',tend2-tstart2
  ENDIF
ENDIF
CALL FLUSH(6)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!
!!
!!  Need this call to update obj%res for the MAP model; necessary forcov estimation when ICOVest=1
!!

IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)  
IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)

!!  Updating covariance matrices       
!!
IF (icovIter > 0_IB)  THEN

    IF (ICOVest==1) THEN
        IF(I_RV==-1) sampleDres(1, 1:NRF2*NTIME) = obj%DresR(NRF1,:)
        IF(I_SWD==1) sampleDres(1, NRF2*NTIME+1:NRF2*NTIME+NMODE2*NDAT_SWD) = obj%DresSWD(NMODE,:) 
        IF(I_ELL==1) sampleDres(1, NRF2*NTIME+NMODE2*NDAT_SWD+1:NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL) = obj%DresELL(NMODE_ELL,:)
        IF(I_MT==1)  sampleDres(1, NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+1:NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+NMT2*2*NDAT_MT) = obj%DresMT
        sampleDres(1, ncount3) = 1._RP

    ELSEIF (ICOVest==2) THEN
        OPEN(UNIT=usample_res_covIter, FILE=samplefile_res_covIter, STATUS='OLD', ACTION='READ')
        DO ismp = 1, NsampleDres
            READ(usample_res_covIter,*) sampleDres(ismp, 1:ncount3)
        END DO
        CLOSE(usample_res_covIter)
        CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    END IF !! ICOVest
    
    IF (rank==src) THEN
        WRITE(*,*) ''
        WRITE(*,*) 'Start updating covariance matrices'
    END IF
    CALL UpdateCOV(obj) !! if covariance iterations are converged, cov_converged=.TRUE. inside UpdateCOV()
    IF (rank==src) THEN
        WRITE(*,*) 'Done updating covariance matrices'
        WRITE(*,*) ''
        IF ( (ICOVest==2) .AND. (.NOT.cov_converged) ) OPEN(UNIT=usample_res_covIter, FILE=samplefile_res_covIter, STATUS='REPLACE', ACTION='WRITE')
    END IF

END IF  !! icovIter>0

!!
!!  Need this call to initiate all chains with updated obj%logL by the new covariance matrix
!!
    IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
    IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)

IF (icovIter == MAXcovIter) cov_converged = .TRUE.  !! icovIter == MAXcovIter .AND. cov_converged = .FALSE.

icount(rank+1) = icount(rank+1) + 1

!------------------------------------------------------------------------
!
!          ************ RJMCMC Sampling ************
!
! -----------------------------------------------------------------------
IF(rank == src) THEN 

  WRITE(*,*) ' '
  IF (.NOT.cov_converged) THEN

    !IF (icovIter<MAXcovIter) THEN
    WRITE(6,*) 'At covariance iteration: ', icovIter 
    IF (icovIter==MAXcovIter) THEN
      WRITE(*,*) 'Covariance did not converge within the MAX number of iterations' 
      WRITE(*,*) 'Posterior sampling is done via the last covariance estimate (last iteration)'
      !cov_converged = .TRUE. 
    END IF

  ELSE
   
    IF (ICOV_iterUpdate==1) THEN
      WRITE(*,*) 'Covariance converged after ', icovIter, 'number of iterations'
      WRITE(*,*) 'Now posterior sampling is started'
    END IF

  END IF  !! cov_converged
  WRITE(*,*) ' '
  !! print some info
  IF(icovIter>0_IB) THEN
      CALL PRINTPAR2()      
      IF(ICOV >= 1)THEN
          IF(I_RV == -1 .OR. I_RV ==1) WRITE(6,*) 'SD parameters:'
          IF(I_RV == -1 .OR. I_RV ==1) WRITE(6,206) 'sigma H   = ',obj%sdparR
          IF(I_RV == -1 .OR. I_RV ==1) WRITE(6,206) 'sigma V   = ',obj%sdparV
          IF(I_RV == -1 .OR. I_RV ==1) WRITE(6,206) 'sigma T   = ',obj%sdparT
      END IF
      IF(ICOV_SWD >= 1) THEN
          IF(I_SWD == 1) WRITE(6,206) 'sigma SWD = ',obj%sdparSWD
      END IF
      IF(ICOV_ELL >= 1) THEN
          IF(I_ELL == 1) WRITE(6,206) 'sigma ELL = ',obj%sdparELL
      END IF
      IF(ICOV_MT >= 1) THEN
          IF(I_MT == 1) WRITE(6,206) 'sigma MT = ',obj%sdparMT
      ENDIF
  END IF
  206 FORMAT(a,128ES12.3)

  WRITE(6,*) 'Starting RJMCMC sampling...'

END IF !! rank

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN

IF ( (ICOV_iterUpdate==1) .AND. cov_converged ) THEN
    IF (ALLOCATED(sample)) DEALLOCATE(sample)
    NKEEP2 = NKEEP                                
    ALLOCATE( sample(NKEEP2,ncount1) , STAT=stat )  
    CLOSE(units, IOSTAT=ios)
    IF (ALLOCATED(sample2)) DEALLOCATE(sample2)
    units = usample
END IF

IF (stat/=0) THEN   
    WRITE(*,*) 'Memory allocation for the array of posterior samples failed:'
    WRITE(*,*) 'at covariance iteration: ', icovIter 
    WRITE(*,*) 'while covariance just converged.'  
    STOP
END IF

imcmc1 = 1

ncswap = 0
ncswapprop = 0
ikeep = 1
ikeep2 = 1
is = 1
ie = NKEEP2
tsave1 = MPI_WTIME()
DO imcmc = 1,NCHAIN

  !!
  !!  Receiving samples from any two slave (no particular order -> auto load
  !!  balancing) and propose swap:
  !!
  DO ithin = 1,ICHAINTHIN
    CALL MPI_RECV(objm(1), 1, objtype1, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
    isource1 = status(MPI_SOURCE)  !! This saves the slave id for the following communication
    tstart = MPI_WTIME()
    CALL MPI_RECV(objm(2), 1, objtype2, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
    isource2 = status(MPI_SOURCE)  !! This saves the slave id for the following communication
    IF(IEXCHANGE == 1)CALL TEMPSWP_MH(objm(1),objm(2))
    !IF(IAR == 1)CALL UDATE_SDAVE(objm(1),objm(2))
    CALL MPI_SEND(objm(1), 1,objtype1, isource1, rank, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(objm(2), 1,objtype2, isource2, rank, MPI_COMM_WORLD, ierr)
    tend = MPI_WTIME()
  ENDDO

  CALL SAVESAMPLE(objm,ikeep,ikeep2,covIter_nsamples,CHAINTHIN_COVest_period,is,ie,isource1,isource2,REAL(tend-tstart,RP))

  IF (.NOT.cov_converged) THEN

      IF ( imcmc1 > covIter_nsamples ) THEN
          DO iprocess = 1,NTHREAD-1
              CALL MPI_RECV(objm(1), 1, objtype1, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
              isource1 = status(MPI_SOURCE)  
              CALL MPI_SEND(objm(1), 1,objtype1, isource1, 1000, MPI_COMM_WORLD, ierr)
          END DO
          WRITE(*,*) ' '
          WRITE(*,*) 'Covariance iteration ', icovIter, 'Finished' 
          WRITE(*,*) ' '
          EXIT
      END IF
 
  END IF

ENDDO
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!    WORKER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ELSE
tstart = MPI_WTIME()
DO imcmc = 1,NCHAIN

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!  EXPLORE CALLS

  DO ithin = 1,1
    !! BD MH for fixed complexity:
    IF(I_VARPAR == 0)CALL EXPLORE_MH_NOVARPAR(obj,objnew,obj%beta)

    !! BD MH for variable complexity:
    IF(I_VARPAR == 1)CALL EXPLORE_MH_VARPAR(obj,objnew,obj%beta)
  ENDDO
  !! MH for non partition parameters:
  !!IF (ISD == 1) THEN
    CALL EXPLORE_MH(obj,objnew,obj%beta)
  !!END IF

  !!  END EXPLORE CALLS
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

  tend = MPI_WTIME()
  obj%tcmp = REAL(tend-tstart,RP)
  CALL MPI_SEND(obj, 1,objtype3, src, rank, MPI_COMM_WORLD, ierr)
  CALL MPI_RECV(obj, 1,objtype3, src, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
  IF (status(MPI_TAG)==1000) EXIT
  tstart = MPI_WTIME()
  

!  t_chckpt2 = MPI_WTIME()
  !! Checkpoint every once in a while...
!  IF(NINT(t_chckpt2-t_chckpt1)>TCHCKPT)THEN
!    PRINT*,rank,'dumps object at',obj%logL
!    !PRINT*,NINT(t_chckpt2-t_chckpt1),TCHCKPT
!    CALL SAVEOBJECT(obj)
!    t_chckpt1 = MPI_WTIME()
!  ENDIF

  !! Scale proposals (diminishing adaptation after Handbook of MCMC)
!  obj%iacceptvoro(:,1:NCMT)  = SUM(iacceptcmt_buf,3)
!  obj%iproposevoro(:,1:NCMT) = SUM(iproposecmt_buf,3)
!  obj%iacceptvoro(:,NCMT+1:NPL)  = SUM(iacceptloc_buf,3)
!  obj%iproposevoro(:,NCMT+1:NPL) = SUM(iproposeloc_buf,3)
!  !IF(rank == 1)PRINT*, 'propose:',obj%iproposevoro
!  !IF(rank == 1)PRINT*, 'accept:',obj%iacceptvoro
!  IF(IADAPT == 1)THEN
!    !IF(MOD(imcmc,NBUF) == 0)THEN
!    IF(iadaptcmt == 1)THEN
!      iadaptcmt = 0
!      DO ivo = 1,obj%k
!      DO ipar = 1,NCMT
!        IF(REAL(SUM(iacceptcmt_buf(ivo,ipar,:)),RP)/REAL(SUM(iproposecmt_buf(ivo,ipar,:)),RP) < 0.20_RP)THEN
!          obj%pertsd(ivo,ipar) = obj%pertsd(ivo,ipar)*MAX(0.90_RP,1._RP-1._RP/SQRT(REAL(imcmc2,RP)))
!        ENDIF
!        IF(REAL(SUM(iacceptcmt_buf(ivo,ipar,:)),RP)/REAL(SUM(iproposecmt_buf(ivo,ipar,:)),RP) > 0.30_RP)THEN
!          obj%pertsd(ivo,ipar) = obj%pertsd(ivo,ipar)/MAX(0.90_RP,1._RP-1._RP/SQRT(REAL(imcmc2,RP)))
!        ENDIF
!      ENDDO
!      ENDDO
!    ENDIF
!  ENDIF
  !! MCMC counter for diminishing adaptation
!  imcmc2 = imcmc2 + 1_IB
ENDDO !! imcmc
ENDIF !! MPI ENDIF process rank

IF (cov_converged) EXIT   !! it contains (icovIter>MAXcovIter) due to the earlier condition that when  (icovIter==MAXcovIter) cov_converged=.True.

IF (rank==src) THEN
    WRITE(*,*), ''
    WRITE(*,*), 'Start updating MAP file'
    !WRITE(*,*) (sample(i,1), i=1,5)
    CALL UPDATE_MAPfile(sample, covIter_nsamples, ncount1)
    WRITE(*,*), 'Done updating MAP file'
    WRITE(*,*), ''
    IF (ICOVest==2) CLOSE(usample_res_covIter)
END IF
CALL MPI_BARRIER( MPI_COMM_WORLD,ierr )

icovIter = icovIter + 1_IB
!!IF (icovIter > MAXcovIter) EXIT     !! also could DO WHILE (icovIter<=MAXcovIter) for the loop on cov iterations in the start
!CLOSE(6)
END DO !! cov iterations
CALL MPI_FINALIZE( ierr )

201 FORMAT(200F12.4)
202 FORMAT(I8,2F16.6,I8,F16.6,2I8,1F13.3,I8,1F13.3)
203 FORMAT(A119)
204 FORMAT(A56)
205 FORMAT(10F10.4)
209 FORMAT(A26,A40)
210 FORMAT(A26,20I4)
211 FORMAT(20000ES12.4)
212 FORMAT(I4,A8,F12.4,A5,I3)
213 FORMAT(I4,I6,A8,F12.4,A5,I3,A14,F12.4,A11,I3)
214 FORMAT(a21,100I4)
215 FORMAT(a21,F8.4)
216 FORMAT(a21,F8.4)
217 FORMAT(a21,F8.4)

END PROGRAM RJMCMC_RF
!!==============================================================================

SUBROUTINE EXPLORE_MH(obj,objnew1,beta_mh)
!!==============================================================================
!!
!! This is MH for all hierarchical and non partition parameters (no birth or 
!! death)
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iwhich,idel,ipar,ilay,ivo
INTEGER(KIND=IB)                 :: ncra,iparst
TYPE(objstruc)                   :: obj,objnew1,objnew2
!TYPE (covstruc),DIMENSION(NRF1)  :: cov
REAL(KIND=RP)                    :: logPLratio,logy,logq1_1,logq1_2
REAL(KIND=RP)                    :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                    :: znew,beta_mh,Lr1,Lr2,PROP
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxrand
INTEGER(KIND=IB)                 :: ipick
INTEGER(KIND=IB),DIMENSION(NLMX) :: idxpick
REAL(KIND=RP)                    :: logarp
INTEGER(KIND=IB)                 :: arptype,choose

objnew1 = obj
!!
!! Do Metropolis-Hastings on data-error standard deviations
!!
IF(1 == 1)THEN
IF(ICOV >= 1)THEN
  !!
  !! Sample standard deviation of R and V components
  !!
  IF(I_RV >= 1 .OR. I_RV == -1)THEN
  IF (ISD_RV == 1)THEN
  !! Perturb std devs with .10 probability
  CALL RANDOM_NUMBER(ran_uni_ar)
  IF(ran_uni_ar>=0.50_RP)THEN
    DO ipar = 1,NRF1
      CALL PROPOSAL_SDH(obj,objnew1,ipar)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
    ENDDO
    IF(I_RV >= 1)THEN
      DO ipar = 1,NRF1
        CALL PROPOSAL_SDV(obj,objnew1,ipar)
        IF(ioutside == 0)THEN
          IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
          IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
          logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
            objnew1 = obj
            ireject = ireject + 1
          ELSE
            obj = objnew1
            iaccept = iaccept + 1
          ENDIF
        ELSE
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
        ENDIF
      ENDDO
    ENDIF
    IF(I_T == 1)THEN
      DO ipar = 1,NRF1
        CALL PROPOSAL_SDT(obj,objnew1,ipar)
        IF(ioutside == 0)THEN
          IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
          IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
          logPLratio = (objnew1%logL - obj%logL)*beta_mh
          CALL RANDOM_NUMBER(ran_uni)
          IF(ran_uni >= EXP(logPLratio))THEN
            objnew1 = obj
            ireject = ireject + 1
          ELSE
            obj = objnew1
            iaccept = iaccept + 1
          ENDIF
        ELSE
          objnew1 = obj
          ireject = ireject + 1
          ioutside = 0
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  ENDIF  !! ISD_RV
  ENDIF  !! IRV
END IF   !! ICOV>=1

IF(ICOV_SWD >= 1) THEN
  !!
  !! Sample standard deviation of SWD
  !!
  IF(I_SWD == 1)THEN
  IF (ISD_SWD == 1)THEN
  CALL RANDOM_NUMBER(ran_uni_ar)
  IF(ran_uni_ar>=0.10_RP)THEN
    DO ipar = 1,NMODE
      CALL PROPOSAL_SDSWD(obj,objnew1,ipar)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
    ENDDO
  ENDIF
  ENDIF !! ISD_SWD
  ENDIF !! I_SWD
END IF  !! ICOV_SWD >=1

IF (ICOV_ELL >= 1) THEN
  !!
  !! Sample standard deviation of ELL
  !!
  IF(I_ELL == 1)THEN
  IF (ISD_ELL == 1)THEN
  CALL RANDOM_NUMBER(ran_uni_ar)
  IF(ran_uni_ar>=0.10_RP)THEN
    DO ipar = 1,NMODE_ELL
      CALL PROPOSAL_SDELL(obj,objnew1,ipar)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
    ENDDO
  ENDIF
  ENDIF !! ISD_ELL
  ENDIF !! I_ELL
END IF  !! ICOV_ELL >= 1

IF(ICOV_MT >= 1) THEN
  !!
  !! sample standard deviation of MT
  !!
  IF(I_MT == 1)THEN
  IF (ISD_MT == 1)THEN
  CALL RANDOM_NUMBER(ran_uni_ar)
  IF(ran_uni_ar>=0.10_RP)THEN
      CALL PROPOSAL_SDMT(obj,objnew1,1)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
  ENDIF  
  ENDIF !! ISD_MT
  ENDIF !! I_MT
ENDIF ! ICOV_MT >= 1
!!
!! Do Metropolis-Hastings on autoregressive model
!!
IF(IAR == 1)THEN
!IF(I_RV == 1)THEN
!  !! Perturb AR model with .25 probability
!  CALL RANDOM_NUMBER(ran_uni_ar)
!  IF(ran_uni_ar>=0.25_RP)THEN
!    DO ipar = 1,3*NRF1
!      IF(obj%idxar(ipar) == 0)THEN
!        !! Propose birth
!        arptype = 1
!        logarp = LOG(0.5_RP)
!      ELSE
!        CALL RANDOM_NUMBER(ran_uni_ar)
!        IF(ran_uni_ar>=0.5_RP)THEN
!          !! Propose death
!          arptype = 2
!          logarp = LOG(2._RP)
!        ELSE
!          !! Propose perturb
!          arptype = 3
!          logarp = 0._RP
!        ENDIF
!      ENDIF
!      CALL PROPOSAL_AR(obj,objnew1,ipar,arptype)
!      IF(ioutside == 0)THEN
!        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
!        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
!        !!logPLratio = (objnew1%logL - obj%logL)*beta_mh
!        !! Input Birth Death AR here:
!        logPLratio = logarp + (objnew1%logL - obj%logL)*beta_mh
!        CALL RANDOM_NUMBER(ran_uni)
!        IF(ran_uni >= EXP(logPLratio))THEN
!          objnew1 = obj
!          ireject = ireject + 1
!        ELSE
!          obj = objnew1
!          iaccept = iaccept + 1
!        ENDIF
!      ELSE
!        objnew1 = obj
!        ireject = ireject + 1
!        ioutside = 0
!      ENDIF
!      i_sdpert = 0
!    ENDDO
!  ENDIF
!ENDIF ! I_RV if
!!
!! Sampling RF case:
!!
IF(I_RV == -1)THEN
  !! Perturb AR model with .25 probability
  CALL RANDOM_NUMBER(ran_uni_ar)
  IF(ran_uni_ar>=0._RP)THEN
    DO ipar = 1,NRF1
      IF(obj%idxar(ipar) == 0)THEN
        !! Propose birth
        arptype = 1
        logarp = LOG(0.5_RP)
      ELSE
        CALL RANDOM_NUMBER(ran_uni_ar)
        IF(ran_uni_ar>=0.5_RP)THEN
          !! Propose death
          arptype = 2
          logarp = LOG(2._RP)
        ELSE
          !! Propose perturb
          arptype = 3
          logarp = 0._RP
        ENDIF
      ENDIF
      CALL PROPOSAL_ARRF(obj,objnew1,ipar,arptype)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        !!logPLratio = (objnew1%logL - obj%logL)*beta_mh
        !! Input Birth Death AR here:
        logPLratio = logarp + (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
      i_sdpert = 0
    ENDDO
  ENDIF
ENDIF ! I_RV if
ENDIF ! AR if
!!
!! Do Metropolis-Hastings on autoregressive model SWD
!!
IF(IAR == 1)THEN
IF(I_SWD == 1)THEN
  !! Perturb AR model with .25 probability
  CALL RANDOM_NUMBER(ran_uni_ar)
  !IF(ran_uni_ar>=0.25_RP)THEN
    DO ipar = 1,NMODE
      IF(obj%idxarSWD(ipar) == 0)THEN
        !! Propose birth
        arptype = 1
        logarp = LOG(0.5_RP)
      ELSE
        CALL RANDOM_NUMBER(ran_uni_ar)
        IF(ran_uni_ar>=0.5_RP)THEN
          !! Propose death
          arptype = 2
          logarp = LOG(2._RP)
        ELSE
          !! Propose perturb
          arptype = 3
          logarp = 0._RP
        ENDIF
      ENDIF
      CALL PROPOSAL_ARSWD(obj,objnew1,ipar,arptype)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        !!logPLratio = (objnew1%logL - obj%logL)*beta_mh
        !! Input Birth Death AR here:
        logPLratio = logarp + (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
      i_sdpert = 0
    ENDDO
  !ENDIF ! ran_uni_ar if
ENDIF ! I_SWD if
ENDIF ! AR if
!!
!! Do Metropolis-Hastings on autoregressive model ELL
!!
IF(IAR == 1)THEN
IF(I_ELL == 1)THEN
  !! Perturb AR model with .25 probability
  CALL RANDOM_NUMBER(ran_uni_ar)
  !IF(ran_uni_ar>=0.25_RP)THEN
    DO ipar = 1,NMODE_ELL
      IF(obj%idxarELL(ipar) == 0)THEN
        !! Propose birth
        arptype = 1
        logarp = LOG(0.5_RP)
      ELSE
        CALL RANDOM_NUMBER(ran_uni_ar)
        IF(ran_uni_ar>=0.5_RP)THEN
          !! Propose death
          arptype = 2
          logarp = LOG(2._RP)
        ELSE
          !! Propose perturb
          arptype = 3
          logarp = 0._RP
        ENDIF
      ENDIF
      CALL PROPOSAL_ARELL(obj,objnew1,ipar,arptype)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        !!logPLratio = (objnew1%logL - obj%logL)*beta_mh
        !! Input Birth Death AR here:
        logPLratio = logarp + (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
      i_sdpert = 0
    ENDDO
  !ENDIF ! ran_uni_ar if
ENDIF ! I_SWD if
ENDIF ! AR if

ENDIF !! 1 == 2 if
END SUBROUTINE EXPLORE_MH
!!==============================================================================

SUBROUTINE EXPLORE_MH_NOVARPAR(obj,objnew1,beta_mh)
!!==============================================================================
!!
!! This is MH for partition parameters with Birth and Death for the case 
!! of fixed complexity in layers
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iwhich,idel,ipar,ilay,ivo
INTEGER(KIND=IB)                 :: ncra,iparst
TYPE(objstruc)                   :: obj,objnew1,objnew2
!TYPE (covstruc),DIMENSION(NRF1)  :: cov
REAL(KIND=RP)                    :: logPLratio,logPratio
REAL(KIND=RP)                    :: logy,logq1_1,logq1_2
REAL(KIND=RP)                    :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                    :: znew,beta_mh,Lr1,Lr2,PROP
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxrand
INTEGER(KIND=IB)                 :: ipick
INTEGER(KIND=IB),DIMENSION(NLMX) :: idxpick
REAL(KIND=RP)                    :: logarp
INTEGER(KIND=IB)                 :: arptype,choose

objnew1 = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)
i_bd = 0
!! Do BIRTH-DEATH MCMC with 0.5 probability
IF(kmin /= kmax)THEN
  !! Perturbing k:
  CALL RANDOM_NUMBER(ran_unik)
  IF(obj%k == kmax)THEN  
    !! If k == kmax, no birth allowed, so 1/3 death, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 2
  ELSEIF(obj%k == kmin)THEN  
    !! If k == kmin, no death allowed, so 1/3 birth, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 1
  ELSE
    !! 1/3 birth, 1/3 death, 1/3 stay local
    IF((ran_unik <= 0.3333_RP))i_bd = 1
    IF(ran_unik > 0.6666_RP)i_bd = 2
  ENDIF
  IF(i_bd == 1) CALL BIRTH_FULL(obj,objnew1)
  IF(i_bd == 2) CALL DEATH_FULL(obj,objnew1)
  IF(obj%k /= objnew1%k)THEN
    !!
    !! If k changed, check BD acceptance
    !!
    CALL CHECKBOUNDS(objnew1)
    IF(ioutside == 0)THEN
      IF(ISMPPRIOR == 0) CALL LOGLHOOD(objnew1,1)
      IF(ISMPPRIOR == 1) CALL LOGLHOOD2(objnew1)
      !logPLratio =  (objnew1%logL - obj%logL)*beta_mh
      logPratio  = objnew1%logPr
      logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh

      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
        objnew1 = obj
        obj%ireject_bd = obj%ireject_bd + 1
      ELSE
        obj = objnew1
        obj%iaccept_bd = obj%iaccept_bd + 1
      ENDIF
    ELSE
      objnew1 = obj
      obj%ireject_bd = obj%ireject_bd + 1
      ioutside = 0
    ENDIF
  ENDIF  ! k-change if
ENDIF
!!
!! Carry out MH sweep
!!
!!
!! Do Metropolis-Hastings update on c, rho, alpha
!!
IF(1 == 1)THEN
DO ivo = 1,obj%k
  idxrand = 0
  idxrand(1:NPL) = RANDPERM(NPL)
  DO ipar = 1,NPL
    iwhich = idxrand(ipar)
    IF(ivo == 1 .AND. iwhich == 1) CYCLE
    IF(obj%voroidx(ivo,iwhich) == 1)THEN
      CALL PROPOSAL(obj,objnew1,ivo,iwhich,1._RP)
      CALL CHECKBOUNDS2(objnew1,ivo,iwhich)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,1)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        !logPLratio = (objnew1%logL - obj%logL)*beta_mh
        logPratio  = objnew1%logPr
        logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF  !! MH if
      ELSE !! outside before delayed rejection
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF !! ioutside if
    ENDIF !!  voroidx if
  ENDDO
ENDDO
ENDIF !! 1 == 2 if
END SUBROUTINE EXPLORE_MH_NOVARPAR
!!==============================================================================
SUBROUTINE EXPLORE_MH_VARPAR(obj,objnew1,beta_mh)
!!==============================================================================
!!
!! This is MH for partition parameters with Birth and Death for the case 
!! of variable complexity in layers (most general case with 2 different 
!! birth and death steps.
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iwhich,idel,ipar,ilay,ivo
TYPE(objstruc)                   :: obj,objnew1,objnew2
!TYPE (covstruc),DIMENSION(NRF1)  :: cov
REAL(KIND=RP)                    :: logPLratio
REAL(KIND=RP)                    :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                    :: beta_mh
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxrand
INTEGER(KIND=IB)                 :: ipick
INTEGER(KIND=IB),DIMENSION(NLMX) :: idxpick
REAL(KIND=RP)                    :: logarp
INTEGER(KIND=IB)                 :: arptype,choose

objnew1 = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)
i_bd = 0
!! Do BIRTH-DEATH MCMC with 0.5 probability
IF(kmin /= kmax)THEN
  !! Perturbing k:
  CALL RANDOM_NUMBER(ran_unik)
  IF(obj%k == kmax)THEN
    !! If k == kmax, no birth allowed, so 1/3 death, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 2
  ELSEIF(obj%k == kmin)THEN
    !! If k == kmin, no death allowed, so 1/3 birth, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 1
  ELSE
    !! 1/3 birth, 1/3 death, 1/3 stay local
    IF(ran_unik <= 0.3333_RP)i_bd = 1
    IF(ran_unik > 0.6666_RP)i_bd = 2
  ENDIF
  IF(i_bd == 1) CALL BIRTH_VARPAR(obj,objnew1)
  IF(i_bd == 2) CALL DEATH_VARPAR(obj,objnew1)
  IF(obj%k /= objnew1%k)THEN
    !!
    !! If k changed, check BD acceptance
    !!
    CALL CHECKBOUNDS(objnew1)
    IF(ioutside == 0)THEN
      IF(ISMPPRIOR == 0) CALL LOGLHOOD(objnew1,1)
      IF(ISMPPRIOR == 1) CALL LOGLHOOD2(objnew1)
      logPLratio =  (objnew1%logL - obj%logL)*beta_mh

      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
        objnew1 = obj
        obj%ireject_bd = obj%ireject_bd + 1
      ELSE
        obj = objnew1
        obj%iaccept_bd = obj%iaccept_bd + 1
      ENDIF
    ELSE
      objnew1 = obj
      obj%ireject_bd = obj%ireject_bd + 1
      ioutside = 0
    ENDIF
  ENDIF  ! k-change if
ENDIF
IF(IBD_SINGLE == 1)THEN
IF(obj%k > 1)THEN
  objnew1 = obj
  i_bds = 0
  !! Pick random node (except for first one):
  idxpick = 0
  idxpick(1:obj%k-1) = RANDPERM(obj%k-1)
  ipick = idxpick(1)+1

  CALL RANDOM_NUMBER(ran_unik)
  IF(SUM(obj%voroidx(ipick,:)) == 2)THEN
    !! If node has min parameters (2), 1/3 birth, 2/3 stay local
    IF(ran_unik<=0.3333_RP) i_bds = 1
  ELSEIF(SUM(obj%voroidx(ipick,:)) == NPL)THEN
    !! If node has max parameters, 1/3 death, 2/3 stay local
    IF(ran_unik<=0.3333_RP) i_bds = 2
  ELSE
    !! 1/3 birth, 1/3 death, 1/3 stay local
    IF(ran_unik <= 0.3333_RP) i_bds = 1
    IF(ran_unik >  0.6666_RP) i_bds = 2
  ENDIF
  IF(i_bds == 1) CALL BIRTH_SINGLE(obj,objnew1,ipick)
  IF(i_bds == 2) CALL DEATH_SINGLE(obj,objnew1,ipick)
  CALL CHECKBOUNDS(objnew1)
  IF(ioutside == 0)THEN
    IF(ISMPPRIOR == 0) CALL LOGLHOOD(objnew1,1)
    IF(ISMPPRIOR == 1) CALL LOGLHOOD2(objnew1)
    !! Likelihood: ratio
    logPLratio =  (objnew1%logL - obj%logL)*beta_mh

    CALL RANDOM_NUMBER(ran_uni)
    IF(ran_uni >= EXP(logPLratio))THEN
      objnew1 = obj
      obj%ireject_bds = obj%ireject_bds + 1
    ELSE
      obj = objnew1
      obj%iaccept_bds = obj%iaccept_bds + 1
    ENDIF
  ELSE
    objnew1 = obj
    obj%ireject_bds = obj%ireject_bds + 1
    ioutside = 0
  ENDIF
ENDIF ! k > 1 if
ENDIF ! BD single if
!!
!! Carry out MH sweep
!!
!!
!! Do Metropolis-Hastings update on c, rho, alpha
!!
IF(1 == 1)THEN
DO ivo = 1,obj%k
  idxrand = 0
  idxrand(1:NPL) = RANDPERM(NPL)
  DO ipar = 1,NPL
    iwhich = idxrand(ipar)
    !! Skip the topmost node position, since it's fixed at 0.
    IF(ivo == 1 .AND. iwhich == 1) CYCLE
    IF(obj%voroidx(ivo,iwhich) == 1)THEN
      CALL PROPOSAL(obj,objnew1,ivo,iwhich,1._RP)
      CALL CHECKBOUNDS2(objnew1,ivo,iwhich)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,1)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF  !! MH if
      ELSE !! outside before delayed rejection
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF !! ioutside if
    ENDIF !!  voroidx if
  ENDDO
ENDDO
ENDIF !! 1 == 2 if
END SUBROUTINE EXPLORE_MH_VARPAR
!=======================================================================

SUBROUTINE DEATH_FULL(obj,objnew)
!!=======================================================================
!!
!! This DEATH is the reverse for Birth from prior.
!!
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB)                 :: idel,ivo,ipar
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxdeath
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort
REAL(KIND=RP)                    :: zdel,zj,zjp1

objnew = obj
objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!! 1) Pick random node:
idxdeath = 0
idxdeath(1:obj%k-1) = RANDPERM(obj%k-1)
idxdeath = idxdeath + 1
idel = idxdeath(1)
zdel = objnew%voro(idel,1)
zj   = objnew%voro(idel-1,1)
IF(idel == obj%k)THEN
  zjp1 = hmx
ELSE
  zjp1 = objnew%voro(idel+1,1)
ENDIF

!! 2) Save location and parameters of node
!objnew%gvoroidx = 1

!! 3) delete node and re-interpolate Voronoi model
objnew%voro(idel,:) = 0._RP
objnew%voroidx(idel,:) = 0
tmpsort = 0._DRP
tmpsort = REAL(objnew%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),objnew%voroidx(1:obj%k,:))
objnew%voro(1:obj%k,:) = REAL(tmpsort,RP)
!CALL QSORTC2D(objnew%voro(1:obj%k,:),objnew%voroidx(1:obj%k,:))
!! 4) Quicksort stores deleted array (all zeros) first.
!!    Hence, move all one up unless deleted node is last node.
objnew%voro(1:objnew%k,:) = objnew%voro(2:obj%k,:)
objnew%voroidx(1:objnew%k,:) = objnew%voroidx(2:obj%k,:)
objnew%voro(objnew%k+1,:) = 0._RP
objnew%voroidx(objnew%k+1,:) = 0

IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)

!PRINT*,'zj',zj
!PRINT*,'znew',zdel,idel
!PRINT*,'zjp1',zjp1
!!
!! Compute prior for even-numbered order statistics (Green 1995)
!! znew is new interface, zj (j) is interface above and zjp1 is interface below (j+1)
!!
IF(ENOS == 0 .AND. IPOIPR == 0)THEN
  objnew%logPr = 0._RP
ELSEIF(ENOS == 1 .AND. IPOIPR == 0)THEN
!! Apply only even-numbered order statistics in prior:
  objnew%logPr = 2._RP*LOG(hmx-hmin)-LOG(2._RP*obj%k*(2._RP*obj%k+1._RP)) + &
                 LOG(zjp1-zj)-LOG(zdel-zj)-LOG(zjp1-zdel)
ELSEIF(ENOS == 0 .AND. IPOIPR == 1)THEN
!! Apply only Poisson prior:
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))
ELSEIF(ENOS == 1 .AND. IPOIPR == 1)THEN
!! Apply Poisson prior with ENOS:
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))+2._RP*LOG(hmx-hmin)- &
                 LOG(2._RP*obj%k*(2._RP*obj%k+1._RP)) + &
                 LOG(zjp1-zj)-LOG(zdel-zj)-LOG(zjp1-zdel)
ENDIF

END SUBROUTINE DEATH_FULL
!=======================================================================

SUBROUTINE BIRTH_FULL(obj,objnew)
!!
!! This birthes a full node with all parameters. Birth is based on the 
!! background value at the position of the new node which are then 
!! perturbed with a Gaussian proposal.
!!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iznew, ipar, ivo, nnew
TYPE(objstruc)                   :: obj,objnew
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=RP),DIMENSION(obj%k+1) :: ztmp
REAL(KIND=RP)                    :: znew,zj,zjp1
INTEGER(KIND=IB),DIMENSION(NPL-1):: idxran
INTEGER(KIND=IB),DIMENSION(NPL)  :: voroidx

objnew = obj
objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!!
!! Sample No. new parameters for node from [1,NPL-1]
!! THIS CODE IS FIXED TO ALWAYS BITH ALL PARAMETERS!
!!
!idxran = RANDPERM(NPL-1)
!nnew = idxran(1)
nnew = NPL-1
idxran = 0
idxran = RANDPERM(NPL-1)+1  !! First nnew elements give index of 
                            !! parameters to be perturbed.
!!
!! voroidx identifies live parameters
!! THIS CODE IS FIXED TO ALWAYS BITH ALL PARAMETERS!
!!
voroidx = 1
!voroidx(1) = 1
!voroidx(idxran(1:nnew)) = 1

!!
!! Draw new z
!!
CALL RANDOM_NUMBER(ran_uni)
znew = maxpert(1)*ran_uni
!!
!! Insert new node at bottom of stack...
!!
!! 1) Sample depth from uniform prior:
objnew%voro(objnew%k,1)   = znew
objnew%voroidx(objnew%k,:)= voroidx
!objnew%gvoroidx           = voroidx(2:NPL)

!! 2) Sample new node parameters from prior
DO ipar = 2,NPL
  !!
  !! Gaussian proposal
  !!
  IF(objnew%voroidx(objnew%k,ipar) == 1)THEN
    CALL RANDOM_NUMBER(ran_uni)
    objnew%voro(objnew%k,ipar) = minlim(ipar) + maxpert(ipar)*ran_uni
  ENDIF
ENDDO
!IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
!IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)

CALL INTERPLAYER_novar(objnew)

ztmp = objnew%ziface(1:objnew%k) - znew
!!
!! Find new interface index
!!
iznew = 0
DO ivo = 1,obj%k
   IF(ztmp(ivo) == 0._RP) iznew = ivo+1
ENDDO
zj = objnew%voro(iznew-1,1)
IF(iznew > obj%k)THEN
  zjp1 = hmx
ELSE
  zjp1 = objnew%voro(iznew+1,1)
ENDIF

!PRINT*,'zj',zj
!PRINT*,'znew',znew,iznew
!PRINT*,'zjp1',zjp1
!iznew = iznew + 1
!!
!! Compute prior for even-numbered order statistics (Green 1995)
!! znew is new interface, zj (j) is interface above and zjp1 is interface below (j+1)
!!
IF(ENOS == 0 .AND. IPOIPR == 0)THEN
  objnew%logPr = 0._RP
ELSEIF(ENOS == 1 .AND. IPOIPR == 0)THEN
  objnew%logPr = LOG(2._RP*obj%k+2._RP)+LOG(2._RP*obj%k+3._RP)-2._RP*LOG(hmx-hmin) + &
                 LOG(znew-zj)+LOG(zjp1-znew)-LOG(zjp1-zj)
ELSEIF(ENOS == 0 .AND. IPOIPR == 1)THEN
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))
ELSEIF(ENOS == 1 .AND. IPOIPR == 1)THEN
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))+LOG(2._RP*obj%k+2._RP)+ &
                 LOG(2._RP*obj%k+3._RP)-2._RP*LOG(hmx-hmin) + &
                 LOG(znew-zj)+LOG(zjp1-znew)-LOG(zjp1-zj)
ENDIF


RETURN
END SUBROUTINE BIRTH_FULL
!=======================================================================

SUBROUTINE DEATH_VARPAR(obj,objnew)
!!=======================================================================
!!
!! This DEATH is the reverse for Birth from prior.
!!
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB)                 :: idel,ivo,ipar
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxdeath
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort

objnew = obj
objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!! 1) Pick random node:
idxdeath = 0
idxdeath(1:obj%k-1) = RANDPERM(obj%k-1)
idxdeath = idxdeath + 1
idel = idxdeath(1)

!! 2) Save location and parameters of node
objnew%gvoroidx = obj%voroidx(idel,2:NPL)

!! 3) delete node and re-interpolate Voronoi model
objnew%voro(idel,:) = 0._RP
objnew%voroidx(idel,:) = 0

tmpsort = 0._DRP
tmpsort = REAL(objnew%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),objnew%voroidx(1:obj%k,:))
objnew%voro(1:obj%k,:) = REAL(tmpsort,RP)
!CALL QSORTC2D(objnew%voro(1:obj%k,:),objnew%voroidx(1:obj%k,:))

!! 4) Quicksort stores deleted array (all zeros) first.
!!    Hence, move all one up unless deleted node is last node.
objnew%voro(1:objnew%k,:) = objnew%voro(2:obj%k,:)
objnew%voroidx(1:objnew%k,:) = objnew%voroidx(2:obj%k,:)
objnew%voro(objnew%k+1,:) = 0._RP
objnew%voroidx(objnew%k+1,:) = 0

!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
END SUBROUTINE DEATH_VARPAR
!=======================================================================

SUBROUTINE BIRTH_VARPAR(obj,objnew)
!!
!! This birthes a full node with variable No. parameters. Birth is 
!! based on the background value at the position of the new node which 
!! are then perturbed with a Gaussian proposal.
!!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)              :: iznew, ipar, ivo, nnew
TYPE(objstruc)                :: obj,objnew
REAL(KIND=RP)                 :: znew,ran_uni
REAL(KIND=RP),DIMENSION(obj%nunique):: ztmp
INTEGER(KIND=IB),DIMENSION(NPL-1):: idxran
INTEGER(KIND=IB),DIMENSION(NPL):: voroidx

objnew = obj
objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!!
!! Sample No. new parameters for node from [1,NPL-1]
!!
idxran = RANDPERM(NPL-1)
nnew = idxran(1)
idxran = 0
idxran = RANDPERM(NPL-1)+1  !! First nnew elements give index of 
                            !! parameters to be perturbed.
!!
!! voroidx identifies live parameters
!!
voroidx = 0
voroidx(1) = 1
voroidx(idxran(1:nnew)) = 1

!!
!! Draw new z
!!
CALL RANDOM_NUMBER(ran_uni)
znew = maxpert(1)*ran_uni
ztmp = obj%ziface(1:obj%nunique) - znew
!!
!! Insert new node at bottom of stack...
!!
!! 1) Sample depth from uniform prior:
objnew%voro(objnew%k,1)   = znew
objnew%voroidx(objnew%k,:)= voroidx
objnew%gvoroidx           = voroidx(2:NPL)

!! 2) Sample new node parameters from prior
DO ipar = 2,NPL
  !!
  !! Gaussian proposal
  !!
  IF(objnew%voroidx(objnew%k,ipar) == 1)THEN
    CALL RANDOM_NUMBER(ran_uni)
    objnew%voro(objnew%k,ipar) = minlim(ipar) + maxpert(ipar)*ran_uni
  ENDIF
ENDDO
!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
RETURN
END SUBROUTINE BIRTH_VARPAR
!!=======================================================================

SUBROUTINE DEATH_SINGLE(obj,objnew,ipick)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB)                 :: ipar,itmp,ipick,idel
INTEGER(KIND=IB)                 :: NPLIVE
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NPL)  :: idxtmp,idxran

!! 1) Pick random live parameter
idxtmp = 0
itmp = 1
DO ipar = 2,NPL
  IF(obj%voroidx(ipick,ipar) == 1)THEN
    idxtmp(itmp) = ipar
    NPLIVE = itmp
    itmp = itmp + 1
  ENDIF
ENDDO
idxran = 0
idxran(1:NPLIVE) = RANDPERM(NPLIVE)
idel = idxtmp(idxran(1))

!! 2) delete parameter on node
objnew%voro(ipick,idel) = 0._RP
objnew%voroidx(ipick,idel) = 0

!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
!PRINT*,'DEATH:',idel,ipick
!PRINT*,obj%voro(ipick,:)
!PRINT*,obj%voroidx(ipick,:)
!PRINT*,objnew%voro(ipick,:)
!PRINT*,objnew%voroidx(ipick,:)
!PRINT*,''

END SUBROUTINE DEATH_SINGLE
!!=======================================================================

SUBROUTINE BIRTH_SINGLE(obj,objnew,ipick)
!!
!! This birthes a single parameter onto node. Birth is based on the 
!! background value at the position of the new node which is then
!! perturbed with a Gaussian proposal.
!!
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)               :: ipar,itmp,ipick,inew
INTEGER(KIND=IB)               :: NPFREE
TYPE(objstruc)                 :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NPL):: idxtmp,idxran
REAL(KIND=RP)                  :: ran_uni

!! 1) Randomly pick free parameter type
idxtmp = 0
itmp = 1
DO ipar = 2,NPL
  IF(obj%voroidx(ipick,ipar) == 0)THEN
    idxtmp(itmp) = ipar
    NPFREE = itmp
    itmp = itmp + 1
  ENDIF
ENDDO
idxran = 0
idxran(1:NPFREE) = RANDPERM(NPFREE)
inew = idxtmp(idxran(1))

!! 2) Sample new value from prior
CALL RANDOM_NUMBER(ran_uni)
objnew%voro(ipick,inew) = minlim(inew) + maxpert(inew)*ran_uni
objnew%voroidx(ipick,inew) = 1

!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
!PRINT*,'BIRTH:',inew
!PRINT*,obj%voro(ipick,:)
!PRINT*,obj%voroidx(ipick,:)
!PRINT*,objnew%voro(ipick,:)
!PRINT*,objnew%voroidx(ipick,:)
!PRINT*,''
RETURN
END SUBROUTINE BIRTH_SINGLE
!=======================================================================
SUBROUTINE TEMPSWP_MH(obj1,obj2)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
TYPE(objstruc)              :: obj1
TYPE(objstruc)              :: obj2
TYPE(objstruc)              :: objtmp1,objtmp2
REAL(KIND=RP)               :: ran_uni1
REAL(KIND=RP)               :: logratio
REAL(KIND=RP)               :: betaratio

betaratio = obj2%beta-obj1%beta
logratio  = betaratio*(obj1%logL-obj2%logL)
CALL RANDOM_NUMBER(ran_uni1)
IF(ran_uni1 <= EXP(logratio))THEN
  !! ACCEPT SWAP
!  IF(obj(ic1idx)%k /= obj(ic2idx)%k)THEN
!    PRINT*,'ACCEPTED',EXP(logratio*betaratio),ran_uni2
!    WRITE(ulog,204)'ic1:',ic1idx,obj(ic1idx)%k,logP1,obj(ic1idx)%logL,logratio,betaratio,EXP(logratio*betaratio),ran_uni2
!    WRITE(ulog,204)'ic2:',ic2idx,obj(ic2idx)%k,logP2,obj(ic2idx)%logL
!  ENDIF
  objtmp1       = obj1
  objtmp2       = obj2
  obj1   = objtmp2
  obj2   = objtmp1
  !! Temperature does not swap
  obj1%beta = objtmp1%beta
  obj2%beta = objtmp2%beta
  !! BD acceptance counters do not swap
!  obj1%iaccept_bd = objtmp1%iaccept_bd
!  obj2%iaccept_bd = objtmp2%iaccept_bd
!  obj1%ipropose_bd = objtmp1%ipropose_bd
!  obj2%ipropose_bd = objtmp2%ipropose_bd
  !! Voro acceptance counters do not swap
!  obj1%iacceptvoro = objtmp1%iacceptvoro
!  obj2%iacceptvoro = objtmp2%iacceptvoro
!  obj1%iproposevoro = objtmp1%iproposevoro
!  obj2%iproposevoro = objtmp2%iproposevoro
  !! Step sizes do not swap
!  obj1%pertsd = objtmp1%pertsd
!  obj2%pertsd = objtmp2%pertsd
  ncswap    = ncswap+1_IB
ENDIF
ncswapprop    = ncswapprop+1_IB

RETURN
END SUBROUTINE TEMPSWP_MH
!=======================================================================

SUBROUTINE PROPOSAL(obj,objnew,ivo,iwhich,factor)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,iwhich,ipar
TYPE(objstruc) :: obj,objnew,objtmp
REAL(KIND=RP)  :: ran_uni, factor
REAL(KIND=RP)                    :: zp,zj,zjp1,zjm1

objnew = obj
objnew%logPr = 0._RP
!!
!! CAUCHY proposal
!!
CALL RANDOM_NUMBER(ran_uni)
!! Cauchy:
IF(iwhich /= 1)THEN
  objnew%voro(ivo,iwhich) = obj%voro(ivo,iwhich) + fact/factor*pertsd(iwhich)*TAN(PI2*(ran_uni-0.5_RP))
ELSE
  !!
  !! Compute prior for even-numbered order statistics (Green 1995) where 
  !! znew is perturbed interface, zj (j) is original interface, 
  !! zjp1 is interface below (j+1) and zjm1 is interface above.
  !!
  IF(ENOS == 0)THEN
    !! CAUCHY proposal
    CALL RANDOM_NUMBER(ran_uni)
    objnew%voro(ivo,iwhich) = obj%voro(ivo,iwhich) + fact/factor*pertsd(iwhich)*TAN(PI2*(ran_uni-0.5_RP))
    objnew%logPr = 0._RP
  ELSE
    CALL RANDOM_NUMBER(ran_uni)
    zj = obj%voro(ivo,iwhich)
    zjm1 = obj%voro(ivo-1,iwhich)
    IF(ivo == obj%k)THEN
      zjp1 = hmx
    ELSE
      zjp1 = obj%voro(ivo+1,iwhich)
    ENDIF
    !! sample uniform:
    zp = zjm1+ran_uni*(zjp1-zjm1)
    objnew%voro(ivo,iwhich) = zp
    !! Apply even-numbered order statistics in prior:
    objnew%logPr = LOG(zjp1-zp)+LOG(zp-zjm1)-LOG(zjp1-zj)-LOG(zj-zjm1)
!    PRINT*,''
!    PRINT*,'ic',iwhich
!    PRINT*,'zjm1',zjm1
!    PRINT*,'zj',zj
!    PRINT*,'zp',zp
!    PRINT*,'zjp1',zjp1
!    PRINT*,'logPr',objnew%logPr
!    PRINT*,''
  ENDIF
ENDIF
!!
!! If node position comes up as negative, INTERPLAYER breaks...
!!
IF(iwhich == 1)THEN
  objnew%voro(ivo,iwhich) = ABS(objnew%voro(ivo,iwhich))
ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
RETURN
END SUBROUTINE PROPOSAL
!=======================================================================

SUBROUTINE PROPOSAL_AR(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arpar(iwhich) = ran_uni*(maxlimar(iwhich)-minlimar(iwhich))+minlimar(iwhich)
  objnew%idxar(iwhich) = 1
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arpar(iwhich) = minlimar(iwhich)-1._RP
  objnew%idxar(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arpar(iwhich) = obj%arpar(iwhich) + pertarsd(iwhich)*ran_nor
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_AR
!=======================================================================

SUBROUTINE PROPOSAL_ARRF(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arpar(iwhich) = ran_uni*(maxlimar(iwhich)-minlimar(iwhich))+minlimar(iwhich)
  objnew%idxar(iwhich) = 1
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arpar(iwhich) = minlimar(iwhich)-1._RP
  objnew%idxar(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arpar(iwhich) = obj%arpar(iwhich) + pertarsd(iwhich)*ran_nor
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_ARRF
!=======================================================================

SUBROUTINE PROPOSAL_ARSWD(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arparSWD(iwhich) = ran_uni*(maxlimarSWD(iwhich)-minlimarSWD(iwhich))+minlimarSWD(iwhich)
  objnew%idxarSWD(iwhich) = 1
  IF(((objnew%arparSWD(iwhich) - minlimarSWD(iwhich)) < 0._RP).OR. &
     ((maxlimarSWD(iwhich) - objnew%arparSWD(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arparSWD(iwhich) = minlimarSWD(iwhich)-1._RP
  objnew%idxarSWD(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arparSWD(iwhich) = obj%arparSWD(iwhich) + pertarsdSWD(iwhich)*ran_nor
  IF(((objnew%arparSWD(iwhich) - minlimarSWD(iwhich)) < 0._RP).OR. &
     ((maxlimarSWD(iwhich) - objnew%arparSWD(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_ARSWD
!=======================================================================

SUBROUTINE PROPOSAL_ARELL(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arparELL(iwhich) = ran_uni*(maxlimarELL(iwhich)-minlimarELL(iwhich))+minlimarELL(iwhich)
  objnew%idxarELL(iwhich) = 1
  IF(((objnew%arparELL(iwhich) - minlimarELL(iwhich)) < 0._RP).OR. &
     ((maxlimarELL(iwhich) - objnew%arparELL(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arparELL(iwhich) = minlimarELL(iwhich)-1._RP
  objnew%idxarELL(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arparELL(iwhich) = obj%arparELL(iwhich) + pertarsdELL(iwhich)*ran_nor
  IF(((objnew%arparELL(iwhich) - minlimarELL(iwhich)) < 0._RP).OR. &
     ((maxlimarELL(iwhich) - objnew%arparELL(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_ARELL
!=======================================================================

SUBROUTINE PROPOSAL_SDH(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparR(iwhich) = obj%sdparR(iwhich) + pertsdsd(iwhich)*ran_nor
IF(((objnew%sdparR(iwhich) - minlimsd(iwhich)) < 0._RP).OR.((maxlimsd(iwhich) - objnew%sdparR(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDH
!=======================================================================

SUBROUTINE PROPOSAL_SDV(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparV(iwhich) = obj%sdparV(iwhich) + pertsdsd(iwhich)*ran_nor
IF(((objnew%sdparV(iwhich) - minlimsd(iwhich)) < 0._RP).OR. & 
   ((maxlimsd(iwhich) - objnew%sdparV(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDV
!=======================================================================

SUBROUTINE PROPOSAL_SDT(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparT(iwhich) = obj%sdparT(iwhich) + pertsdsd(iwhich)*ran_nor
IF(((objnew%sdparT(iwhich) - minlimsd(iwhich)) < 0._RP).OR.((maxlimsd(iwhich) - objnew%sdparT(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDT
!=======================================================================

SUBROUTINE PROPOSAL_SDSWD(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparSWD(iwhich) = obj%sdparSWD(iwhich) + pertsdsdSWD(iwhich)*ran_nor
IF(((objnew%sdparSWD(iwhich) - minlimsdSWD(iwhich)) < 0._RP).OR. &
   ((maxlimsdSWD(iwhich) - objnew%sdparSWD(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDSWD

!=======================================================================

SUBROUTINE PROPOSAL_SDELL(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparELL(iwhich) = obj%sdparELL(iwhich) + pertsdsdELL(iwhich)*ran_nor
IF(((objnew%sdparELL(iwhich) - minlimsdELL(iwhich)) < 0._RP).OR. &
   ((maxlimsdELL(iwhich) - objnew%sdparELL(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDELL

!==============================================================================
!==============================================================================

SUBROUTINE PROPOSAL_SDMT(obj,objnew,iwhich)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparMT(iwhich) = obj%sdparMT(iwhich) + pertsdsdMT(iwhich)*ran_nor
IF(((objnew%sdparMT(iwhich) - minlimsdMT(iwhich)) < 0._RP).OR. &
   ((maxlimsdMT(iwhich) - objnew%sdparMT(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDMT

!==============================================================================
SUBROUTINE CHECKBOUNDS(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)               :: ivo,ipar,ilay,ncra
INTEGER(KIND=IB)               :: ih,OS,ihamrej
TYPE(objstruc)                 :: obj
REAL(KIND=RP)                  :: vspmin,vspmax,zhere
REAL(KIND=RP),DIMENSION(NFPMX2):: par

DO ilay = 1,obj%nunique
   IF(hmin > obj%hiface(ilay))ioutside = 1
   IF(maxlim(1) < obj%ziface(ilay))ioutside = 1
ENDDO

DO ivo = 1,obj%k
!! Starts at ivo = 2, since 1st voro node is fixed and gives half space if k = 0
  IF(ivo > 1)THEN
    IF((obj%voro(ivo,1) < 0._RP).OR.(obj%voro(ivo,1) > maxlim(1)))THEN
      ioutside = 1
      IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=', &
                     ivo,1,minlim(1),maxlim(1),obj%voro(ivo,1)
    ENDIF
  ENDIF
  DO ipar = 2,NPL
    IF(obj%voroidx(ivo,ipar) == 1)THEN
    IF(((obj%voro(ivo,ipar) - minlim(ipar)) < 0._RP).OR. & 
      ((maxlim(ipar) - obj%voro(ivo,ipar)) < 0._RP))THEN
      ioutside = 1
      IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=',ivo,ipar,&
                     minlim(ipar),maxlim(ipar),obj%voro(ivo,ipar)
    ENDIF
    ENDIF
  ENDDO 
ENDDO 
201 FORMAT(a27,2I4,3F18.4)

RETURN
END SUBROUTINE CHECKBOUNDS
!=======================================================================

SUBROUTINE CHECKBOUNDS2(obj,ivo,ipar)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)               :: ivo,ipar,ilay,ncra
INTEGER(KIND=IB)               :: ih,OS,ihamrej
TYPE(objstruc)                 :: obj
REAL(KIND=RP)                  :: vspmin,vspmax,zhere
REAL(KIND=RP),DIMENSION(NFPMX2):: par

DO ilay = 1,obj%nunique
   IF(hmin > obj%hiface(ilay))ioutside = 1
   IF(maxlim(1) < obj%ziface(ilay))ioutside = 1
ENDDO

!! Starts at ivo = 2, since 1st voro node is fixed and gives half space if k = 0
IF(ivo > 1)THEN
  IF((obj%voro(ivo,1) < 0._RP).OR.(obj%voro(ivo,1) > maxlim(1)))THEN
    ioutside = 1
    IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=', &
                   ivo,1,minlim(1),maxlim(1),obj%voro(ivo,1)
  ENDIF
ENDIF
IF(obj%voroidx(ivo,ipar) == 1)THEN
  IF(((obj%voro(ivo,ipar) - minlim(ipar)) < 0._RP).OR. & 
    ((maxlim(ipar) - obj%voro(ivo,ipar)) < 0._RP))THEN
    ioutside = 1
    IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=',ivo,ipar,&
                   minlim(ipar),maxlim(ipar),obj%voro(ivo,ipar)
  ENDIF
ENDIF
201 FORMAT(a27,2I4,3F18.4)

RETURN
END SUBROUTINE CHECKBOUNDS2
!=======================================================================
!SUBROUTINE UPDATE_SDAVE(obj1,obj2)
!!==============================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: irf
!TYPE(objstruc)   :: obj1
!TYPE(objstruc)   :: obj2
!
!DO irf = 1,NRF1
!  sdbuf(1,isdbuf)  = obj1%sdparR()
!ENDDO
!DO irf = 1,NRF1
!  obj1%sdaveH(irf)   = SUM(sdbuf(1,:))/REAL(NBUF,RP)
!  obj1%sdaveV(irf)   = SUM(sdbuf(2,:))/REAL(NBUF,RP)
!  obj1%sdaveSWD(irf) = SUM(sdbuf(3,:))/REAL(NBUF,RP)
!  obj2%sdaveH(irf)   = SUM(sdbuf(1,:))/REAL(NBUF,RP)
!  obj2%sdaveV(irf)   = SUM(sdbuf(2,:))/REAL(NBUF,RP)
!  obj2%sdaveSWD(irf) = SUM(sdbuf(3,:))/REAL(NBUF,RP)
!ENDDO
!
!RETURN
!END SUBROUTINE UPDATE_SDAVE
!=======================================================================

SUBROUTINE PARALLEL_SEED()
!!
!!  Ensure unique random seed for each CPU
!!
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed1
REAL(KIND=RP) :: ran_uni

CALL RANDOM_SEED
CALL RANDOM_SEED(SIZE=iseedsize)
ALLOCATE( iseed1(iseedsize) )
IF(ISETSEED == 1)THEN
   iseed1 = (/2303055,     2435432,     5604058,     4289794,     3472290, &
      7717070,      141180,     3783525,     3087889,     4812786,     3028075, &
      3712062,     6316731,      436800,     7957708,     2055697,     1944360, &
      1222992,     7537775,     7769874,     5588112,     7590383,     1426393, &
      1753301,     7681841,     2842400,     4411488,     7304010,      497639, &
      4978920,     5345495,      754842,     7360599,     5776102/)
   CALL RANDOM_SEED(PUT=iseed1)
ELSE
   CALL RANDOM_SEED(GET=iseed1)
!   IF(rank==src)WRITE(6,*) 'Master seed:',iseed1
ENDIF

ALLOCATE( iseed(iseedsize), rseeds(iseedsize,NTHREAD), iseeds(iseedsize,NTHREAD) )

iseed = 0
rseeds = 0._RP
iseeds = 0
IF(rank == src)THEN
   CALL RANDOM_NUMBER(rseeds)
   iseeds = -NINT(rseeds*1000000._RP)
ENDIF
DO i = 1,iseedsize
   CALL MPI_BCAST( iseeds(i,:), NTHREAD, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
ENDDO
iseed = iseeds(:,rank+1)

!!
!! Write seeds to seed logfile:
!!
IF(rank == src)THEN
   OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
   ACTION='WRITE',POSITION='REWIND',RECL=1024)
   WRITE(50,*) 'Rank: ',rank
   WRITE(50,201) iseed
   WRITE(50,*) ''
   CLOSE(50)
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
DO i = 1,NTHREAD-1
   IF(rank == i)THEN
      OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
      ACTION='WRITE',POSITION='APPEND',RECL=1024)
      WRITE(50,*) 'Rank: ',rank
      WRITE(50,201) iseed
      WRITE(50,*) ''
      CLOSE(50)
   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
ENDDO
CALL RANDOM_SEED(PUT=iseed)

DO i = 1,100
   CALL RANDOM_NUMBER(ran_uni)
ENDDO
DO i = 1,CEILING(ran_uni*10000._RP)
   CALL RANDOM_NUMBER(ran_uni)
ENDDO

201   FORMAT(50I10)
END SUBROUTINE PARALLEL_SEED
!=======================================================================

SUBROUTINE SAVESAMPLE(objm,ikeep,ikeep2,covIter_nsamples,CHAINTHIN_COVest_period,is,ie,isource1,isource2,tmaster)
!!=======================================================================
!!
!! Exchanging and saving posterior samples
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc),DIMENSION(2)     :: objm     ! Objects on master node
INTEGER(KIND=IB)                 :: ipar,ipar2,ivo,ic,imcmc,ikeep,ikeep2,is,ie,jidx
INTEGER(KIND=IB)                 :: CHAINTHIN_COVest_period
INTEGER(KIND=IB)                 :: isource,isource1,isource2
REAL(KIND=RP),DIMENSION(NFPMX)   :: tmpvoro  ! Temporary Voronoi cell array for saving
!INTEGER(KIND=RP),DIMENSION(NFPMX):: tmpprop,tmpacc
REAL(KIND=RP),DIMENSION(NPL)     :: tmp  ! Temporary Voronoi cell array for saving
REAL(KIND=RP) :: tmaster
INTEGER(KIND=IB), INTENT(IN)     :: covIter_nsamples
LOGICAL                          :: write_to_file

!!
!! Save object into global sample array
!!
DO ic = 1,2 
  IF ( (.NOT.cov_converged) .AND. (imcmc1 > (covIter_nsamples)) )  RETURN !! this is required in addition to the similar statement in the main function due to the loop on ic (to protect out of range subscript for the array sample)
  IF(ic == 1)isource = isource1
  IF(ic == 2)isource = isource2
  IF(objm(ic)%beta > 1._RP/dTlog)THEN
    !!
    !! Write to stdout
    !!
    IF(MOD(imcmc1,20*NKEEP2) == 0)THEN
      tsave2 = MPI_WTIME()
                        WRITE(*,213)'          iaccept_bd = ',objm(ic)%iaccept_bd,objm(ic)%beta
      IF(IEXCHANGE == 1)WRITE(*,215)'       T swap accept = ',REAL(ncswap,RP)/REAL(ncswapprop,RP),isource1,isource2
                        WRITE(*,216)'Total time for block = ',tsave2-tsave1,'s'
      WRITE(*,*) ''
      WRITE(*,203) '     imcmc1,         logL,  k,  iacc_bd, iacc_bds, source, time(slave)'
      WRITE(*,203) '----------------------------------------------------------------------'
      tsave1 = MPI_WTIME()
    ENDIF
    tmpvoro = 0._RP
     DO ivo = 1,objm(ic)%k
       ipar = (ivo-1)*NPL+1
       DO ipar2 = 1,NPL
         IF(objm(ic)%voroidx(ivo,ipar2)  == 1)THEN
           tmpvoro(ipar) = objm(ic)%voro(ivo,ipar2)
         ELSE
           tmpvoro(ipar) = -100._RP
         ENDIF
         ipar = ipar + 1
       ENDDO
!      tmpprop(ipar:ipar+NPL-1) = objm(ic)%iproposevoro(ivo,:)
!      tmpacc(ipar:ipar+NPL-1)  = objm(ic)%iacceptvoro(ivo,:)
     ENDDO
!    sample(ikeep,:) =  (/ objm(ic)%logL,objm(ic)%logPr,objm(ic)%lognorm,REAL(objm(ic)%k,RP),tmpvoro,&
!                          objm(ic)%arpar,objm(ic)%sdpar,&
!                          REAL(tmpacc,RP)/REAL(tmpprop,RP),REAL(objm(ic)%iaccept_bd,RP),&
!                          REAL(objm(ic)%ipropose_bd,RP),REAL(i_bd,RP),objm(ic)%tcmp,REAL(isource,RP) /)
     !WRITE(*,*) 'logL= ', objm(ic)%logL
     sample(ikeep,:) =  (/ objm(ic)%logL, objm(ic)%logPr, objm(ic)%tcmp, REAL(objm(ic)%k,RP), & ! 4 parameters
                           tmpvoro,objm(ic)%sdparR,objm(ic)%sdparV,objm(ic)%sdparT,objm(ic)%sdparSWD,objm(ic)%sdparELL,objm(ic)%sdparMT,  &
                           objm(ic)%arpar,objm(ic)%arparSWD,objm(ic)%arparELL, &
                           REAL(iaccept,RP)/REAL(iaccept+ireject,RP),REAL(objm(ic)%iaccept_bd,RP),&
                           REAL(objm(ic)%ireject_bd,RP),REAL(objm(ic)%iaccept_bds,RP),REAL(ic,RP),REAL(rank,RP) /)
    !!----------------------------------------------------------------------------------------
    IF (.NOT.cov_converged) THEN
    IF ( (ICOVest==2) .AND. (MOD(imcmc1,CHAINTHIN_COVest_period)==0) ) THEN
        IF(I_RV==-1) sample2(ikeep2, 1:NRF2*NTIME) = objm(ic)%DresR(NRF1,:)
        IF(I_SWD==1) sample2(ikeep2, NRF2*NTIME+1:NRF2*NTIME+NMODE2*NDAT_SWD) = objm(ic)%DresSWD(NMODE,:) 
        IF(I_ELL==1) sample2(ikeep2, NRF2*NTIME+NMODE2*NDAT_SWD+1:NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL) = objm(ic)%DresELL(NMODE_ELL,:)
        IF(I_MT==1)  sample2(ikeep2, NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+1:NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+NMT2*2*NDAT_MT) = objm(ic)%DresMT
        sample2(ikeep2, ncount3) = LOG(pk(objm(ic)%k)) + objm(ic)%logL
        ikeep2 = ikeep2 + 1_IB
        IF(ikeep2 > NKEEP3)THEN
            DO jidx=1,NKEEP3
            WRITE(usample_res_covIter,207) sample2(jidx,:)
            ENDDO
            CALL FLUSH(usample_res_covIter)
            sample2  = 0._RP
            ikeep2   = 1
        END IF
    END IF
    END IF
    !!------------------------------------------------------------------------------------------
    IF(MOD(imcmc1,NKEEP2)==0)THEN
      WRITE(*,202) '       ',imcmc1,objm(ic)%logL,objm(ic)%k,objm(ic)%ireject_bd,objm(ic)%iaccept_bds,isource,objm(ic)%tcmp
    ENDIF
    ikeep = ikeep + 1_IB
    imcmc1 = imcmc1 + 1_IB
    !! during covariance iterations (i.e., while cov_converged==.FALSE.) ikeep and imcmc1 are the same

    !!
    !! Write to file
    !!
    IF (cov_converged) THEN
      write_to_file = (ikeep > NKEEP2)
      !IF(MOD(imcmc1-1,NKEEP2)==0) WRITE(*,*) 'write_to_file: ', write_to_file
      !IF(MOD(imcmc1-1,NKEEP2)==0) WRITE(*,*) 'is, ie : ', is, ie
      !IF(MOD(imcmc1-1,NKEEP2)==0) WRITE(*,*) 'units : ', units
    ELSE
      write_to_file = ( (iSAVEsample_covIter==1) .AND. (ikeep>1) .AND. (MOD(ikeep-1,NKEEP2)==0)  ) !!ikeep>1 is unnecessary due to previousely ikeep=ikeep+1 
      !IF(MOD(ikeep-1,NKEEP2)==0) WRITE(*,*) 'write_to_file: ', write_to_file
      !IF(MOD(ikeep-1,NKEEP2)==0) WRITE(*,*) 'is, ie : ', is, ie
      !IF(MOD(ikeep-1,NKEEP2)==0) WRITE(*,*) 'units : ', units
    END IF 

    IF(write_to_file)THEN
      DO jidx=is,ie
        WRITE(units,207) sample(jidx,:)
      ENDDO
      CALL FLUSH(units)
!      WRITE(ustep,207) (objm(ic)%pertsd(ivo,:),ivo=1,objm(ic)%k)
!      CALL FLUSH(ustep)
      IF (cov_converged) THEN
        sample  = 0._RP
        ikeep   = 1
      ELSE
        is = is + NKEEP2  
        ie = ie + NKEEP2
      END IF
    ENDIF !! write_to_file

  ENDIF !! (objm(ic)%beta > 1._RP/dTlog)
ENDDO   !!ic
CALL FLUSH(6)
202 FORMAT(A3,I8,1F14.4,I4,I10,I10,I8,1f12.4)
203 FORMAT(A69)
207 FORMAT(500ES18.8)
213 FORMAT(a23,I,F8.4)
214 FORMAT(a23,I6,a9,I4)
215 FORMAT(a23,F8.4,I4,I4)
216 FORMAT(a23,F8.2,a)
RETURN
END SUBROUTINE 
!=======================================================================
SUBROUTINE SAVEREPLICA(obj)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i, ilayer
TYPE (objstruc)  :: obj      ! Best object

WRITE(6,*) 'Global best model:'
CALL PRINTPAR(obj)
WRITE(6,*) 'Global best logL = ',obj%logL

IF(I_RV >= 1)THEN
  OPEN(UNIT=50,FILE=predfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NRF1
    WRITE(50,208) obj%DpredR(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(50,208) obj%DpredV(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(50,208) obj%DpredT(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(50,208) obj%S(i,:)
  ENDDO
  WRITE(6,*)'Done writing predicted V and H components.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=obsfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NRF1
    WRITE(50,208) obj%DobsR(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(50,208) obj%DobsV(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(50,208) obj%DobsT(i,:)
  ENDDO
  WRITE(6,*)'Done writing observed V and H components.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=arfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NRF1
    WRITE(50,208) obj%DarR(i,:)
    WRITE(50,208) obj%DarV(i,:)
    WRITE(50,208) obj%DarT(i,:)
  ENDDO
  CLOSE(50)
ELSEIF(I_RV == -1)THEN
  OPEN(UNIT=50,FILE=predfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NRF1
    WRITE(50,208) obj%DpredR(i,:)
  ENDDO
  WRITE(6,*)'Done writing predicted RF.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=obsfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NRF1
    WRITE(50,208) obj%DobsR(i,:)
  ENDDO
  WRITE(6,*)'Done writing observed RF.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=arfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NRF1
    WRITE(50,208) obj%DarR(i,:)
  ENDDO
  CLOSE(50)

ENDIF
IF(I_SWD == 1)THEN
  OPEN(UNIT=50,FILE=predfileSWD,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE
    WRITE(50,208) obj%DpredSWD(i,:)
  ENDDO
  WRITE(6,*)'Done writing predicted SWD curve.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=obsfileSWD,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE
    WRITE(50,208) obj%DobsSWD(i,:)
  ENDDO
  WRITE(6,*)'Done writing observed SWD curve.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=arfileSWD,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE
    WRITE(50,208) obj%DarSWD(i,:)
  ENDDO
  CLOSE(50)

ENDIF
208 FORMAT(5000ES20.10)

IF(I_ELL == 1)THEN
  OPEN(UNIT=50,FILE=predfileELL,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE_ELL
    WRITE(50,208) obj%DpredELL(i,:)
  ENDDO
  WRITE(6,*)'Done writing predicted ELL curve.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=obsfileELL,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE_ELL
    WRITE(50,208) obj%DobsELL(i,:)
  ENDDO
  WRITE(6,*)'Done writing observed ELL curve.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=arfileELL,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE_ELL
    WRITE(50,208) obj%DarELL(i,:)
  ENDDO
  CLOSE(50)

ENDIF

IF(I_MT == 1)THEN
   OPEN(UNIT=50,FILE=predfileMT,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=8192)

   WRITE(50,208) obj%DpredMT
   
   WRITE(6,*)'Done writing predicted MT curve.'
   CLOSE(50)
 
   OPEN(UNIT=50,FILE=obsfileMT,FORM='formatted',STATUS='REPLACE', &
   ACTION='WRITE',RECL=8192)
   WRITE(50,208) obj%DobsMT
 
   WRITE(6,*)'Done writing observed MT curve.'
   CLOSE(50)
 
   !OPEN(UNIT=50,FILE=arfileSWD,FORM='formatted',STATUS='REPLACE', &
   !ACTION='WRITE',RECL=8192)
   !DO i = 1,NMODE
   !  WRITE(50,208) obj%DarSWD(i,:)
   !ENDDO
   !CLOSE(50)
 
   !!!!!!!!!!!!!!!  Pejman: training dataset for training the surrogate

   !OPEN(UNIT=500,FILE='/data/pejman/trainMT.dat',FORM='formatted',STATUS='REPLACE', &
   !ACTION='WRITE',RECL=8192)

   !WRITE(500,208) LOG10(obj%freqMT), cur_modMT(1:obj%nunique+1, 2), cur_modMT(1:obj%nunique, 1), obj%DpredMT  !!! thicknesses are in km
    !WRITE(500,208) LOG10(obj%freqMT), ( cur_modMT(ilayer,:), ilayer=1,obj%nunique  ), cur_modMT(obj%nunique+1,2:NPL), obj%DpredMT  !!!thicknesses are in km
   !!!!!!!!!!!!!!

ENDIF

RETURN
END SUBROUTINE SAVEREPLICA
!=======================================================================

FUNCTION RANDPERM(num)
!==============================================================================
USE data_type, ONLY : IB, RP
IMPLICIT NONE
INTEGER(KIND=IB), INTENT(IN) :: num
INTEGER(KIND=IB) :: numb, i, j, k
INTEGER(KIND=IB), DIMENSION(num) :: randperm
REAL(KIND=RP), DIMENSION(num) :: rand2
INTRINSIC RANDOM_NUMBER
CALL RANDOM_NUMBER(rand2)
DO i=1,num
   numb=1
   DO j=1,num
      IF (rand2(i) > rand2(j)) THEN
           numb=numb+1
      END IF
   END DO
   DO k=1,i-1
      IF (rand2(i) <= rand2(k) .AND. rand2(i) >= rand2(k)) THEN
           numb=numb+1
      END IF
   END DO
   randperm(i)=numb
END DO
RETURN
END FUNCTION RANDPERM
!====================================================================

SUBROUTINE GASDEVJ(harvest)
!====================================================================
USE nrtype
USE nr
IMPLICIT NONE
REAL(DP), INTENT(OUT) :: harvest
REAL(DP) :: rsq,v1,v2
REAL(DP), SAVE :: g
LOGICAL, SAVE :: gaus_stored=.FALSE.
IF (gaus_stored) THEN
   harvest=g
   gaus_stored=.FALSE.
ELSE
   DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1=2.0_DP*v1-1.0_DP
      v2=2.0_DP*v2-1.0_DP
      rsq=v1**2+v2**2
      IF (rsq > 0.0_DP .AND. rsq < 1.0_DP) EXIT
   END DO
   rsq=SQRT(-2.0_DP*LOG(rsq)/rsq)
   harvest=v1*rsq
   g=v2*rsq
   gaus_stored=.TRUE.
END IF
END SUBROUTINE GASDEVJ
!====================================================================

Function ASINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,ii,asinc
ii    = cmplx(0.,1.)
asinc = -ii*LOG(ii*z+SQRT(1.-z**2))
RETURN
END FUNCTION
!=======================================================================

Function COSC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,cosc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
cosc = CMPLX(COS(x)*COSH(y),-SIN(x)*SINH(y),RP)
RETURN
END FUNCTION
!=======================================================================

Function SINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,sinc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
sinc = CMPLX(SIN(x)*COSH(y),COS(x)*SINH(y),RP)
RETURN
END FUNCTION
!==============================================================================

FUNCTION CACOS(z)
!==============================================================================

USE DATA_TYPE
COMPLEX(KIND=RP) :: CACOS
COMPLEX(KIND=RP) :: z
REAL(KIND=RP) :: zrp1,zrm1,zi,zizi,a1,a2,a,b

!CACOS = -CMPLX(0._RP,1._RP,RP)*LOG(z+CMPLX(0._RP,1._RP,RP)*SQRT(1._RP-z*z))
!!
!! This version from IDL; much faster than above
!!
zrp1 = REAL(z,RP)+1._RP
zrm1 = zrp1-2._RP
zi = AIMAG(z)
zizi = zi*zi
a1 = 0.5_RP*SQRT(zrp1*zrp1 + zizi)
a2 = 0.5_RP*SQRT(zrm1*zrm1 + zizi)
a = a1+a2
b = a1- a2
IF(zi >= 0._RP)THEN
   CACOS = ACOS(b) - CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ELSE
   CACOS = ACOS(b) + CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ENDIF

RETURN
END FUNCTION CACOS
!==============================================================================

FUNCTION ASINH(x)
!==============================================================================
USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) :: ASINH
REAL(KIND=RP) :: x

ASINH = LOG(x+SQRT(x**2._RP+1))

RETURN
END FUNCTION ASINH
!==============================================================================

FUNCTION CSIN(z)
!==============================================================================
!! Complex sine (Jan's version)

USE DATA_TYPE
COMPLEX(KIND=RP) :: CSIN
COMPLEX(KIND=RP) :: z

CSIN =  (EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /CMPLX(0._RP,2._RP,RP)
RETURN
END FUNCTION CSIN
!==============================================================================

FUNCTION CTAN(z)
!==============================================================================
!! Complex TAN

USE DATA_TYPE
COMPLEX(KIND=RP) :: CTAN
COMPLEX(KIND=RP) :: z

CTAN =  -CMPLX(0._RP,1._RP,RP)*(EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /(EXP( CMPLX(0._RP,1._RP,RP)*z)+EXP(-CMPLX(0._RP,1._RP,RP)*z))
RETURN
END FUNCTION CTAN

!==============================================================================
function factorial (n) result (res)
!==============================================================================
 
implicit none
integer, intent (in) :: n
integer :: res
integer :: i
res = product ((/(i, i = 1, n)/))
end function factorial
!==============================================================================
function choose (n, k) result (res)
!==============================================================================
implicit none
integer, intent (in) :: n
integer, intent (in) :: k
integer :: res,factorial
res = factorial (n) / (factorial (k) * factorial (n - k))
end function choose
!=======================================================================

LOGICAL FUNCTION ISNAN(a)
!=======================================================================
USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) a
IF (a.NE.a) THEN
ISNAN = .TRUE.
ELSE
ISNAN = .FALSE.
END IF
RETURN
END
!=======================================================================

SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
!!=======================================================================
!!
!! Time domain convolution
!!
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: Nx,Ny,Nz
INTEGER(KIND=IB):: irow,icol1,icol2,ix1,ix2
REAL(KIND=RP)   :: x(Nx), y(Ny), z(Nz),xr(Nx)
REAL(KIND=RP)   :: XMAT(Nz,Ny)

XMAT = 0._RP
xr = x(Nx:1:-1)
DO irow=1,Nz
  icol1 = MAX(1,irow-Nx+1)
  icol2 = MIN(Ny,irow)
  ix1   = MAX(1,Nx-irow+1)
  ix2   = MIN(Nz-irow+1,Nx)
  XMAT(irow,icol1:icol2) = xr(ix1:ix2)
ENDDO

z = MATMUL(XMAT,y)
RETURN
END SUBROUTINE CONVT
!=======================================================================

SUBROUTINE DCONVT(z,x,y,Nz,Nx,Ny)
!!=======================================================================
!!
!! Time-domain deconvolution in matrix formulation:
!!    If z=x*y is the  Nz=Nx+Ny-1 length 
!!    convolution, compute y by deconvolving z 
!!    by x. Considering the convolution via
!!    matrix multiplication z=Xy, where X is an
!!    Nz by Ny matrix, the deconvolution is 
!!    carried out using Lapack DGELSS
!!
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: Nx,Ny,Nz
INTEGER(KIND=IB):: irow,icol1,icol2,ix1,ix2
REAL(KIND=RP)   :: x(Nx), y(Ny), z(Nz),xr(Nx)
REAL(KIND=RP)   :: XMAT(Nz,Ny)
!! Lapack variables:
INTEGER(KIND=IB)          :: MRANK, LWORK, INFO,LDA, LDB
INTEGER(KIND=IB),PARAMETER:: LWMAX = 4000, NRHS = 1

!! These are always double precision for DGELLS to work
REAL(KIND=DRP)          :: SV(Ny),WORK(LWMAX)
REAL(KIND=DRP)          :: XMAT2(Nz,Ny),b(NZ)
REAL(KIND=DRP),PARAMETER:: RCOND = 1.e-12_DRP

b = REAL(z,DRP)
LDA = Nz
LDB = Nz
!!
!! Build linear system of equations:
XMAT = 0._RP
xr = x(Nx:1:-1)
DO irow=1,Nz
  icol1 = MAX(1,irow-Nx+1)
  icol2 = MIN(Ny,irow)
  ix1   = MAX(1,Nx-irow+1)
  ix2   = MIN(Nz-irow+1,Nx)
  XMAT(irow,icol1:icol2) = xr(ix1:ix2)
ENDDO

!! In Matlab, can use Moore-Penrose pseudo inverse:
!! y = pinv(X)*z;
!! Here, use Lapack to solve LLS via SGELSS:
!! SGELSS computes the minimum norm solution to a real linear least
!! squares problem:
!!
!! Minimize 2-norm(| b - A*x |).
!!
!! using the singular value decomposition (SVD) of A. A is an M-by-N
!! matrix which may be rank-deficient.
!!
!! Several right hand side vectors b and solution vectors x can be
!! handled in a single call; they are stored as the columns of the
!! M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
!! X.
!!
!! The effective rank of A is determined by treating as zero those
!! singular values which are less than RCOND times the largest singular
!! value.
!! subroutine sgelss(integer M,integer N,integer NRHS,real, dimension( lda, * ) A,
!!                   integer LDA,real, dimension( ldb, * ) B,integer LDB,
!!                   real, dimension( * ) S,real RCOND,integer RANK,
!!                   real, dimension( * ) WORK,integer LWORK,integer INFO )

!! first, figure out optimal work length...
XMAT2 = REAL(XMAT,DRP)
LWORK = -1
CALL DGELSS(Nz, Ny, NRHS, XMAT2, LDA, b, LDB, SV, RCOND, MRANK, WORK, LWORK, INFO)
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!! Carry out DGELSS with optimal length
CALL DGELSS(Nz, Ny, NRHS, XMAT2, LDA, b, LDB, SV, RCOND, MRANK, WORK, LWORK, INFO)
IF(INFO /= 0)WRITE(*,*)'WARNING DGELSS unstable!'

y = REAL(b(1:Ny),RP)
RETURN
END SUBROUTINE DCONVT
!=======================================================================
! This is the end my fiend...
! EOF
