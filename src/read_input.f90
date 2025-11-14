!==============================================================================

SUBROUTINE READPARFILE()
!==============================================================================
!USE MPI
USE RJMCMC_COM
USE RAY3D_COM
IMPLICIT NONE
INTEGER(KIND=RP):: ntmp,iaz,ik,ilim
REAL(KIND=RP)   :: vref, dVs, VpVsmin, VpVsmax, dVpVs
INTERFACE
   FUNCTION LOGFACTORIAL(n)
     USE DATA_TYPE
     REAL(KIND=RP) :: LOGFACTORIAL
     REAL(KIND=RP),INTENT(IN):: n
   END FUNCTION LOGFACTORIAL
END INTERFACE

parfile        = filebase(1:filebaselen) // '_parameter.dat'!! Inversion parameter file
!!
!! Read parameter file
!!
OPEN(UNIT=20,FILE=parfile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) IMAP       !! 1 
READ(20,*) ICOV       !! 2
READ(20,*) IMAGSCALE  !! 3
READ(20,*) ENOS       !! 4 Even numbered order statistic on k (avoids thin layers)
READ(20,*) IPOIPR     !! 5 Applies Poisson prior on k
READ(20,*) IAR        !! 6
READ(20,*) I_VARPAR   !! 7
READ(20,*) IBD_SINGLE !! 8
READ(20,*) I_RV       !! 9Invert Radial and Vertical seismogram components or RF
READ(20,*) I_T        !! 10nvert Transverse component as well
READ(20,*) I_SWD      !!11Invert SWD data
READ(20,*) I_ELL      !!12Invert ELL data
READ(20,*) I_MT       !!13Invert MT data
READ(20,*) I_ZMT      !!14Invert MT complex impedance tensor
READ(20,*) I_VREF     !!15 
READ(20,*) I_VPVS     !!16 
READ(20,*) ISMPPRIOR  !!17 
READ(20,*) ISETSEED   !!18 
READ(20,*) IEXCHANGE  !!19 
READ(20,*) IDIP       !!20 
READ(20,*) NDAT_SWD   !!21 No. SWD data
READ(20,*) NMODE      !!22 No. SWD modes
READ(20,*) NDAT_ELL   !!23 No. ELL data
READ(20,*) NMODE_ELL  !!24 No. ELL modes
READ(20,*) NDAT_MT    !!25 No. MT data
READ(20,*) NTIME      !!26 No. time samples
READ(20,*) NSRC       !!27 No. time samples in source-time function
READ(20,*) NLMN       !!28 Max number of layers
READ(20,*) NLMX       !!29 Max number of layers
READ(20,*) ICHAINTHIN !!30 Chain thinning interval
READ(20,*) NKEEP      !!31 No. models to keep before writing
READ(20,*) NPTCHAINS1 !!32 Chain thinning interval
READ(20,*) dTlog      !!33 Temperature increment
READ(20,*) lambda     !!34 Lambda for Poisson prior on k
READ(20,*) hmx        !!35 Max. partition depth
READ(20,*) hmin       !!36 Min. layer thickness (must be small enough to not violate detailed balance)
READ(20,*) armxH      !!37 Max. AR prediction size
READ(20,*) armxV      !!38 Max. AR prediction size
READ(20,*) armxSWD    !!39 Max. AR prediction size
READ(20,*) armxELL    !!40 Max. AR prediction size
READ(20,*) TCHCKPT    !!41 Checkpointing interval in s
READ(20,*) shift2     !!42 time series shift for raysum
READ(20,*) width2     !!43 peak width for raysum (-1 returns impulse response)
READ(20,*) wl         !!44 water level for ray3d
READ(20,*) sampling_dt!!45 time series sampling rate raysum
READ(20,*) dVs        !!46 Vs one sided prior width (relative to background model)
READ(20,*) dVpVs      !!47 VpVs ratio one sided prior width
READ(20,*) sigmamin   !!48 ///
READ(20,*) sigmamax   !!49 ///
READ(20,*) sdmn       !!50 data (residual) error standard deviation prior lower limit
READ(20,*) sdmx       !!51 data (residual) error standard deviation prior upper limit
READ(20,*) iraysum    !!52 Use raysum (1) or ray3d (0)?
READ(20,*) mults      !!53 0 = no multiples; 2 = Josip's RF version (may include PmP depending on version of phaselist.f
READ(20,*) directs    !!54 0 = no multiples; 2 = Josip's RF version (may include PmP depending on version of phaselist.f !!!Pejman
READ(20,*) IDECON     !!55 0 = water level; 1 = ITD !!!Pejman
READ(20,*) tolerance  !!56 tol for ITD !!!Pejman
READ(20,*) maxbumps   !!57 for ITD !!!Pejman
READ(20,*) MAXG       !!58 for ITD!!!Pejman
READ(20,*) norm_ITD   !!59 for ITD!!!Pejman
READ(20,*) ISD_RV     !!60 1 = hierarchical sd
READ(20,*) ISD_SWD    !!61 1 = hierarchical sd
READ(20,*) ISD_ELL    !!62 1 = hierarchical sd
READ(20,*) ISD_MT     !!63 1 = hierarchical sd
READ(20,*) ITAPER     !!64 1 = taper pred data
READ(20,*) ICOV_SWD   !!65                     
READ(20,*) ICOV_ELL   !!66                      
READ(20,*) ICOV_MT    !!67                          
READ(20,*) ELL_verbose!!68 
READ(20,*) ELL_prec   !!69
READ(20,*) I_ABS_ELL  !!70 
READ(20,*) I_LOG10_ELL!!71 
READ(20,*) I_SAMPLING_TYPE_ELL !!72
READ(20,*) I_SET_STEP_ELL !!73
READ(20,*) STEP_SIZE_ELL  !!74
READ(20,*) I_SET_COUNT_ELL !!75
READ(20,*) COUNT_ELL  !!76 
READ(20,*) I_SET_RANGE_ELL !!77
!READ(20,*) VpVsmin   !! minimum VpVs ratio
!READ(20,*) VpVsmax   !! maximum VpVs ratio
CLOSE(20)

!! Allocate raysum oarameters
!! (these parameters are also used in ray3d)
CALL ALLOC_RAYSUM()
!! Read geometry file:
CALL readgeom(modname,baz2,slow2,sta_dx,sta_dy,ntr)
NRF1 = ntr
  
IF(iraysum /= 1)THEN
  slow2 = slow2 * 1000._RP
  baz2 = baz2 * 180._RP / PI2
ENDIF

IF (I_RV==-1) THEN
    NRF2 = NRF1
ELSE
    NRF2 = 0_IB
ENDIF

IF (I_SWD==1) THEN
    NMODE2 = NMODE
ELSE
    NMODE2 = 0_IB
ENDIF

IF (I_ELL==1) THEN
    NMODE_ELL2 = NMODE_ELL
ELSE
    NMODE_ELL2 = 0_IB
ENDIF

IF (I_MT==1) THEN
    NMT2 = 1_IB
ELSE
    NMT2 = 0_IB
ENDIF


kmin = NLMN
kmax = NLMX
!kmin = 1
!kmax = 8

!! Poisson Prior on k:
ALLOCATE(pk(kmax))
DO ik = kmin,kmax
  pk(ik)  = EXP(-lambda)*lambda**REAL(ik,RP)/EXP(LOGFACTORIAL(REAL(ik,RP)))
ENDDO

infileV        = filebase(1:filebaselen) // '_V_b.txt'
IF(I_RV >= 1)THEN
  infileR      = filebase(1:filebaselen) // '_R_b.txt'
ELSEIF(I_RV == -1)THEN
  infileR      = filebase(1:filebaselen) // '_RF.txt'
ENDIF
infileT        = filebase(1:filebaselen) // '_T_b.txt'
infileSWD      = filebase(1:filebaselen) // '_SWD.dat'
infile_sdSWD   = filebase(1:filebaselen) // '_sdSWD.dat'
infileELL      = filebase(1:filebaselen) // '_ELL.dat'
infileMT       = filebase(1:filebaselen) // '_MT.dat'
infileMT_ZVAR  = filebase(1:filebaselen) // '_MT_ZVAR.dat' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
infileCdi      = filebase(1:filebaselen) // '_Cdi.dat'
infileCdiSWD   = filebase(1:filebaselen) // '_CdiSWD.dat'
infileCdiELL   = filebase(1:filebaselen) // '_CdiELL.dat'
infileCdiMT    = filebase(1:filebaselen) // '_CdiMT.dat'
infileref      = filebase(1:filebaselen) // '_vel_ref.txt'
logfile        = filebase(1:filebaselen) // '_RJMH.log'
seedfile       = filebase(1:filebaselen) // '_seeds.log'
mapfile        = filebase(1:filebaselen) // '_map_voro.dat'
arfile         = filebase(1:filebaselen) // '_ar.dat'
predfile       = filebase(1:filebaselen) // '_mappred.dat'
obsfile        = filebase(1:filebaselen) // '_obs.dat'
arfileSWD      = filebase(1:filebaselen) // '_maparSWD.dat'
predfileSWD    = filebase(1:filebaselen) // '_mappredSWD.dat'
obsfileSWD     = filebase(1:filebaselen) // '_obsSWD.dat'
arfileELL      = filebase(1:filebaselen) // '_maparELL.dat'
predfileELL    = filebase(1:filebaselen) // '_mappredELL.dat'
obsfileELL     = filebase(1:filebaselen) // '_obsELL.dat'
predfileMT     = filebase(1:filebaselen) // '_mappredMT.dat'
obsfileMT      = filebase(1:filebaselen) // '_obsMT.dat'
!covfile        = filebase(1:filebaselen) // '_cov.txt'
sdfile         = filebase(1:filebaselen) // '_sigma.txt'
samplefile     = filebase(1:filebaselen) // '_voro_sample.txt'
stepsizefile     = filebase(1:filebaselen) // '_stepsize.txt'
IF(I_RV == 1 .OR. I_RV == -1 .OR. I_SWD == 1 .OR. I_ELL == 1) THEN
  IF(I_VPVS == 1) THEN
    NPL = 3
  ELSEIF(I_VPVS == -1) THEN
    NPL = 2
  END IF
  !!accounting for dipping layers
  IF(IDIP == 1) NPL=NPL+1      ! No. parameters per layer with dip
  IF(IDIP == 2) NPL=NPL+2      ! No. parameters per layer with strike and dip
ELSE
  NPL = 1
END IF


!! add MT
IF(I_MT == 1) NPL = NPL + 1

NTIME2 = NTIME+NSRC-1  ! Number time samples for zero padded observations
NRRG   = NTIME2 + NTIME - 1 ! No. points for convolution RRG
NRGRG  = NTIME  + NTIME - 1 ! No. points for convolution RGRG
NFPMX  = NLMX * NPL
NFPMX2 = NLMX * NPL * NPL

!! No. time samples
IF(rank == src) PRINT*,'NTIME =',NTIME,'NTIME2 =',NTIME2

ioutside = 0;ireject  = 0; iaccept = 0; iaccept_delay = 0; ireject_delay = 0
i_sdpert = 0;ishearfail = 0 ;i_ref_nlay = 0

!!
!!
!! Read velocity reference file
!!
201 FORMAT(a64)
202 FORMAT(a28,I4,I4)
203 FORMAT(4F12.3)
IF(I_VREF == 1)THEN
  OPEN(UNIT=20,FILE=infileref,FORM='formatted',STATUS='OLD',ACTION='READ')
  READ(20,*) NVELREF, ntmp
  NPREM = ntmp - NVELREF
  ALLOCATE(vel_ref(4,NVELREF),vel_prem(4,NPREM))
  IF(rank == src)WRITE(6,201) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
  IF(rank == src)WRITE(6,202) ' Velocity reference model:  ',NVELREF,NPREM
  IF(rank == src)WRITE(6,201) ' Reference:                                                     '
  IF(rank == src)WRITE(6,201) '   Depth(km)   Vs (km/s)      VpVs        Density               '
  !!  Read reference model to max sampling depth
  DO iaz = 1,NVELREF
    READ(20,*) vel_ref(:,iaz)
    IF(rank == src)WRITE(6,203) vel_ref(:,iaz)
  ENDDO
  !!  Read PREM beyond that
  IF(rank == src)WRITE(6,201) ' Prem (deep reference):                                         '
  DO iaz = 1,NPREM
    READ(20,*) vel_prem(:,iaz)
    IF(rank == src)WRITE(6,203) vel_prem(:,iaz)
  ENDDO
  CLOSE(20)
  IF(rank == src)WRITE(6,201) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
ENDIF

!!
!! Set width of Gaussian for raysum. -1. means get impulse response.
!!
IF(I_RV >= 0)THEN
  width2 = -1._SP
!! The width for the RF case is read from input parameter file
!ELSE
!  width = 2._SP
ENDIF

ALLOCATE( sdbuf(3,NRF1,NBUF) ) 
ALLOCATE(idxpar(NPL),sdevm((NLMX*NPL)+NPL-1,NLMX))
sdevm = 0._RP

ALLOCATE(minlim(NPL),maxlim(NPL),maxpert(NPL),pertsd(NPL),pertsdsc(NPL))
minlim   = 0._RP;maxlim   = 0._RP;maxpert  = 0._RP
pertsd   = 0._RP;pertsdsc = 30._RP

ALLOCATE(minlimar(3*NRF1),maxlimar(3*NRF1),maxpertar(3*NRF1))
ALLOCATE(pertarsd(3*NRF1),pertarsdsc(3*NRF1))
minlimar  = 0._RP;maxlimar  = 0._RP;maxpertar = 0._RP
pertarsd  = 0._RP;pertarsdsc= 18._RP

ALLOCATE(minlimarSWD(NMODE),maxlimarSWD(NMODE),maxpertarSWD(NMODE))
ALLOCATE(pertarsdSWD(NMODE),pertarsdscSWD(NMODE))
minlimarSWD  = 0._RP;maxlimarSWD  = 0._RP;maxpertarSWD = 0._RP
pertarsdSWD  = 0._RP;pertarsdscSWD= 18._RP    !! better to pertarSWDsd and pertarSWDsdsc, change them in the future

ALLOCATE(minlimarELL(NMODE_ELL),maxlimarELL(NMODE_ELL),maxpertarELL(NMODE_ELL))
ALLOCATE(pertarsdELL(NMODE_ELL),pertarsdscELL(NMODE_ELL))
minlimarELL  = 0._RP;maxlimarELL  = 0._RP;maxpertarELL = 0._RP
pertarsdELL  = 0._RP;pertarsdscELL= 18._RP

ALLOCATE(minlimsd(NRF1),maxlimsd(NRF1),maxpertsd(NRF1))
ALLOCATE(pertsdsd(NRF1),pertsdsdsc(NRF1))
minlimsd  = 0._RP;maxlimsd  = 0._RP;maxpertsd = 0._RP
pertsdsd  = 0._RP;pertsdsdsc= 18._RP

ALLOCATE(minlimsdSWD(NMODE),maxlimsdSWD(NMODE),maxpertsdSWD(NMODE))
ALLOCATE(pertsdsdSWD(NMODE),pertsdsdscSWD(NMODE))
minlimsdSWD  = 0._RP;maxlimsdSWD  = 0._RP;maxpertsdSWD = 0._RP
pertsdsdSWD  = 0._RP;pertsdsdscSWD= 18._RP    !! better to pertsdSWDsd and pertsdSWDsdsc. change them in the future

ALLOCATE(minlimsdELL(NMODE_ELL),maxlimsdELL(NMODE_ELL),maxpertsdELL(NMODE_ELL))
ALLOCATE(pertsdsdELL(NMODE_ELL),pertsdsdscELL(NMODE_ELL))
minlimsdELL  = 0._RP;maxlimsdELL  = 0._RP;maxpertsdELL = 0._RP
pertsdsdELL  = 0._RP;pertsdsdscELL= 18._RP

ALLOCATE(minlimsdMT(1),maxlimsdMT(1),maxpertsdMT(1))
ALLOCATE(pertsdsdMT(1),pertsdsdscMT(1))
minlimsdMT  = 0._RP;maxlimsdMT  = 0._RP;maxpertsdMT = 0._RP
pertsdsdMT  = 0._RP;pertsdsdscMT= 18._RP

IF(I_RV == -1 .OR. I_RV == 1 .OR. I_SWD == 1 .OR. I_ELL == 1)THEN
IF(I_VPVS == 1)THEN
  IF(iraysum == 1)THEN
    idxpar   = (/ 1, 4, 3 /)        !! Sample Vs
    IF(IDIP == 1)idxpar   = (/ 1, 4, 3 , 10 /)     !! Sample dip
    IF(IDIP == 2)idxpar   = (/ 1, 4, 3 ,  9, 10 /) !! Sample strike and dip
  ELSE
    !idxpar   = (/ 3, 4, 1 /)        !! Sample Vs
    !IF(IDIP == 1)idxpar   = (/ 3, 4, 1 , 7 /)    !! Sample dip
    !IF(IDIP == 2)idxpar   = (/ 3, 4, 1 , 6, 7 /) !! Sample strike and dip
    idxpar   = (/ 1, 4, 3 /)        !! Sample Vs
    IF(IDIP == 1)idxpar   = (/ 1, 4, 3, 7 /)    !! Sample dip
    IF(IDIP == 2)idxpar   = (/ 1, 4, 3, 6, 7 /) !! Sample strike and dip
  ENDIF
ELSEIF(I_VPVS == -1)THEN
  idxpar(1:2)   = (/ 1, 4 /)        !! Sample Vs
ELSE
  idxpar   = (/ 1, 3, 4 /)        !! Sample Vp
  IF(IDIP == 1)idxpar   = (/ 1, 3, 4 , 10 /)     !! Sample dip
  IF(IDIP == 2)idxpar   = (/ 1, 3, 4 ,  9, 10 /) !! Sample strike and dip
ENDIF
ENDIF
!!
!!  Prior bounds
!! (Note: Density is empirical through Birch's Law)
!!
!! Without dip:
!!           h     vs     vp/vs
ilim = 1
minlim(ilim) = hmin
maxlim(ilim) = hmx
ilim = ilim + 1
IF(I_RV == 1 .OR. I_RV == -1 .OR. I_SWD == 1 .OR. I_ELL == 1)THEN
  minlim(ilim) = -dVs
  maxlim(ilim) = dVs
  ilim = ilim + 1
  IF(I_VPVS == 1)THEN
    !! Sample Vs and VpVs ratio
    minlim(ilim) =  -dVpVs
    maxlim(ilim) =  dVpVs
    ilim = ilim + 1
  ENDIF
  IF(IDIP == 1)THEN
    minlim(ilim) = 1.0_RP
    maxlim(ilim) = 45.0_RP !! Frederiksen paper suggest method only stable for dip < 50 deg.
    ilim = ilim + 1
  ENDIF
  IF(IDIP == 2)THEN
    minlim(ilim) =    0.0_RP
    maxlim(ilim) =    100.0_RP
    ilim = ilim + 1
  ENDIF
END IF

IF(I_MT == 1)THEN
  minlim(ilim) = sigmamin 
  maxlim(ilim) = sigmamax
END IF

maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc

IF(IAR == 1)THEN
  !! Set prior and proposal scaling for AR model:
  minlimar   = -0.5000_RP
  maxlimar   =  0.90_RP
  pertarsdsc =  10._RP
  maxpertar  = maxlimar-minlimar
  pertarsd   = maxpertar/pertarsdsc
  minlimarSWD   = -0.5000_RP
  maxlimarSWD   =  0.90_RP
  pertarsdscSWD =  10._RP
  maxpertarSWD  = maxlimarSWD-minlimarSWD
  pertarsdSWD   = maxpertarSWD/pertarsdscSWD
  minlimarELL   = -0.5000_RP
  maxlimarELL   =  0.90_RP
  pertarsdscELL =  10._RP
  maxpertarELL  = maxlimarELL-minlimarELL
  pertarsdELL   = maxpertarELL/pertarsdscELL
ENDIF

IF(ICOV >= 1)THEN
  !! Set prior and proposal scaling for data error standard deviations:
  minlimsd   = sdmn(1)
  maxlimsd   = sdmx(1)
  pertsdsdsc = 10._RP  !! also it is set in UpdateCOV()
  maxpertsd  = maxlimsd-minlimsd
  pertsdsd   = maxpertsd/pertsdsdsc
ENDIF
IF(ICOV_SWD >= 1)THEN
  minlimsdSWD   = sdmn(2)
  maxlimsdSWD   = sdmx(2)
  pertsdsdscSWD = 10._RP  !! also it is set in UpdateCOV()
  maxpertsdSWD  = maxlimsdSWD-minlimsdSWD
  pertsdsdSWD   = maxpertsdSWD/pertsdsdscSWD
ENDIF
IF(ICOV_ELL >= 1)THEN
  minlimsdELL   = sdmn(3)
  maxlimsdELL   = sdmx(3)
  pertsdsdscELL = 10._RP  !! also it is set in UpdateCOV()
  maxpertsdELL  = maxlimsdELL-minlimsdELL
  pertsdsdELL   = maxpertsdELL/pertsdsdscELL
ENDIF
IF(ICOV_MT >= 1)THEN
  minlimsdMT   = sdmn(4)
  maxlimsdMT   = sdmx(4)
  pertsdsdscMT = 10._RP  !! also it is set in UpdateCOV()
  maxpertsdMT  = maxlimsdMT-minlimsdMT
  pertsdsdMT   = maxpertsdMT/pertsdsdscMT
ENDIF

END SUBROUTINE READPARFILE
!==============================================================================
!!
!! Write some info:
!!
!!=============================================================================

SUBROUTINE PRINTPAR2()

USE RJMCMC_COM
IMPLICIT NONE

!IF(rank == src)THEN
  WRITE(6,*) 'IMAP      = ', IMAP
  WRITE(6,*) 'ICOV      = ', ICOV
  WRITE(6,*) 'ICOV_SWD  = ', ICOV_SWD
  WRITE(6,*) 'ICOV_ELL  = ', ICOV_ELL
  WRITE(6,*) 'ICOV_MT   = ', ICOV_MT
  WRITE(6,*) 'I_RV      = ', I_RV
  WRITE(6,*) 'I_T       = ', I_T
  WRITE(6,*) 'I_SWD     = ', I_SWD
  WRITE(6,*) 'I_ELL     = ', I_ELL
  WRITE(6,*) 'I_MT      = ', I_MT
  WRITE(6,*) 'I_ZMT     = ', I_ZMT
  WRITE(6,*) 'IAR       = ', IAR
  WRITE(6,*) 'I_VPVS    = ', I_VPVS
  WRITE(6,*) 'ISMPPRIOR = ', ISMPPRIOR
  WRITE(6,*) 'ISETSEED  = ', ISETSEED
  WRITE(6,*) 'IEXCHANGE = ', IEXCHANGE
  WRITE(6,*) 'IDIP      = ', IDIP
  WRITE(6,*) 'NDAT_MT   = ', NDAT_MT
  WRITE(6,*) 'NTIME     = ', NTIME      !! No. time samples
  WRITE(6,*) 'NSRC      = ', NSRC       !! No. time samples in source-time function
  WRITE(6,*) 'NLMN      = ', NLMN       !! Max number of layers
  WRITE(6,*) 'NLMX      = ', NLMX       !! Max number of layers
  WRITE(6,*) 'ICHAINTHIN= ', ICHAINTHIN !! Chain thinning interval
  WRITE(6,*) 'NKEEP     = ', NKEEP      !! No. models to keep before writing
  IF (icovIter>0_IB) WRITE(6,*) 'NKEEP2    = ', NKEEP2      !! No. models to keep before writing
  WRITE(6,*) 'NPTCHAINS1= ', NPTCHAINS1 !! Chain thinning interval
  WRITE(6,*) 'dTlog     = ', dTlog      !! Temperature increment
  WRITE(6,*) 'hmx       = ', hmx        !! Max. partition depth
  WRITE(6,*) 'hmin      = ', hmin       !! Min. layer thickness (must be small enough to not violate detailed balance)
  WRITE(6,*) 'armxH     = ', armxH      !! Max. AR prediction size
  WRITE(6,*) 'armxV     = ', armxV      !! Max. AR prediction size
  WRITE(6,*) 'TCHCKPT   = ', TCHCKPT    !! Checkpointing interval in s
  WRITE(6,*) 'shift2    = ', shift2     !! time series shift for raysum
  WRITE(6,*) 'width2    = ', width2     !! time series shift for raysum
  WRITE(6,*) 'wter level= ', wl         !! water level for ray 3D
  WRITE(6,*) 'sampling_dt= ', sampling_dt!! time series sampling rate raysum
  WRITE(6,*) 'NRF1       = ', NRF1                         
  WRITE(6,*) 'ISD_RV     = ', ISD_RV                       
  WRITE(6,*) 'ISD_SWD    = ', ISD_SWD                      
  WRITE(6,*) 'ISD_ELL    = ', ISD_ELL                      
  WRITE(6,*) 'ISD_MT     = ', ISD_MT                       
  WRITE(6,*) 'Sample file: ',samplefile
  IF(ICOV_iterUpdate==1) WRITE(6,*) 'Sample file covIter   : ', samplefile_covIter
  WRITE(6,*) ''
  WRITE(6,*) 'minlim:  '
  WRITE(6,203) minlim
  WRITE(6,*) 'maxlim:  '
  WRITE(6,203) maxlim
  WRITE(6,*) 'pertsdsc:'
  WRITE(6,203) pertsdsc
  WRITE(6,*) 'minlim sigma(std):  '
  WRITE(6,203) minlimsd
  WRITE(6,*) 'maxlim sigma(std):  '
  WRITE(6,203) maxlimsd
  WRITE(6,*) 'pertsdsdsc:'
  WRITE(6,203) pertsdsdsc
  WRITE(6,*) 'minlim sigma(std) SWD:  '
  WRITE(6,203) minlimsdSWD
  WRITE(6,*) 'maxlim sigma(std) SWD:  '
  WRITE(6,203) maxlimsdSWD
  WRITE(6,*) 'pertsdsdscSWD:'
  WRITE(6,203) pertsdsdscSWD
  WRITE(6,*) 'minlim sigma(std) ELL:  '
  WRITE(6,203) minlimsdELL
  WRITE(6,*) 'maxlim sigma(std) ELL:  '
  WRITE(6,203) maxlimsdELL
  WRITE(6,*) 'pertsdsdscELL:'
  WRITE(6,203) pertsdsdscELL
  WRITE(6,*) 'minlim sigma(std) MT:  '
  WRITE(6,203) minlimsdMT
  WRITE(6,*) 'maxlim sigma(std) MT:  '
  WRITE(6,203) maxlimsdMT
  WRITE(6,*) 'pertsdsdscMT:'
  WRITE(6,203) pertsdsdscMT
  !WRITE(6,*) 'Done reading data.'
  IF (icovIter==0_IB) WRITE(6,*) 'Done reading parameter file.'
  WRITE(6,*) ''
  WRITE(6,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
  WRITE(6,*) ''
!ENDIF
CALL FLUSH(6)

203 FORMAT(100ES12.3)

RETURN
END SUBROUTINE PRINTPAR2   
!==============================================================================

SUBROUTINE READDATA(obj)
!==============================================================================
!USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=RP):: iaz,idat,io
TYPE(objstruc)  :: obj
REAL(KIND=RP), DIMENSION(NDAT_MT) :: CdiMT_R, CdiMT_I

IF(I_RV >= 1)THEN
  !! Radial and Vertical components:
  OPEN(UNIT=20,FILE=infileR,FORM='formatted',STATUS='OLD',ACTION='READ')
  OPEN(UNIT=21,FILE=infileV,FORM='formatted',STATUS='OLD',ACTION='READ')
  !!
  !! Observed are NTIME long but need to zero pad to NTIME2 for convolution
  !!
  obj%DobsR(1,1:NTIME2) = 0._RP
  obj%DobsV(1,1:NTIME2) = 0._RP
  DO iaz = 1,NRF1
    READ(20,*) obj%DobsR(iaz,1:NTIME)
    READ(21,*) obj%DobsV(iaz,1:NTIME)
  ENDDO
  CLOSE(20)
  CLOSE(21)
  IF(I_T == 1)THEN
    !! Transverse components:
    OPEN(UNIT=20,FILE=infileT,FORM='formatted',STATUS='OLD',ACTION='READ')
    obj%DobsT(1,1:NTIME2) = 0._RP
    DO iaz = 1,NRF1
      READ(20,*) obj%DobsT(iaz,1:NTIME)
    ENDDO
    CLOSE(20)
  ENDIF
ELSEIF(I_RV == -1)THEN
  !! RF case, read RF from file:
  OPEN(UNIT=20,FILE=infileR,FORM='formatted',STATUS='OLD',ACTION='READ')
  !!
  !! Observed RF are NTIME long
  !!
  obj%DobsR(1,:) = 0._RP
  DO iaz = 1,NRF1
    READ(20,*) obj%DobsR(iaz,1:NTIME)
    READ(20,*) taper_dpred(1:NTIME)
  ENDDO
  CLOSE(20)
  IF(ICOV == 2)THEN
    !! RF case, read Cd^-1 from file:
    OPEN(UNIT=30,FILE=infileCdi,FORM='formatted',STATUS='OLD',ACTION='READ')
    !!
    !! Cov is NTIME by NTIME
    !!
    !!ALLOCATE(Cdi(NTIME,NTIME))
    Cdi = 0._RP
    DO iaz = 1,NTIME
      READ(30,*) Cdi(iaz,1:NTIME)
    ENDDO
    CLOSE(30)
  ENDIF
ENDIF
IF(I_SWD == 1)THEN
  !!
  !! Surface wave dispersion data:
  !!
  OPEN(20,FILE=infileSWD,FORM='formatted',STATUS='OLD',ACTION='READ')
  !ndat = 0
  DO idat=1,NDAT_SWD
     READ(20,*,IOSTAT=io) obj%periods(1,idat),obj%DobsSWD(1,idat)
     IF (io > 0) THEN
       STOP "Check input.  Something was wrong"
     ELSEIF (io < 0) THEN
       EXIT
     ELSE
     !  ndatad=ndatad+1
     ENDIF
     !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
  ENDDO
  CLOSE(20)! close the file
  IF(ICOV_SWD == 2)THEN
    !! SWD case, read CdSWD^-1 from file:
    OPEN(UNIT=30,FILE=infileCdiSWD,FORM='formatted',STATUS='OLD',ACTION='READ')
    !!
    !! Cov is NDAT_SWD by NDAT_SWD
    !!
    !!ALLOCATE(CdiSWD(NDAT_SWD,NDAT_SWD))
    CdiSWD = 0._RP
    DO iaz = 1,NDAT_SWD
      READ(30,*) CdiSWD(iaz,1:NDAT_SWD)
    ENDDO
    CLOSE(30)
  ELSEIF (ICOV_SWD == 3) THEN
    ALLOCATE(sdSWD(NMODE,NDAT_SWD))
    sdSWD = 0._RP
    OPEN(20,FILE=infile_sdSWD,FORM='formatted',STATUS='OLD',ACTION='READ')
    DO idat=1,NDAT_SWD
       READ(20,*,IOSTAT=io) sdSWD(1,idat)
       IF (io > 0) THEN
         STOP "Check input.  Something was wrong"
       ELSEIF (io < 0) THEN
         EXIT
       ELSE
       !  ndatad=ndatad+1
       ENDIF
       !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
    ENDDO
    CLOSE(20)
  ENDIF!!ICOV_SWD
ENDIF
IF(I_ELL == 1)THEN
  !!
  !! Surface wave dispersion data:
  !!
  OPEN(20,FILE=infileELL,FORM='formatted',STATUS='OLD',ACTION='READ')
  !ndat = 0
  DO idat=1,NDAT_ELL
     READ(20,*,IOSTAT=io) obj%periods_ELL(1,idat),obj%DobsELL(1,idat)
     IF (io > 0) THEN
       STOP "Check input.  Something was wrong"
     ELSEIF (io < 0) THEN
       EXIT
     ELSE
     !  ndatad=ndatad+1
     ENDIF
     !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
  ENDDO
  CLOSE(20)! close the file
  IF(ICOV_ELL == 2)THEN
    !! ELL case, read CdELL^-1 from file:
    OPEN(UNIT=30,FILE=infileCdiELL,FORM='formatted',STATUS='OLD',ACTION='READ')
    !!
    !! Cov is NDAT_ELL by NDAT_ELL
    !!
    !!ALLOCATE(CdiELL(NDAT_ELL,NDAT_ELL))
    CdiELL = 0._RP
    DO iaz = 1,NDAT_ELL
      READ(30,*) CdiELL(iaz,1:NDAT_ELL)
    ENDDO
    CLOSE(30)
  ENDIF
ENDIF
IF(I_MT == 1)THEN
  !!
  !! MT data:
  !!
  OPEN(20,FILE=infileMT,FORM='formatted',STATUS='OLD',ACTION='READ')
  !READ(20,*,IOSTAT=io) NDAT_MT
  DO idat=1,NDAT_MT
     READ(20,*,IOSTAT=io) obj%freqMT(idat),obj%DobsMT(idat),obj%DobsMT(idat+NDAT_MT)
     IF (io > 0) THEN
       STOP "Check input.  Something was wrong"
     ELSEIF (io < 0) THEN
       EXIT
     ELSE
     !  ndatad=ndatad+1
     ENDIF
     !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
  ENDDO
  CLOSE(20)! close the file
  
  !!
  !! MT_ZVAR data:
  !!
  IF (I_ZMT == 1)THEN
  !IF (I_COV_MT == 1) THEN     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OPEN(20,FILE=infileMT_ZVAR,FORM='formatted',STATUS='OLD',ACTION='READ')
    !READ(20,*,IOSTAT=io) NDAT_MT
    ALLOCATE(MTZVAR(NDAT_MT)) !!!!!!!!!!!!!!!
    DO idat=1,NDAT_MT
       READ(20,*,IOSTAT=io) MTZVAR(idat)
       IF (io > 0) THEN
         STOP "Check input.  Something was wrong"
       ELSEIF (io < 0) THEN
         EXIT
       ELSE
       !  ndatad=ndatad+1
       ENDIF
       !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
    ENDDO
    CLOSE(20)
  !END IF !!I_COV_MT
  END IF !!I_ZMT

  IF(ICOV_MT == 2)THEN
    !! MT case, read CdMT^-1 from file:
    OPEN(UNIT=30,FILE=infileCdiMT,FORM='formatted',STATUS='OLD',ACTION='READ')
    !!
    !! Cov is NDAT_MT by NDAT_MT
    !!
    !!ALLOCATE(CdiMT(NDAT_MT,NDAT_MT))
    CdiMT = CMPLX(0._RP, 0._RP, RP)
    DO iaz = 1,NDAT_MT
      READ(30,*) CdiMT_R(1:NDAT_MT), CdiMT_I(1:NDAT_MT)
      CdiMT(iaz,1:NDAT_MT) = CMPLX(CdiMT_R(1:NDAT_MT), CdiMT_I(1:NDAT_MT), RP)
    ENDDO
    CLOSE(30)
  ENDIF
ENDIF !I_MT
RETURN
END SUBROUTINE READDATA
!==============================================================================

SUBROUTINE PRINTPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: ivo
TYPE(objstruc) :: obj

WRITE(6,*) 'Voronoi nodes:',obj%k
WRITE(6,*) 'NFP:',obj%NFP

DO ivo=1,obj%k
   WRITE(6,201) obj%voro(ivo,1:NPL)
   WRITE(6,205) obj%voroidx(ivo,1:NPL)
ENDDO
WRITE(6,*) 'Layer parameter vector:',obj%nunique,'layers'
DO ivo=1,obj%nunique
   WRITE(6,201) obj%par((ivo-1)*NPL+1:ivo*NPL)
ENDDO
WRITE(6,202) '            ',obj%par(obj%nunique*NPL+1:obj%nunique*NPL+(NPL-1))
WRITE(6,*) 'Partition:'
WRITE(6,203) obj%ziface(1:obj%nunique)
WRITE(6,*) 'Layers:'
WRITE(6,203) obj%hiface(1:obj%nunique)
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
IF(IAR == 1)THEN
   WRITE(6,*) 'AR parameters:'
   WRITE(6,206) 'R and V: ',obj%arpar
   WRITE(6,206) 'SWD:',obj%arparSWD
   WRITE(6,206) 'ELL:',obj%arparELL
ENDIF

201 FORMAT(6F12.4)
205 FORMAT(6I12)
202 FORMAT(A12,8F12.4)
203 FORMAT(11F12.4)
204 FORMAT(4F16.4)
!206 FORMAT(a,128F12.4)
206 FORMAT(a,128ES12.3)
END SUBROUTINE PRINTPAR
!!=======================================================================
RECURSIVE FUNCTION LOGFACTORIAL(n)  RESULT(fact)
!-----Factorial------------------------------------------------------
!!=======================================================================

USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) :: fact
REAL(KIND=RP), INTENT(IN) :: n

IF (n == 0) THEN
   fact = 0
ELSE
   fact = LOG(n) + LOGFACTORIAL(n-1)
END IF

END FUNCTION LOGFACTORIAL
!==============================================================================

SUBROUTINE READCOVPARFILE()
!!=============================================================================

USE RJMCMC_COM
IMPLICIT NONE


covparfile = filebase(1:filebaselen) // '_covparameter.dat'!! Inversion iterative cov parameter file
!!
!! Read parameter file
!!
OPEN(UNIT=20,FILE=covparfile,FORM='FORMATTED',STATUS='OLD',ACTION='READ')
READ(20,*) Icov_iterUpdate_RV
READ(20,*) Icov_iterUpdate_SWD
READ(20,*) Icov_iterUpdate_ELL
READ(20,*) Icov_iterUpdate_MT
READ(20,*) covIter_zero_nsamples
READ(20,*) covIter_period 
READ(20,*) MAXcovIter
READ(20,*) ICOVest    
READ(20,*) CHAINTHIN_COVest_period_zeroIter
READ(20,*) CHAINTHIN_COVest_period_nonzeroIter
READ(20,*) ISD_RV_covIter
READ(20,*) ISD_SWD_covIter
READ(20,*) ISD_ELL_covIter
READ(20,*) ISD_MT_covIter 
READ(20,*) sdmn_covIter 
READ(20,*) sdmx_covIter 
READ(20,*) sdpar_covIter 
READ(20,*) NKEEP_covIter 
READ(20,*) NKEEP_covIter_res 
READ(20,*) iSAVEsample_covIter 
READ(20,*) iSAVEsample_only_zeroIter 
READ(20,*) iMAP_calc 
READ(20,*) iconverge_criterion 
READ(20,*) iconverge_criterion_RV 
READ(20,*) iconverge_criterion_SWD 
READ(20,*) iconverge_criterion_ELL 
READ(20,*) iconverge_criterion_MT 
READ(20,*) converge_threshold_RV 
READ(20,*) converge_threshold_SWD 
READ(20,*) converge_threshold_ELL 
READ(20,*) converge_threshold_MT 
READ(20,*) nfrac_RV 
READ(20,*) MAX_NAVE_RV 
READ(20,*) inonstat_RV 
READ(20,*) iunbiased_RV 
READ(20,*) imr_RV 
READ(20,*) damp_power_RV 
READ(20,*) nfrac_SWD
READ(20,*) MAX_NAVE_SWD
READ(20,*) inonstat_SWD
READ(20,*) iunbiased_SWD
READ(20,*) imr_SWD
READ(20,*) damp_power_SWD
READ(20,*) nfrac_ELL
READ(20,*) MAX_NAVE_ELL
READ(20,*) inonstat_ELL
READ(20,*) iunbiased_ELL
READ(20,*) imr_ELL
READ(20,*) damp_power_ELL
READ(20,*) nfrac_MT 
READ(20,*) MAX_NAVE_MT 
READ(20,*) inonstat_MT 
READ(20,*) iunbiased_MT 
READ(20,*) imr_MT 
READ(20,*) damp_power_MT 
CLOSE(20) 

samplefile_covIter         = filebase(1:filebaselen) // '_voro_sample_covIter.txt'
samplefile_res_covIter     = filebase(1:filebaselen) // '_voro_sample_res_covIter.txt'
infileCd                   = filebase(1:filebaselen) // '_Cd.dat'
infileCdSWD                = filebase(1:filebaselen) // '_CdSWD.dat'
infileCdELL                = filebase(1:filebaselen) // '_CdELL.dat'
infileCdMT                 = filebase(1:filebaselen) // '_CdMT.dat'

!!
!! Write some info:
!!
IF(rank == src)THEN
  WRITE(6,*) 'Icov_iterUpdate_RV            = ', Icov_iterUpdate_RV
  WRITE(6,*) 'Icov_iterUpdate_SWD           = ', Icov_iterUpdate_SWD
  WRITE(6,*) 'Icov_iterUpdate_ELL           = ', Icov_iterUpdate_ELL
  WRITE(6,*) 'Icov_iterUpdate_MT            = ', Icov_iterUpdate_MT
  WRITE(6,*) 'covIter_zero_nsamples         = ', covIter_zero_nsamples
  WRITE(6,*) 'covIter_period                = ', covIter_period
  WRITE(6,*) 'MAXcovIter                    = ', MAXcovIter 
  WRITE(6,*) 'ICOVest                       = ', ICOVest
  WRITE(6,*) 'CHAINTHIN_COVest_period_zeroIter     = ', CHAINTHIN_COVest_period_zeroIter
  WRITE(6,*) 'CHAINTHIN_COVest_period_nonzeroIter  = ', CHAINTHIN_COVest_period_nonzeroIter
  WRITE(6,*) 'NKEEP_covIter                 = ', NKEEP_covIter
  WRITE(6,*) 'NKEEP_covIter_res             = ', NKEEP_covIter_res
  WRITE(6,*) 'iSAVEsample_covIter           = ', iSAVEsample_covIter
  WRITE(6,*) 'iSAVEsample_only_zeroIter     = ', iSAVEsample_only_zeroIter
  WRITE(6,*) 'iMAP_calc                     = ', iMAP_calc 
  WRITE(6,*) 'iconverge_criterion           = ', iconverge_criterion
  WRITE(6,*) 'iconverge_criterion_RV        = ', iconverge_criterion_RV 
  WRITE(6,*) 'iconverge_criterion_SWD       = ', iconverge_criterion_SWD
  WRITE(6,*) 'iconverge_criterion_ELL       = ', iconverge_criterion_ELL
  WRITE(6,*) 'iconverge_criterion_MT        = ', iconverge_criterion_MT
  WRITE(6,202) 'converge_threshold_RV          = ', converge_threshold_RV
  WRITE(6,202) 'converge_threshold_SWD         = ', converge_threshold_SWD
  WRITE(6,202) 'converge_threshold_ELL         = ', converge_threshold_ELL
  WRITE(6,202) 'converge_threshold_RV          = ', converge_threshold_MT
  WRITE(6,*) 'nfrac_RV                      = ', nfrac_RV
  WRITE(6,*) 'MAX_NAVE_RV                   = ', MAX_NAVE_RV
  WRITE(6,*) 'inonstat_RV                   = ', inonstat_RV
  WRITE(6,*) 'iunbiased_RV                  = ', iunbiased_RV
  WRITE(6,*) 'imr_RV                        = ', imr_RV   
  WRITE(6,202) 'damp_power_RV                  = ', damp_power_RV    
  WRITE(6,*) 'nfrac_SWD                     = ', nfrac_SWD
  WRITE(6,*) 'MAX_NAVE_SWD                  = ', MAX_NAVE_SWD
  WRITE(6,*) 'inonstat_SWD                  = ', inonstat_SWD
  WRITE(6,*) 'iunbiased_SWD                 = ', iunbiased_SWD
  WRITE(6,*) 'imr_SWD                       = ', imr_SWD  
  WRITE(6,202) 'damp_power_SWD                 = ', damp_power_SWD   
  WRITE(6,*) 'nfrac_ELL                     = ', nfrac_ELL
  WRITE(6,*) 'MAX_NAVE_ELL                  = ', MAX_NAVE_ELL
  WRITE(6,*) 'inonstat_ELL                  = ', inonstat_ELL
  WRITE(6,*) 'iunbiased_ELL                 = ', iunbiased_ELL
  WRITE(6,*) 'imr_ELL                       = ', imr_ELL  
  WRITE(6,202) 'damp_power_ELL                 = ', damp_power_ELL   
  WRITE(6,*) 'nfrac_MT                      = ', nfrac_MT
  WRITE(6,*) 'MAX_NAVE_MT                   = ', MAX_NAVE_MT
  WRITE(6,*) 'inonstat_MT                   = ', inonstat_MT
  WRITE(6,*) 'iunbiased_MT                  = ', iunbiased_MT
  WRITE(6,*) 'imr_MT                        = ', imr_MT   
  WRITE(6,202) 'damp_power_MT                  = ', damp_power_MT    
  WRITE(6,*) 'ISD_RV_covIter     = ', ISD_RV_covIter
  WRITE(6,*) 'ISD_SWD_covIter    = ', ISD_SWD_covIter
  WRITE(6,*) 'ISD_ELL_covIter    = ', ISD_ELL_covIter
  WRITE(6,*) 'ISD_MT_covIter     = ', ISD_MT_covIter
  WRITE(6,*) ''
  WRITE(6,*) 'min scaling factor:  '
  WRITE(6,203) sdmn_covIter
  WRITE(6,*) 'max scaling factor:  '
  WRITE(6,203) sdmx_covIter
  WRITE(6,*) 'starting scaling factor:  '
  WRITE(6,203) sdpar_covIter

  WRITE(6,*) 'Done reading covariance parameter file.'
  WRITE(6,*) ''
  WRITE(6,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
  WRITE(6,*) ''

END IF

202 FORMAT(a,100ES12.3)
203 FORMAT(100ES12.3)

END SUBROUTINE READCOVPARFILE
!!============================================================================

!EOF
