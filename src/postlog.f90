PROGRAM POSTLOG 

USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE

CHARACTER(LEN=100)    :: infilesample, postlogfile
TYPE (objstruc)       :: obj
REAL(KIND=RP), ALLOCATABLE :: tmppar(:), Dpred(:)
!INTEGER(KIND=IB)      :: NRF2, NMODE2, NMODE_ELL2
INTEGER(KIND=IB)      :: NSMP, ismp, stat
INTEGER(KIND=IB)      :: ivo, ipar, ipar2, irow              


OPEN(UNIT=20,FILE=filebasefile,STATUS='OLD',ACTION='READ')
READ(20,*) filebaselen
READ(20,*) filebase
CLOSE(20)

!PRINT*, 'Enter the name of the station(filebase): '
!READ(*,*) filebase
!filebaselen = LEN(filebase)
  
!NSMP = 100000
infilesample = filebase(1:filebaselen) // '_' // 'sample_postpred.dat'
postlogfile = filebase(1:filebaselen) // '_' // 'postlog.dat'
OPEN(UNIT=193,FILE=postlogfile,STATUS='REPLACE',ACTION='WRITE')

CALL READPARFILE()

!I_RV = 0
ICOV = 2
!I_SWD = 1
!ICOV_SWD = 2
!I_ELL = 0
!I_MT = 0
!I_ZMT = 0
!ICOV_MT = 2

IMAP = 0
ICOV_iterUpdate_RV = 0
ICOV_iterUpdate_SWD = 0
ICOV_iterUpdate_ELL = 0
ICOV_iterUpdate_MT = 0

ncount1 = NFPMX+NAP+3*NRF1+3*NRF1+2*NMODE+2*NMODE_ELL+1
ALLOCATE(tmppar(ncount1))

IF (I_RV == -1) THEN
    NRF2 = NRF1
ELSE
    NRF2 = 0
END IF
IF (I_SWD == 1) THEN
    NMODE2 = NMODE
ELSE
    NMODE2 = 0
END IF
IF (I_ELL == 1) THEN
    NMODE_ELL2 = NMODE_ELL
ELSE
    NMODE_ELL2 = 0
END IF
ALLOCATE( Dpred(NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+2*NDAT_MT)  )

CALL ALLOC_OBJ(obj)
ALLOCATE( taper_dpred(NTIME) )
CALL ALLOC_COVmat()
CALL READDATA(obj)

NSMP = 0_IB
ismp = 1_IB
OPEN(UNIT=293,FILE=infilesample,STATUS='OLD',ACTION='READ')
DO
    READ(293,*, IOSTAT=stat) tmppar
    IF (stat>0) THEN
        WRITE(*,*) 'Something went wrong reading ', ismp, 'th sample'
    ELSEIF (stat<0) THEN
        WRITE(*,*) 'Reached end of file'
        EXIT
    ELSE
        IF (MOD(ismp,100)==0) PRINT*, 'Starting reading sample NO.: ', ismp 
        obj%k       = INT(tmppar(4),IB)
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
          ipar = (ivo-1)*NPL+5
          obj%voro(ivo,:) = tmppar(ipar:ipar+NPL-1)
          obj%voroidx(ivo,:) = 1
          DO ipar2 = 1,NPL
            IF(obj%voro(ivo,ipar2) < -99._RP)THEN
              obj%voroidx(ivo,ipar2) = 0
            ELSE 
              obj%voroidx(ivo,ipar2) = 1
            ENDIF
          ENDDO
        ENDDO  !!ivo

        IF(ICOV >= 1)THEN
        IF(I_RV == -1 .OR. I_RV == 1) obj%sdparR = tmppar(NFPMX+4+1:NFPMX+4+NRF1)
        IF(I_RV == 1) obj%sdparV   = tmppar(NFPMX+4+1+NRF1:NFPMX+4+2*NRF1)
        IF(I_T == 1) obj%sdparT   = tmppar(NFPMX+4+1+2*NRF1:NFPMX+4+3*NRF1)
        END IF
        IF(ICOV_SWD >= 1) THEN 
        IF(I_SWD == 1) obj%sdparSWD = tmppar(NFPMX+4+1+3*NRF1:NFPMX+4+3*NRF1+NMODE)
        END IF
        IF(ICOV_ELL >= 1) THEN 
          IF(I_ELL == 1) obj%sdparELL = tmppar(NFPMX+4+3*NRF1+NMODE+1:NFPMX+4+3*NRF1+NMODE+NMODE_ELL)
        END IF
        IF(ICOV_MT >= 1) THEN 
          IF(I_MT == 1) obj%sdparMT = tmppar(NFPMX+4+3*NRF1+NMODE+NMODE_ELL+1)
        END IF

        IF(I_VARPAR == 0)CALL INTERPLAYER_novar(obj)
        IF(I_VARPAR == 1)CALL INTERPLAYER(obj)
        IF (MOD(ismp,100)==0) PRINT*, 'Done reading sample NO.: ', ismp 
        
        IF (MOD(ismp,100)==0) PRINT*, 'Starting forward calculation for sample NO.: ', ismp 
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)
        IF (MOD(ismp,100)==0) PRINT*, 'Done forward calculation for sample NO.: ', ismp 

        !IF (I_RV==-1) Dpred(1:NRF2*NTIME) = (/ ( obj%DpredR(irow,1:NTIME), irow=1,NRF1 ) /)
        !IF (I_SWD==1) Dpred(NRF2*NTIME+1:NRF2*NTIME+NMODE2*NDAT_SWD) = (/ ( obj%DpredSWD(irow,1:NDAT_SWD), irow=1,NMODE ) /)
        !IF (I_ELL==1) Dpred(NRF2*NTIME+NMODE2*NDAT_SWD+1:NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL) = (/ ( obj%DpredELL(irow,1:NDAT_ELL), irow=1,NMODE_ELL ) /)
        !IF (I_MT==1) Dpred(NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+1:NRF2*NTIME+NMODE2*NDAT_SWD+NMODE_ELL2*NDAT_ELL+2*NDAT_MT) = obj%DpredMT
 
        !WRITE(193,208) obj%logL, NDAT_MT*LOG(obj%sdparMT)
        !WRITE(193,208) obj%logL, NDAT_SWD*LOG(obj%sdparSWD) / 2._RP
        WRITE(193,208) obj%logL, NTIME*LOG(obj%sdparR) / 2._RP
        IF (MOD(ismp,100)==0) PRINT*, 'Done writing logL for sample NO.: ', ismp 
        IF (MOD(ismp,100)==0) PRINT*, ''
    END IF !! IOSTAT

    IF(stat==0) NSMP = NSMP + 1_IB
    IF(stat>=0) ismp = ismp + 1_IB
END DO
CLOSE(293)
CLOSE(193)
PRINT*, ''
WRITE(*,*) 'number of samples in the input sample file: ', ismp-1
WRITE(*,*) 'number of samples in the output predictive file: ', NSMP
208 FORMAT(5000ES20.10)

END PROGRAM POSTLOG 
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
