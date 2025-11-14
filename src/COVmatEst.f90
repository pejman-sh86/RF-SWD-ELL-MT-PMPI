SUBROUTINE COVmatEst(Cd, Cdi, res, logW, NSMP, NDAT, nfrac, MAX_NAVE, damp_power, inonstat, iunbiased, imr, ierr)

USE DATA_TYPE
IMPLICIT NONE

INTEGER(KIND=IB), INTENT(IN)                     :: NSMP, NDAT 
REAL(KIND=RP), DIMENSION(NSMP, NDAT), INTENT(IN) :: res
REAL(KIND=RP), DIMENSION(NSMP), INTENT(IN)       :: logW
INTEGER(KIND=IB), INTENT(IN)                     :: nfrac, MAX_NAVE, inonstat, iunbiased, imr
REAL(KIND=RP), INTENT(IN)                        :: damp_power
REAL(KIND=RP), DIMENSION(NDAT,NDAT), INTENT(OUT) :: Cd, Cdi
INTEGER(KIND=IB), INTENT(OUT) :: ierr 

REAL(KIND=RP), DIMENSION(NDAT,NDAT)              :: CdSMP
REAL(KIND=RP), DIMENSION(NSMP)                   :: WCd
INTEGER(KIND=IB)                                 :: ismp, istart, iend, idat, irow, icol
INTEGER(KIND=IB)                                 :: NAVE, NAVE1, NAVE2
REAL(KIND=RP), DIMENSION(NDAT)                   :: res_mr, sd_nonstat, res_ns, res_s, res_ns_mr
REAL(KIND=RP)                                    :: sd_stat
REAL(KIND=RP), DIMENSION(2*NDAT-1)               :: Acovn  !! normalized autocovariance
REAL(KIND=RP), DIMENSION(2*NDAT-1)               :: taper_Acov  !! normalized autocovariance
INTEGER(KIND=IB), DIMENSION(2*NDAT-1)            :: lag

REAL(KIND=RP) :: MEAN, STD

istart = 1
iend = NDAT

!!-------
!! Estimate non stationary standard deviation

NAVE = NINT( REAL(NDAT,RP)/REAL(nfrac,RP), IB)
IF (NAVE > MAX_NAVE) NAVE = MAX_NAVE

WCd = EXP( logW - MAXVAL(logW) )
Cd = 0._RP
DO ismp = 1, NSMP

IF (imr==1) THEN
    res_mr = res(ismp,:) - MEAN(res(ismp,:), NDAT) !! residuals mean removed
ELSE
    res_mr = res(ismp,:)
END IF
!WRITE(*,10) res_mr

sd_stat = STD(res_mr, NDAT, 1) 
!WRITE(*,10) sd_stat

IF (inonstat==1) THEN
    DO idat = istart, iend

        IF (idat <= istart+NAVE) THEN
            NAVE1 = istart
        ELSE
            NAVE1 = idat - NAVE
        END IF
        IF (idat >= iend-NAVE-1) THEN
            NAVE2 = iend
        ELSE
            NAVE2 = idat + NAVE
        END IF

        sd_nonstat(idat) = SQRT( MEAN(res_mr(NAVE1:NAVE2)**2,NAVE2-NAVE1+1)  )

    END DO
ELSE !! inonstat

    sd_nonstat = sd_stat

END IF !! inonstat
!WRITE(*,10) sd_nonstat
!!-----------

res_ns = res(ismp,:) / sd_nonstat   !! standardize residuals by non-stationary stds
res_s  = res(ismp,:) / sd_stat
!WRITE(*,10) res_ns    

!!------------
! compute autocovariance

res_ns_mr = res_ns - MEAN(res_ns, NDAT)
CALL XCORR( res_ns_mr, res_ns_mr, Acovn, lag, NDAT, NDAT, 2*NDAT-1, iunbiased ) 
!WRITE(*,10) Acovn     
!WRITE(*,*) lag       

!!-----------
! tapering to help positive definiteness

taper_Acov = COS( PI2*REAL(lag,RP)/(2*REAL(NDAT,RP)-1) ) ** damp_power 
Acovn = Acovn * taper_Acov

!!----------
!construct covariance matrix

DO irow = 1, NDAT
    CdSMP(irow, :) = Acovn( NDAT-(irow-1) : SIZE(Acovn)-(irow-1) )
    DO icol = 1, NDAT
        CdSMP(irow,icol) = CdSMP(irow,icol) * sd_nonstat(irow) * sd_nonstat(icol) 
    END DO
END DO

!DO irow = 1, NDAT
!    DO icol = 1, NDAT
!        CdSMP(irow,icol) = CdSMP(irow,icol) * sd_nonstat(irow) * sd_nonstat(icol) 
!    END DO
!END DO

Cd = Cd + WCd(ismp)*CdSMP 
END DO !!ismp 
!IF (NSMP>1_IB) Cd = Cd / REAL(NSMP,KIND=RP)
IF (NSMP>1_IB) Cd = Cd / SUM(WCd)

!!------------
!Calculate the inverse covariance matrix

Cdi = Cd  !! just copying the lower triangular part of Cd into Cdi is sufficient
!! cholesky factorization of Cd
CALL DPOTRF('L', NDAT, Cdi, NDAT, ierr)
IF (ierr/=0) RETURN
!! Inverse computation: stores lower triangular part of inverse covariance into Cdi
CALL DPOTRI('L', NDAT, Cdi, NDAT, ierr)  
IF (ierr/=0) RETURN

DO icol = 1, NDAT
    DO irow = 1, icol-1
        Cdi(irow,icol) = Cdi(icol,irow)
    END DO
END DO

10 FORMAT(20000ES16.7)

END SUBROUTINE COVmatEst

!!---------------------------------------------------------------------------------------------------------------

FUNCTION MEAN(A, N)

USE DATA_TYPE
IMPLICIT NONE

INTEGER(KIND=IB) :: N
REAL(KIND=RP), DIMENSION(N) :: A
REAL(KIND=RP) :: MEAN

MEAN = SUM(A) / REAL(N,RP)

RETURN
END FUNCTION MEAN
!!---------------------------

FUNCTION STD(A, N, w)

USE DATA_TYPE
IMPLICIT NONE


INTEGER(KIND=IB) :: N, w !! w=0: unbiased (divided by N-1), w=1: biased (divided by N)
REAL(KIND=RP), DIMENSION(N) :: A
REAL(KIND=RP) :: STD 

REAL(KIND=RP) :: MEAN
INTEGER(KIND=IB) :: denom

IF (w==0) THEN
    denom = N-1
ELSEIF (w==1) THEN
    denom = N
END IF

STD = SQRT( SUM( (A-MEAN(A,N))**2 ) / denom )

RETURN
END FUNCTION STD
!!------------------------------
