SUBROUTINE UpdateCOV(obj)

USE RJMCMC_COM
IMPLICIT NONE

TYPE(objstruc) :: obj   !!INTENT(INOUT)

COMPLEX(KIND=RP), DIMENSION(NsampleDres, NDAT_MT) :: Complex_Dres
REAL(KIND=RP) :: Ferro_norm, ZFerro_norm
INTEGER(KIND=IB) :: is, ie, is2, ie2, ierror                      


IF ( cov_converged_datasets(1)==0 ) THEN

    IF (ICOV==2) THEN
        Cd_old = Cd       !!!The best approach is to use pointers instead of copying arrays
        Cdi_old = Cdi
    END IF
    is = 1
    ie = NRF2*NTIME
    CALL COVmatEst(Cd, Cdi, sampleDres(:,is:ie), sampleDres(:,ncount3), NsampleDres, NTIME, nfrac_RV, MAX_NAVE_RV, damp_power_RV, inonstat_RV, iunbiased_RV, imr_RV, ierror)
    IF (ierror/=0) THEN
        WRITE(*,*) 'non positive definite RF cov estimate at iter: ',icovIter
        WRITE(*,*) 'algorithm continues with previous cov matrix'
        IF (ICOV==2) THEN
            Cd = Cd_old       !!!The best approach is to use pointers instead of copying arrays
            Cdi = Cdi_old
        END IF
    ELSE IF (ICOV==2) THEN
        IF (iconverge_criterion<=2) THEN
            IF (iconverge_criterion_RV==1) THEN
                covIter_errRV = MAXVAL(ABS(Cd-Cd_old))  / MAXVAL(ABS(Cd_old))
            ELSEIF (iconverge_criterion_RV==2) THEN            
                covIter_errRV =  MAXVAL(ABS(Cd-Cd_old)) 
            ELSEIF (iconverge_criterion_RV==3) THEN            
                covIter_errRV =  Ferro_norm(Cd-Cd_old,NTIME) / Ferro_norm(Cd_old,NTIME)
            ELSEIF (iconverge_criterion_RV==4) THEN            
                covIter_errRV =  Ferro_norm(Cd-Cd_old,NTIME)
            END IF
            IF (rank==src) WRITE(*,*) 'RF covariance convergence criterion value: ', covIter_errRV 
            IF (covIter_errRV < converge_threshold_RV) cov_converged_datasets(1)=1
        END IF
    ELSE !!ierror
        ICOV = 2
        ISD_RV = ISD_RV_covIter
        sdmn(1) = sdmn_covIter(1)
        sdmx(1) = sdmx_covIter(1)
        obj%sdparR = sdpar_covIter(1)
        minlimsd   = sdmn(1)
        maxlimsd   = sdmx(1)
        pertsdsdsc = 10._RP
        maxpertsd  = maxlimsd-minlimsd
        pertsdsd   = maxpertsd/pertsdsdsc
    END IF !!ierror

END IF

IF ( cov_converged_datasets(2)==0 ) THEN
    
    IF (ICOV_SWD==2) THEN
        CdSWD_old = CdSWD
        CdiSWD_old = CdiSWD
    END IF
    is = NRF2*NTIME + 1
    ie = NRF2*NTIME + NMODE2*NDAT_SWD
    CALL COVmatEst(CdSWD, CdiSWD, sampleDres(:,is:ie), sampleDres(:,ncount3), NsampleDres, NDAT_SWD, nfrac_SWD, MAX_NAVE_SWD, damp_power_SWD, inonstat_SWD, iunbiased_SWD, imr_SWD, ierror)
    IF (ierror/=0) THEN
        WRITE(*,*) 'non positive definite SWD cov estimate at iter: ',icovIter
        WRITE(*,*) 'algorithm continues with previous cov matrix'
        IF (ICOV_SWD==2) THEN
            CdSWD = CdSWD_old       !!!The best approach is to use pointers instead of copying arrays
            CdiSWD = CdiSWD_old
        END IF
    ELSE IF (ICOV_SWD==2) THEN
        IF (iconverge_criterion<=2) THEN
            IF (iconverge_criterion_SWD==1) THEN
                covIter_errSWD = MAXVAL(ABS(CdSWD-CdSWD_old))  / MAXVAL(ABS(CdSWD_old))
            ELSEIF (iconverge_criterion_SWD==2) THEN            
                covIter_errSWD =  MAXVAL(ABS(CdSWD-CdSWD_old)) 
            ELSEIF (iconverge_criterion_SWD==3) THEN            
                covIter_errSWD =  Ferro_norm(CdSWD-CdSWD_old,NDAT_SWD) / Ferro_norm(CdSWD_old,NDAT_SWD)
            ELSEIF (iconverge_criterion_SWD==4) THEN            
                covIter_errSWD =  Ferro_norm(CdSWD-CdSWD_old,NDAT_SWD)
            END IF
            IF (rank==src) WRITE(*,*) 'SWD covariance convergence criterion value: ', covIter_errSWD
            IF (covIter_errSWD < converge_threshold_SWD) cov_converged_datasets(2)=1
        END IF
    ELSE !!ierror
        ICOV_SWD = 2
        ISD_SWD = ISD_SWD_covIter
        sdmn(2) = sdmn_covIter(2)
        sdmx(2) = sdmx_covIter(2)
        obj%sdparSWD = sdpar_covIter(2) 
        minlimsdSWD   = sdmn(2)
        maxlimsdSWD   = sdmx(2)
        pertsdsdscSWD = 10._RP
        maxpertsdSWD  = maxlimsdSWD-minlimsdSWD
        pertsdsdSWD   = maxpertsdSWD/pertsdsdscSWD
    END IF  !!ierror

END IF

IF (  cov_converged_datasets(3)==0 ) THEN
    
    IF (ICOV_ELL==2) THEN
        CdELL_old = CdELL
        CdiELL_old = CdiELL
    END IF
    is = NRF2*NTIME + NMODE2*NDAT_SWD + 1
    ie = NRF2*NTIME + NMODE2*NDAT_SWD + NMODE_ELL2*NDAT_ELL
    CALL COVmatEst(CdELL, CdiELL, sampleDres(:,is:ie), sampleDres(:,ncount3), NsampleDres, NDAT_ELL, nfrac_ELL, MAX_NAVE_ELL, damp_power_ELL, inonstat_ELL, iunbiased_ELL, imr_ELL, ierror)
    IF (ierror/=0) THEN
        WRITE(*,*) 'non positive definite ELL cov estimate at iter: ',icovIter
        WRITE(*,*) 'algorithm continues with previous cov matrix'
    IF (ICOV_ELL==2) THEN
        CdELL = CdELL_old
        CdiELL = CdiELL_old
    END IF
    ELSE IF (ICOV_ELL==2) THEN
        IF (iconverge_criterion<=2) THEN
            IF (iconverge_criterion_ELL==1) THEN
                covIter_errELL = MAXVAL(ABS(CdELL-CdELL_old))  / MAXVAL(ABS(CdELL_old))
            ELSEIF (iconverge_criterion_ELL==2) THEN            
                covIter_errELL =  MAXVAL(ABS(CdELL-CdELL_old)) 
            ELSEIF (iconverge_criterion_ELL==3) THEN            
                covIter_errELL =  Ferro_norm(CdELL-CdELL_old,NDAT_ELL) / Ferro_norm(CdELL_old,NDAT_ELL)
            ELSEIF (iconverge_criterion_ELL==4) THEN            
                covIter_errELL =  Ferro_norm(CdELL-CdELL_old,NDAT_ELL)
            END IF
            IF (rank==src) WRITE(*,*) 'ELL covariance convergence criterion value: ', covIter_errELL
            IF (covIter_errELL < converge_threshold_ELL) cov_converged_datasets(3)=1
        END IF
    ELSE !!ierror
        ICOV_ELL = 2
        ISD_ELL = ISD_ELL_covIter
        sdmn(3) = sdmn_covIter(3)
        sdmx(3) = sdmx_covIter(3)
        obj%sdparELL = sdpar_covIter(3)
        minlimsdELL   = sdmn(3)
        maxlimsdELL   = sdmx(3)
        pertsdsdscELL = 10._RP
        maxpertsdELL  = maxlimsdELL-minlimsdELL
        pertsdsdELL   = maxpertsdELL/pertsdsdscELL
    END IF !!ierror

END IF

IF ( cov_converged_datasets(4)==0 ) THEN
    
    IF (ICOV_MT==2) THEN
        CdMT_old = CdMT
        CdiMT_old = CdiMT
    END IF
    is = NRF2*NTIME + NMODE2*NDAT_SWD + NMODE_ELL2*NDAT_ELL + 1
    ie = NRF2*NTIME + NMODE2*NDAT_SWD + NMODE_ELL2*NDAT_ELL + NMT2*NDAT_MT
    is2 = ie + 1_IB
    ie2 = ie + NMT2*NDAT_MT
    Complex_Dres = CMPLX(sampleDres(:,is:ie), sampleDres(:,is2:ie2), RP)
    CALL ZCOVmatEst(CdMT, CdiMT, Complex_Dres, sampleDres(:,ncount3), NsampleDres, NDAT_MT, nfrac_MT, MAX_NAVE_MT, damp_power_MT, inonstat_MT, iunbiased_MT, imr_MT, ierror)
    IF (ierror/=0) THEN
        WRITE(*,*) 'non positive definite MT cov estimate at iter: ',icovIter
        WRITE(*,*) 'algorithm continues with previous cov matrix'
        IF (ICOV_MT==2) THEN
            CdMT = CdMT_old       !!!The best approach is to use pointers instead of copying arrays
            CdiMT = CdiMT_old
        END IF
    ELSE IF (ICOV_MT==2) THEN
        IF (iconverge_criterion<=2) THEN
            IF (iconverge_criterion_MT==1) THEN
                covIter_errMT = MAXVAL( ABS(CdMT-CdMT_old) )  / MAXVAL( ABS(CdMT_old) )
            ELSEIF (iconverge_criterion_MT==2) THEN            
                covIter_errMT =  MAXVAL( ABS(CdMT-CdMT_old) ) 
            ELSEIF (iconverge_criterion_MT==3) THEN            
                covIter_errMT =  ZFerro_norm(CdMT-CdMT_old,NDAT_MT) / ZFerro_norm(CdMT_old,NDAT_MT)
            ELSEIF (iconverge_criterion_MT==4) THEN            
                covIter_errMT =  ZFerro_norm(CdMT-CdMT_old,NDAT_MT)
            END IF
            IF (rank==src) WRITE(*,*) 'MT covariance convergence criterion value: ', covIter_errMT
            IF (covIter_errMT < converge_threshold_MT) cov_converged_datasets(4)=1
        END IF
    ELSE !!ierror
        ICOV_MT = 2
        ISD_MT = ISD_MT_covIter
        sdmn(4) = sdmn_covIter(4)
        sdmx(4) = sdmx_covIter(4)
        obj%sdparMT = sdpar_covIter(4)
        minlimsdMT   = sdmn(4)
        maxlimsdMT   = sdmx(4)
        pertsdsdscMT = 10._RP
        maxpertsdMT  = maxlimsdMT-minlimsdMT
        pertsdsdMT   = maxpertsdMT/pertsdsdscMT
    END IF !!ierror

END IF

!! check for cov convergence
IF (iconverge_criterion<=2) THEN
   
    IF ( ALL( cov_converged_datasets == 1) ) THEN
            cov_converged = .TRUE.
            IF(rank==src) CALL SAVECOVS()
    ELSE 
        IF ( iconverge_criterion==1 ) THEN
            WHERE ( ICOViter_datasets==1 ) cov_converged_datasets=0
        END IF
    END IF

!ELSEIF (iconverge_criterion==3) THEN
   
    
!    covIter_err = ( Ferro_norm(Cd-Cd_old) + Ferro_norm(CdSWD-CdSWD_old) + Ferro_norm(CdELL-CdELL_old) + ZFerro_norm(CdMT-CdMT_old) ) &
 !               / ( Ferro_norm(Cd_old) + Ferro_norm(CdSWD_old) + Ferro_norm(CdELL_old) + ZFerro_norm(CdMT_old ) )
    
END IF  

END SUBROUTINE UpdateCOV
!!----------------------------

FUNCTION Ferro_norm(A, N)

USE DATA_TYPE
IMPLICIT NONE

INTEGER(KIND=IB) :: N
REAL(KIND=RP) :: A(N,N)
REAL(KIND=RP) :: Ferro_norm

Ferro_norm = SQRT( SUM(A**2) / REAL(N,KIND=RP) )

RETURN
END FUNCTION Ferro_norm
!!-----------------------------

FUNCTION ZFerro_norm(A, N)

USE DATA_TYPE
IMPLICIT NONE

INTEGER(KIND=IB) :: N
COMPLEX(KIND=RP) :: A(N,N)
REAL(KIND=RP) :: ZFerro_norm

ZFerro_norm = SQRT( SUM( ABS(A)**2 ) / REAL(N,KIND=RP) )

RETURN
END FUNCTION ZFerro_norm
!!-----------------------------

SUBROUTINE SAVECOVS()

USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: irow                      


IF (ICOV_iterUpdate_RV==1) THEN

    OPEN(UNIT=317, FILE=infileCd, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NTIME
        WRITE(317,100) Cd(irow,:)
    END DO 
    CLOSE(317)
    OPEN(UNIT=318, FILE=infileCdi, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NTIME
        WRITE(318,100) Cdi(irow,:)
    END DO 
    CLOSE(318)

END IF

IF (ICOV_iterUpdate_SWD==1) THEN

    OPEN(UNIT=417, FILE=infileCdSWD, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_SWD
        WRITE(417,100) CdSWD(irow,:)
    END DO 
    CLOSE(417)
    OPEN(UNIT=418, FILE=infileCdiSWD, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_SWD
        WRITE(418,100) CdiSWD(irow,:)
    END DO 
    CLOSE(418)

END IF

IF (ICOV_iterUpdate_ELL==1) THEN

    OPEN(UNIT=517, FILE=infileCdELL, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_ELL
        WRITE(517,100) CdELL(irow,:)
    END DO 
    CLOSE(517)
    OPEN(UNIT=518, FILE=infileCdiELL, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_ELL
        WRITE(518,100) CdiELL(irow,:)
    END DO 
    CLOSE(518)

END IF

IF (ICOV_iterUpdate_MT==1) THEN

    OPEN(UNIT=617, FILE=infileCdMT, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_MT
        WRITE(617,100) REAL(CdMT(irow,:),RP), AIMAG(CdMT(irow,:))
    END DO 
    CLOSE(617)
    OPEN(UNIT=618, FILE=infileCdiMT, STATUS='REPLACE', ACTION='WRITE')
    DO irow = 1, NDAT_MT
        WRITE(618,100) REAL(CdiMT(irow,:),RP), AIMAG(CdiMT(irow,:))
    END DO 
    CLOSE(618)

END IF

100 FORMAT(1000000ES20.10)

END SUBROUTINE SAVECOVS
