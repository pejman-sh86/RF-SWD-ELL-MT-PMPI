SUBROUTINE XCORR(X, Y, Z, lag, Nx, Ny, Nz, inormal)

USE DATA_TYPE
IMPLICIT NONE
!!inormal: 0 biased, 1 unbiased, 2 raw cross correlation

INTEGER(KIND=IB), INTENT(IN) :: Nx, Ny
REAL(KIND=RP), DIMENSION(Nx), INTENT(IN) :: X
REAL(KIND=RP), DIMENSION(Ny), INTENT(IN) :: Y
INTEGER(KIND=IB), INTENT(IN) :: Nz
INTEGER(KIND=IB), INTENT(IN) :: inormal  
REAL(KIND=RP), DIMENSION(Nz), INTENT(OUT) :: Z
INTEGER(KIND=IB), DIMENSION(Nz), INTENT(OUT) :: lag

REAL(KIND=RP), DIMENSION(Nz, Nx) :: YMAT
INTEGER(KIND=IB) :: irow, icol1, icol2, istart, iend

YMAT = 0._RP
DO irow = 1, Nz

    icol1 = MAX(1, Nx-irow+1)
    icol2 = MIN(Nx, Nz-irow+1)
    istart = MAX(1, irow-Nx+1)
    iend = MIN(irow, Ny)
    YMAT(irow, icol1:icol2) = Y(istart:iend)
    IF(Nx==Ny) THEN
        IF ( inormal==0 ) THEN
            YMAT(irow, icol1:icol2) = YMAT(irow, icol1:icol2) / REAL(Nx,RP) 
        ELSEIF ( inormal==1 ) THEN
            YMAT(irow, icol1:icol2) = YMAT(irow, icol1:icol2) / REAL(Nx-ABS(irow-Nx),RP)
        END IF
    END IF
    lag(irow) = irow - Nx

END DO

Z = MATMUL(YMAT, X)

END SUBROUTINE XCORR
