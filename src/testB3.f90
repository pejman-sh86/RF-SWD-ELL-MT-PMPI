PROGRAM testB3

USE DATA_TYPE
IMPLICIT NONE

INTERFACE
    SUBROUTINE MT1DforwardB3(appRes, phase, sigma, h, angFre, Nl)
        USE data_type !module contains integer constants for sigle and double precision (SP and RP)
        IMPLICIT NONE

       REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: sigma    !vector of layers conductivities
       REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: h        !vector of layers thicknesses
       REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: angFre   !vector of ANGULAR frequencies
       INTEGER,                         INTENT(IN)  :: Nl

       REAL(KIND = RP), DIMENSION(:),   INTENT(OUT) :: appRes   !vector of apparent resistivities
       REAL(KIND = RP), DIMENSION(:),   INTENT(OUT) :: phase    !vector of phases in radians
    END SUBROUTINE MT1DforwardB3

END INTERFACE


INTEGER, PARAMETER               :: Nl = 3, Nf = 9
INTEGER                          :: i
REAL(KIND = RP), DIMENSION(Nl)   :: sigma = (/ 1.0D-1, 1.0D-4, 1.0D14 /)
REAL(KIND = RP), DIMENSION(Nl)   :: h = (/ 1.0, 49.0, 0.0 /) !in km
REAL(KIND = RP), DIMENSION(Nf)   :: fre = (/ (10.0D0**i, i = -Nf/2, Nf/2) /)
REAL(KIND = RP), DIMENSION(Nf)   :: appReS
REAL(KIND = RP), DIMENSION(Nf)   :: phase

REAL(KIND = RP), DIMENSION(2*NF) :: Dobs, Dpred, Dres
REAL(KIND = RP)                  :: std = 1.0D-4
REAL(KIND = RP)                  :: logL
!REAL(KIND = RP), PARAMETER       :: PI = DACOS(-1.0D0)
!REAL                             :: SECNDS, startTime, endTime

Dobs = 0.0_RP
!startTime = SECNDS(0.0)
!WRITE(*,*) 'fre = ', fre
CALL MT1DforwardB3(appRes, phase, sigma, 1000.0_RP*h, 2.0_RP * PI2 * fre, Nl)
!endTime = SECNDS(0.0)
Dpred(1:Nf) = appRes
Dpred(Nf+1:2*Nf) = phase
Dres = Dpred - Dobs
!WRITE(*, 1) endTime - startTime
!1 FORMAT ( 1X, 'elapsed time for B3 is : ', ES16.4, 's' )

!IF (endTime - startTime == 0.0) WRITE(*,*) 'yes'

WRITE(*, *) appRes
Write(*, *) phase
logL = LOG(1._RP/(2._RP*PI2)**(REAL(2*Nf,RP)/2._RP)) &
            -(SUM(Dres**2._RP)/(2._RP*std**2._RP)&
            +REAL(2*Nf,RP)*LOG(std))
WRITE(*,*) 'logL = ', logL

END PROGRAM testB3
