
SUBROUTINE MT1DforwardB3(ZMT, appRes, phase, sigma, h, angFre, Nl)

! This subroutine calculates vectors of 
! apparent resistivities and phases for 
! the range of input frequencies

USE data_type !module contains integer constants for sigle and double precision (SP and RP)
IMPLICIT NONE

! declaring inputs and outputs:

REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: sigma    !vector of layers conductivities
REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: h        !vector of layers thicknesses
REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: angFre   !vector of ANGULAR frequencies
INTEGER,                         INTENT(IN)  :: Nl

COMPLEX(KIND = RP), DIMENSION(:),   INTENT(OUT) :: ZMT      !vector of impedances
REAL(KIND = RP), DIMENSION(:),      INTENT(OUT) :: appRes   !vector of apparent resistivities
REAL(KIND = RP), DIMENSION(:),      INTENT(OUT) :: phase    !vector of phases in radians

! declaring local variables and constants:

!REAL(KIND = RP),    PARAMETER                 :: PI = DACOS(-1.0D0)          !DACOS is the specific function of ACOS for RP inputs
REAL(KIND = RP),    PARAMETER                 :: u0 = 4.0_RP * PI2 * 1.0D-7
REAL(KIND = RP),    DIMENSION(SIZE(angFre))   :: g                           !real wave number
COMPLEX(KIND = RP), DIMENSION(SIZE(angFre))   :: kstar                       !complex conjugate of yhe complex wave number
COMPLEX(KIND = RP), DIMENSION(SIZE(angFre))   :: R                           !correction factor for layered Earth
REAL(KIND=RP),      DIMENSION(SIZE(angFre))   :: real_R, imag_r              !real and imaginary parts of R
COMPLEX(KIND = RP), DIMENSION(SIZE(angFre))   :: x, a, b, c, d               !variables for holding intermediate results
REAL(KIND = RP),    DIMENSION(SIZE(angFre))   :: wu0                         !product of u0 and angular ferquency
INTEGER                                       :: I   !loop counter  

!PRINT*, 'initially inside the function '
!WRITE(*,*) 'size of angFre = ', SIZE(angFre)
!WRITE(*,*) 'size of R = ', SIZE(R)
!WRITE(*,*) 'size of appRes = ', SIZE(appRes)
!WRITE(*,*) 'size of phase = ', SIZE(phase)

R = (1.0_RP, 0.0_RP)    
!WRITE(*,*) R
wu0 = u0 * angFre

loop_on_layers: DO I = Nl-1, 1, -1

                   g = SQRT( wu0 * sigma(I) / 2.0_RP )
                   kstar = CMPLX(g, -g, RP)
                   x = SQRT( sigma(I)/sigma(I+1) ) * R
                   a = ( 1.0_RP - x ) / ( 1.0_RP + x )
                   b = a * EXP( -2.0_RP * kstar * h(I) )
                   R = ( 1.0_RP - b ) / ( 1.0_RP + b )

                END DO loop_on_layers  

!WRITE(*,*) R

c = SQRT( wu0 / (2.0_RP*sigma(1)) )
d = CMPLX(c, -c, RP)
ZMT = CONJG(d * R)
real_R = REAL(R)
imag_R = AIMAG(R)
!appRes = ( real_R**2 + imag_R**2 ) / sigma(1)
!phase = ATAN2( imag_R-real_R, imag_R+real_R ) * 180.0 / PI2               
appres = (ZMT*CONJG(ZMT))/wu0
phase = ATAN2( AIMAG(ZMT), REAL(ZMT) ) * 180.0 / PI2

END SUBROUTINE MT1DforwardB3
