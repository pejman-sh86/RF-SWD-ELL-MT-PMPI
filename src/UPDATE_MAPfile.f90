SUBROUTINE UPDATE_MAPfile(sample_arr, nrow, ncol)

USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB), INTENT(IN) :: nrow, ncol
REAL(KIND=RP), DIMENSION(nrow, ncol), INTENT(IN) :: sample_arr

INTEGER(KIND=IB) :: ilayer!, icol
INTEGER(KIND=IB), DIMENSION(NLMX) :: layer_count != 0
INTEGER(KIND=IB) :: MAP_nlayer
!INTEGER(KIND=IB) :: MAP_nlayer_count
INTEGER(KIND=IB) :: MAP_idx
REAl(KIND=RP) :: MAP_logl
!INTEGER(KIND=IB) :: iend
REAL(KIND=RP), DIMENSION(ncol-9) :: MAP 

!WRITE(*,*) (sample_arr(ilayer,1), ilayer=1,5)
IF (iMAP_calc==0) THEN

  !DO ilayer = INT( MIN(sample_arr(:,4)), IB ), INT( MAX(sample_arr(:,4)) ) 
  DO ilayer = 1, NLMX
      layer_count(ilayer) = COUNT( sample_arr(:,4) == REAL(ilayer,RP) )  !!! In SAVESAMPLE subroutine, the layer number is converted to REAL type
  END DO   

  MAP_nlayer = MAXLOC(layer_count, 1, KIND=IB)
  !!If DIM (second argument) is absent, the result is a rank-one array with a length equal to the rank of ARRAY. If DIM is present, the result is an array with a rank one less than the rank of ARRAY, and a size corresponding to the size of ARRAY with the DIM dimension removed. If DIM is present and ARRAY has a rank of one, the result is a scalar. If the optional argument KIND is present, the result is an integer of kind KIND, otherwise it is of default kind.
  !MAP_nlayer_count = layer_count(MAP_nlayer)
  MAP_idx = MAXLOC( sample_arr(:,1), 1, sample_arr(:,4)==MAP_nlayer, KIND=IB )
  MAP_logl = sample_arr(MAP_idx, 1) 
  WRITE(*,*) 'MAP_logL =', MAP_logL

ELSE

  MAP_idx = MAXLOC( sample_arr(:,1), 1, KIND=IB )
  !MAP_logl = sample_arr(MAP_idx, 1)
  !MAP_nlayer = sample_arr(MAP_idx, 4)

END IF

!iend = SIZE(sample_arr, 2)
MAP = sample_arr(MAP_idx, 4:ncol-6) 
OPEN(UNIT=600, FILE=mapfile, STATUS='REPLACE', ACTION='WRITE')
WRITE(600,200) MAP 
CLOSE(600) 
200 FORMAT(20000ES16.7)

END SUBROUTINE UPDATE_MAPfile
