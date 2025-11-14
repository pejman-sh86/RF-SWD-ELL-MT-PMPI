!=======================================================================

SUBROUTINE MAKE_MPI_STRUC_SP(obj,objtype)
!!=======================================================================
!!
!! 
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc) :: obj
INTEGER(KIND=IB) :: ifield,objtype
INTEGER(KIND=IB) :: oldtypes(NFIELD), blockcounts(NFIELD)
INTEGER(KIND=MPI_ADDRESS_KIND) :: offsets(NFIELD)
INTEGER(KIND=IB) :: iextent,rextent,dextent

!  Need to first figure offset by getting size of MPI_REAL etc 
call MPI_TYPE_EXTENT(MPI_INTEGER, iextent, ierr)
call MPI_TYPE_EXTENT(MPI_REAL, rextent, ierr)
call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dextent, ierr)

!!
!! Lengths of arrays
!!
!!                k   voro     voroidx           par          hiface     ziface      
blockcounts = (/ 1,  NLMX*NPL, NLMX*NPL , (NLMX+1)*NPL*NPL,  NPL*NLMX,   NPL*NLMX,&
!!               sdparR    sdparV   sdparT   sdparSWD  sdparELL   sdparMT    sdaveH    sdaveV   sdaveT   sdaveSWD   sdaveELL   sdaveMT
                 NRF1,      NRF1,    NRF1,    NMODE,  NMODE_ELL,      1,       NRF1,     NRF1,    NRF1,    NMODE,    NMODE_ELL,    1,&
!!               arpar   arparSWD arparELL  idxar   idxarSWD  idxarELL  gvoroidx nunique NFP
                 3*NRF1,  NMODE, NMODE_ELL, 3*NRF1,    NMODE,  NMODE_ELL,   NPL-1,     1,    1, & 
!!               beta    logL     logPr tcmp  ireject_bd  iaccept_bd ireject_bds  iaccept_bds 
                  1,       1,       1,    1,      1,          1,         1,          1,&
!!               DobsR           DpredR        DobsV       DpredV         DobsT        DpredT
                 NRF1*NTIME2,   NRF1*NTIME,  NRF1*NTIME2,  NRF1*NTIME,  NRF1*NTIME2,  NRF1*NTIME, &
!!                S            DresR         DresV        DresT        DarR        DarV        DarT
                 NRF1*NSRC, NRF1*NTIME,   NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME, &
!!                DobsSWD           DpredSWD        DresSWD       DarSWD           periods
                 NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD,&
!!                      DobsELL           DpredELL            DresELL             DarELL           periods_ELL
                 NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL, NMODE_ELL*NDAT_ELL,&
!!                    DobsMT       DpredMT      DresMT       freqMT 
                    2*NDAT_MT,    2*NDAT_MT,   2*NDAT_MT,    NDAT_MT /)

!! check that all MT adjustments are added

!  1   INTEGER(KIND=IB)                        :: k          ! No. nodes
!  2   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: voro       ! 1D voronoi nodes
!  3   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:):: voroidx    ! 1D voronoi nodes
!  4   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     ! Forward parameters
!  5   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: hiface
!  6   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: ziface
!  7   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparR        !! Std dev Radial comp
!  8   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparV        !! Std dev Vertical comp
!  9   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparT        !! Std dev Vertical comp
! 10   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparSWD      !! Std dev SWD data
! 11   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparELL      !! Std dev ELL data
! 12   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparMT       !! Std dev MT data
! 13   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveH        !! Std dev Radial comp running average for AR discrimination
! 14   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveV        !! Std dev Vertical comp running average for AR discrimination
! 15   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveT        !! Std dev Vertical comp running average for AR discrimination
! 16   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveSWD      !! Std dev SWD data running average for AR discrimination
! 17   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveELL      !! Std dev ELL data running average for AR discrimination
! 18   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveMT       !! Std dev MT data running average for AR discrimination
! 19   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arpar         !! AR model forward parameters
! 20   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparSWD      !! AR model forward parameters
! 21   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparELL      !! AR model forward parameters
! 22   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxar      !! AR on/off index (1=on)
! 23   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarSWD      !! AR on/off index (1=on)
! 24   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarELL      !! AR on/off index (1=on)
! 25   INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: gvoroidx   !! Index of live parameters on birth/death node
! 26   INTEGER(KIND=IB)                        :: nunique     !! 
! 27   INTEGER(KIND=IB)                        :: NFP         !! Number forward parameters
! 28   REAL(KIND=RP)                           :: beta
! 29   REAL(KIND=RP)                           :: logL        !! log likelihood
! 30   REAL(KIND=RP)                           :: logPr       !! log Prior probability ratio
! 31   REAL(KIND=RP)                           :: tcmp
! 32   INTEGER(KIND=IB)                        :: ireject_bd = 0
! 33   INTEGER(KIND=IB)                        :: iaccept_bd = 0
! 34   INTEGER(KIND=IB)                        :: ireject_bds = 0
! 35   INTEGER(KIND=IB)                        :: iaccept_bds = 0
! 36   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsR       !! Observed data for one logL eval.
! 37   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredR      !! Predicted data for trial model
! 38   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsV       !! Observed data for one logL eval.
! 39   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredV      !! Predicted data for trial model
! 40   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsT       !! Observed data for one logL eval.
! 41   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredT      !! Predicted data for trial model
! 42   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: S           !! Predicted data for trial model
! 43   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresR       !! Data residuals for trial model
! 44   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresV       !! Data residuals for trial model
! 45   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresT       !! Data residuals for trial model
! 46   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarR        !! Autoregressive model predicted data
! 47   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarV        !! Autoregressive model predicted data
! 48   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarT        !! Autoregressive model predicted data
! 49   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsSWD     !! Observed data SWD
! 50   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredSWD    !! Predicted SWD data for trial model
! 51   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresSWD     !! SWD data residuals for trial model
! 52   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarSWD      !! SWD autoregressive model predicted data
! 53   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: periods     !! SWD autoregressive model predicted data
! 54   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsELL     !! Observed data ELL
! 55   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredELL    !! Predicted ELL data for trial model
! 56   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresELL     !! ELL data residuals for trial model
! 57   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarELL      !! ELL autoregressive model predicted data
! 58   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: periods_ELL !! ELL autoregressive model predicted data
! 59   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: DobsMT      !! Observed data MT
! 60   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: DpredMT     !! Predicted MT data for trial model
! 61   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: DresMT      !! MT data residuals for trial model
! 62   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: freqMT      !! MT autoregressive model predicted data
!! Oldtypes
oldtypes(1)     = MPI_INTEGER
oldtypes(2)     = MPI_DOUBLE_PRECISION
oldtypes(3)     = MPI_INTEGER
oldtypes(4:21)   = MPI_DOUBLE_PRECISION
oldtypes(22:27) = MPI_INTEGER
oldtypes(28:31)    = MPI_DOUBLE_PRECISION
oldtypes(32:35) = MPI_INTEGER
oldtypes(36:62) = MPI_DOUBLE_PRECISION

call mpi_get_address(obj%k,offsets(1),ierr)
call mpi_get_address(obj%voro,offsets(2),ierr)
call mpi_get_address(obj%voroidx,offsets(3),ierr)
call mpi_get_address(obj%par,offsets(4),ierr)
call mpi_get_address(obj%hiface,offsets(5),ierr)
call mpi_get_address(obj%ziface,offsets(6),ierr)
call mpi_get_address(obj%sdparR,offsets(7),ierr)
call mpi_get_address(obj%sdparV,offsets(8),ierr)
call mpi_get_address(obj%sdparT,offsets(9),ierr)
call mpi_get_address(obj%sdparSWD,offsets(10),ierr)
call mpi_get_address(obj%sdparELL,offsets(11),ierr)
call mpi_get_address(obj%sdparMT,offsets(12),ierr)
call mpi_get_address(obj%sdaveH,offsets(13),ierr)
call mpi_get_address(obj%sdaveV,offsets(14),ierr)
call mpi_get_address(obj%sdaveT,offsets(15),ierr)
call mpi_get_address(obj%sdaveSWD,offsets(16),ierr)
call mpi_get_address(obj%sdaveELL,offsets(17),ierr)
call mpi_get_address(obj%sdaveMT,offsets(18),ierr)
call mpi_get_address(obj%arpar,offsets(19),ierr)
call mpi_get_address(obj%arparSWD,offsets(20),ierr)
call mpi_get_address(obj%arparELL,offsets(21),ierr)
call mpi_get_address(obj%idxar,offsets(22),ierr)
call mpi_get_address(obj%idxarSWD,offsets(23),ierr)
call mpi_get_address(obj%idxarELL,offsets(24),ierr)
call mpi_get_address(obj%gvoroidx,offsets(25),ierr)
call mpi_get_address(obj%nunique,offsets(26),ierr)
call mpi_get_address(obj%NFP,offsets(27),ierr)
call mpi_get_address(obj%beta,offsets(28),ierr)
call mpi_get_address(obj%logL,offsets(29),ierr)
call mpi_get_address(obj%logPr,offsets(30),ierr)
call mpi_get_address(obj%tcmp,offsets(31),ierr)
call mpi_get_address(obj%ireject_bd,offsets(32),ierr)
call mpi_get_address(obj%iaccept_bd,offsets(33),ierr)
call mpi_get_address(obj%ireject_bds,offsets(34),ierr)
call mpi_get_address(obj%iaccept_bds,offsets(35),ierr)
call mpi_get_address(obj%DobsR,offsets(36),ierr)
call mpi_get_address(obj%DpredR,offsets(37),ierr)
call mpi_get_address(obj%DobsV,offsets(38),ierr)
call mpi_get_address(obj%DpredV,offsets(39),ierr)
call mpi_get_address(obj%DobsT,offsets(40),ierr)
call mpi_get_address(obj%DpredT,offsets(41),ierr)
call mpi_get_address(obj%S,offsets(42),ierr)
call mpi_get_address(obj%DresR,offsets(43),ierr)
call mpi_get_address(obj%DresV,offsets(44),ierr)
call mpi_get_address(obj%DresT,offsets(45),ierr)
call mpi_get_address(obj%DarR,offsets(46),ierr)
call mpi_get_address(obj%DarV,offsets(47),ierr)
call mpi_get_address(obj%DarT,offsets(48),ierr)
!! SWD data
call mpi_get_address(obj%DobsSWD,offsets(49),ierr)
call mpi_get_address(obj%DpredSWD,offsets(50),ierr)
call mpi_get_address(obj%DresSWD,offsets(51),ierr)
call mpi_get_address(obj%DarSWD,offsets(52),ierr)
call mpi_get_address(obj%periods,offsets(53),ierr)
!! ELL data
call mpi_get_address(obj%DobsELL,offsets(54),ierr)
call mpi_get_address(obj%DpredELL,offsets(55),ierr)
call mpi_get_address(obj%DresELL,offsets(56),ierr)
call mpi_get_address(obj%DarELL,offsets(57),ierr)
call mpi_get_address(obj%periods_ELL,offsets(58),ierr)
!! MT data
call mpi_get_address(obj%DobsMT,offsets(59),ierr)
call mpi_get_address(obj%DpredMT,offsets(60),ierr)
call mpi_get_address(obj%DresMT,offsets(61),ierr)
call mpi_get_address(obj%freqMT,offsets(62),ierr)

DO ifield=2,SIZE(offsets)
  offsets(ifield) = offsets(ifield) - offsets(1)
ENDDO
offsets(1) = 0
!IF(rank == src)PRINT*,'offsets:',offsets

!  Now define structured type and commit it
call MPI_TYPE_CREATE_STRUCT( NFIELD, blockcounts, offsets, oldtypes, objtype,ierr)
call MPI_TYPE_COMMIT(objtype, ierr)

RETURN
END SUBROUTINE MAKE_MPI_STRUC_SP
!!==============================================================================

SUBROUTINE ALLOC_OBJ(obj)
!!==============================================================================
!!
!! Allocates memory.
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc) :: obj

ALLOCATE( obj%voro(NLMX,NPL),obj%voroidx(NLMX,NPL),obj%par((NLMX+1)*NPL*NPL) )
ALLOCATE( obj%hiface(NLMX*NPL),obj%ziface(NLMX*NPL) )
ALLOCATE( obj%sdparR(NRF1),obj%sdparV(NRF1),obj%sdparT(NRF1),obj%arpar(3*NRF1),obj%idxar(3*NRF1) )
ALLOCATE( obj%sdaveH(NRF1),obj%sdaveV(NRF1),obj%sdaveT(NRF1),obj%sdaveSWD(NMODE),obj%sdaveELL(NMODE_ELL) )
ALLOCATE( obj%sdparSWD(NMODE),obj%arparSWD(NMODE),obj%idxarSWD(NMODE) )
ALLOCATE( obj%sdparELL(NMODE_ELL),obj%arparELL(NMODE_ELL),obj%idxarELL(NMODE_ELL) )
ALLOCATE( obj%sdparMT(1),obj%sdaveMT(1) )
ALLOCATE( obj%gvoroidx(NPL-1) )
!! R, T, and V seismograms:
ALLOCATE( obj%DobsR(NRF1,NTIME2),obj%DpredR(NRF1,NTIME) )
ALLOCATE( obj%DobsV(NRF1,NTIME2),obj%DpredV(NRF1,NTIME) )
ALLOCATE( obj%DobsT(NRF1,NTIME2),obj%DpredT(NRF1,NTIME) )
ALLOCATE( obj%S(NRF1,NSRC) )
ALLOCATE( obj%DresR(NRF1,NTIME),obj%DresV(NRF1,NTIME),obj%DresT(NRF1,NTIME) )
ALLOCATE( obj%DarR(NRF1,NTIME),obj%DarT(NRF1,NTIME),obj%DarV(NRF1,NTIME) )
!! SWD data:
ALLOCATE( obj%DobsSWD(NMODE,NDAT_SWD),obj%DpredSWD(NMODE,NDAT_SWD),obj%DresSWD(NMODE,NDAT_SWD) )
ALLOCATE( obj%DarSWD(NMODE,NDAT_SWD),obj%periods(NMODE,NDAT_SWD) )
!! ELL data:
ALLOCATE( obj%DobsELL(NMODE_ELL,NDAT_ELL),obj%DpredELL(NMODE_ELL,NDAT_ELL),obj%DresELL(NMODE_ELL,NDAT_ELL) )
ALLOCATE( obj%DarELL(NMODE_ELL,NDAT_ELL),obj%periods_ELL(NMODE_ELL,NDAT_ELL) )
!! MT data:
ALLOCATE( obj%DobsMT(2*NDAT_MT),obj%DpredMT(2*NDAT_MT),obj%DresMT(2*NDAT_MT) )
ALLOCATE( obj%freqMT(NDAT_MT) )

obj%voro     = 0._RP
obj%voroidx  = 0
obj%par      = 0._RP
obj%hiface   = 0._RP
obj%ziface   = 0._RP
obj%sdparR   = 0._RP
obj%sdparV   = 0._RP
obj%sdparT   = 0._RP
obj%sdparSWD = 0._RP
obj%sdparELL = 0._RP
obj%sdparMT  = 0._RP
obj%sdaveH   = 0._RP
obj%sdaveV   = 0._RP
obj%sdaveT   = 0._RP
obj%sdaveSWD = 0._RP
obj%sdaveELL = 0._RP
obj%sdaveMT  = 0._RP
obj%arpar    = 0._RP
obj%idxar    = 0
obj%gvoroidx = 0
obj%DobsR    = 0._RP
obj%DpredR   = 0._RP
obj%DobsV    = 0._RP
obj%DpredV   = 0._RP
obj%DobsT    = 0._RP
obj%DpredT   = 0._RP
obj%S        = 0._RP
obj%DresR    = 0._RP
obj%DresV    = 0._RP
obj%DresT    = 0._RP
obj%DarR     = 0._RP
obj%DarV     = 0._RP
obj%DarT     = 0._RP
obj%DobsSWD  = 0._RP
obj%DpredSWD = 0._RP
obj%DresSWD  = 0._RP
obj%DarSWD   = 0._RP
obj%periods  = 0._RP
obj%DobsELL  = 0._RP
obj%DpredELL = 0._RP
obj%DresELL  = 0._RP
obj%DarELL   = 0._RP
obj%periods_ELL  = 0._RP
obj%DobsMT   = 0._RP
obj%DpredMT  = 0._RP
obj%DresMT   = 0._RP
obj%freqMT   = 0._RP
END SUBROUTINE ALLOC_OBJ
!=======================================================================

SUBROUTINE ALLOC_RAYSUM()
!!==============================================================================
!!
!! Reads observed data.
!!
USE RJMCMC_COM
IMPLICIT NONE
INCLUDE 'raysum/params.h'

ALLOCATE( baz2(maxtr),slow2(maxtr),sta_dx(maxtr),sta_dy(maxtr))
ALLOCATE( phaselist(maxseg,2,maxph),nseg(maxph),isoflag(maxlay) )
ALLOCATE( synth_cart(3,maxsamp,maxtr),synth_ph(3,maxsamp,maxtr) )
isoflag = .FALSE.

END SUBROUTINE ALLOC_RAYSUM
!!==============================================================================

SUBROUTINE ALLOC_COVmat()
!!==============================================================================

USE RJMCMC_COM
IMPLICIT NONE

IF (I_RV == -1) THEN

    IF (ICOV==2) ALLOCATE( Cdi(NTIME,NTIME) )
    IF (ICOV_iterUpdate_RV==1) THEN
        ALLOCATE( Cd(NTIME,NTIME) )
        ALLOCATE( Cd_old(NTIME,NTIME) )
        IF (.NOT.ALLOCATED(Cdi)) ALLOCATE( Cdi(NTIME,NTIME) )
        ALLOCATE( Cdi_old(NTIME,NTIME) )
    END IF

END IF !!I_RV

IF (I_SWD == 1) THEN

    IF (ICOV_SWD==2) ALLOCATE( CdiSWD(NDAT_SWD,NDAT_SWD) )
    IF (ICOV_iterUpdate_SWD==1) THEN
        ALLOCATE( CdSWD(NDAT_SWD,NDAT_SWD) )
        ALLOCATE( CdSWD_old(NDAT_SWD,NDAT_SWD) )
        IF (.NOT.ALLOCATED(CdiSWD)) ALLOCATE( CdiSWD(NDAT_SWD,NDAT_SWD) )
        ALLOCATE( CdiSWD_old(NDAT_SWD,NDAT_SWD) )
    END IF

END IF !!I_SWD

IF (I_ELL == 1) THEN

    IF (ICOV_ELL==2) ALLOCATE( CdiELL(NDAT_ELL,NDAT_ELL) )
    IF (ICOV_iterUpdate_ELL==1) THEN
        ALLOCATE( CdELL(NDAT_ELL,NDAT_ELL) )
        ALLOCATE( CdELL_old(NDAT_ELL,NDAT_ELL) )
        IF (.NOT.ALLOCATED(CdiELL)) ALLOCATE( CdiELL(NDAT_ELL,NDAT_ELL) )
        ALLOCATE( CdiELL_old(NDAT_ELL,NDAT_ELL) )
    END IF

END IF !!I_ELL

IF (I_MT == 1) THEN

    IF (ICOV_MT==2) ALLOCATE( CdiMT(NDAT_MT,NDAT_MT) )
    IF (ICOV_iterUpdate_MT==1) THEN
        ALLOCATE( CdMT(NDAT_MT,NDAT_MT) )
        ALLOCATE( CdMT_old(NDAT_MT,NDAT_MT) )
        IF (.NOT.ALLOCATED(CdiMT)) ALLOCATE( CdiMT(NDAT_MT,NDAT_MT) )
        ALLOCATE( CdiMT_old(NDAT_MT,NDAT_MT) )
    END IF

END IF !!I_MT

END SUBROUTINE ALLOC_COVmat
!!=============================================================================
!EOF
