!==============================================================================
MODULE RJMCMC_COM
   USE MPI
   USE DATA_TYPE
   IMPLICIT NONE

!!
!! General switches
!!
   INTEGER(KIND=IB) :: IMAP       !! WRITE REPLICA AND EXIT
   INTEGER(KIND=IB) :: ICOV       !! 0 = Sample implicit over sigma
                                   !! 1 = Sample over sigma
   INTEGER(KIND=IB) :: ICOV_SWD   !! 0 = Sample implicit over sigma
   INTEGER(KIND=IB) :: ICOV_ELL   !! 0 = Sample implicit over sigma
   INTEGER(KIND=IB) :: ICOV_MT    !! 0 = Sample implicit over sigma
   INTEGER(KIND=IB) :: IMAGSCALE  !! 0 = turn off magnitude scale error
   INTEGER(KIND=IB) :: ENOS       !! 1 = Turn on even numbered order stats
   INTEGER(KIND=IB) :: IPOIPR     !! 1 = Turn on Poisson prior on k
   INTEGER(KIND=IB) :: IAR        !! 1 = Use Autoregressive error model
   INTEGER(KIND=IB) :: I_VARPAR   !! 1 = invert Radial and Vertical components
   INTEGER(KIND=IB) :: IBD_SINGLE !! 1 = include BD for single parameters onto nodes
   INTEGER(KIND=IB) :: I_RV       !! 1 = invert Radial and Vertical components
   INTEGER(KIND=IB) :: I_T        !! 1 = invert Transverse component as well
   INTEGER(KIND=IB) :: I_SWD      !! 1 = invert SWD data
   INTEGER(KIND=IB) :: I_ELL      !! 1 = invert ELL data
   INTEGER(KIND=IB) :: I_MT       !! 1 = invert MT data
   INTEGER(KIND=IB) :: I_ZMT      !! 1 = invert MT complex impedance tensor
   INTEGER(KIND=IB) :: iraysum    !! 1 = invert use raysum 0 = use ray3d
   INTEGER(KIND=IB) :: I_VREF     !! 1 = sample VpVs ratio
   INTEGER(KIND=IB) :: I_VPVS     !! 1 = sample VpVs ratio
   INTEGER(KIND=IB) :: ISMPPRIOR  !! Sample from the prior (sets logL = 1._RP)
   INTEGER(KIND=IB) :: ISETSEED   !! Fix the random seed 
   INTEGER(KIND=IB) :: IEXCHANGE  !! 1 = turn on exchange moves (parallel tempering)
   INTEGER(KIND=IB) :: IDIP       !! 1 = allow dipping layers (dip and strike)
   INTEGER(KIND=IB) :: ISD_RV     !! 1 = allow hiararchical sd
   INTEGER(KIND=IB) :: ISD_SWD    !! 1 = allow hiararchical sd
   INTEGER(KIND=IB) :: ISD_ELL    !! 1 = allow hiararchical sd
   INTEGER(KIND=IB) :: ISD_MT     !! 1 = allow hiararchical sd
   INTEGER(KIND=IB) :: ITAPER     !! 1 = taper pred data
   INTEGER(KIND=IB) :: I_ABS_ELL  !! 1 = inverts abs of ellipticities
   INTEGER(KIND=IB) :: I_LOG10_ELL!! 1 = inverts log10 of absolute ellipticities
   INTEGER(KIND=IB) :: IDECON     !! 0 = water-level with raysum, 1 = ITD with raysum

!!
!! Model and data dimensions
!!
   INTEGER(KIND=IB)            :: NDAT_MT      ! Number of MT data (apparent resistivities and phases)
   INTEGER(KIND=IB)            :: NDAT_SWD     ! Number surface wave dispersion data (phases)
   INTEGER(KIND=IB)            :: NMODE        ! Number surface wave dispersion modes
   INTEGER(KIND=IB)            :: NDAT_ELL     ! Number ellipticity data
   INTEGER(KIND=IB)            :: NMODE_ELL    ! Number ellipticity modes
   INTEGER(KIND=IB)            :: NTIME        ! Number time samples
   INTEGER(KIND=IB)            :: NSRC         ! Number time samples in source-time function
   INTEGER(KIND=IB)            :: NTIME2       ! Number time samples for zero padded observations
   INTEGER(KIND=IB)            :: NRRG         ! No. points for convolution RRG
   INTEGER(KIND=IB)            :: NRGRG        ! No. points for convolution RGRG
   INTEGER(KIND=IB)            :: NLMN         ! Min number of layers
   INTEGER(KIND=IB)            :: NLMX         ! Max number of layers
   INTEGER(KIND=IB)            :: NPL          ! No. parameters per layer

   CHARACTER(len=64) :: filebasefile      = 'filebase.txt'
   INTEGER(KIND=IB)                 :: NRF1

!! 
!! Forward RF specific 
!!
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxpar
  LOGICAL,ALLOCATABLE,DIMENSION(:)         :: isoflag
  !!                                                      thick      rho      alph     beta    %P    %S    tr    pl    st    di
  REAL(KIND=SP),DIMENSION(10),PARAMETER:: curmod_glob = (/10000._RP,2600._RP,6000._RP,3600._RP,0._RP,0._RP,0._RP,0._RP,0._RP,0._RP/)
  INTEGER :: raysumfail = 0_IB
  INTEGER :: NTR
  REAL(KIND=SP),ALLOCATABLE,DIMENSION(:)       :: baz2,slow2,sta_dx,sta_dy
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:)    :: nseg
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:,:):: phaselist
  REAL(KIND=SP),ALLOCATABLE,DIMENSION(:,:,:)   :: synth_cart,synth_ph
  INTEGER(KIND=IB)          :: numph,iphase
  INTEGER(KIND=IB)          :: mults  !! = 0  !! Parameter for multiples in layers
  INTEGER(KIND=IB)          :: directs !!= 0  !! Parameter for directse  in layers !!!Pejman
  INTEGER(KIND=IB),PARAMETER:: out_rot = 1  !! 
  INTEGER(KIND=IB),PARAMETER:: align   = 1  !! 
  REAL(KIND=SP)             :: shift2
  REAL(KIND=SP)             :: sampling_dt
  REAL(KIND=SP),PARAMETER   :: sig   = 0.01_SP
  !REAL(KIND=SP),PARAMETER   :: VPVS  = 1.75_RP
  REAL(KIND=SP)             :: width2             !! Set < 0 to return impulse response
  REAL(KIND=SP)             :: wl                 !! water-level for ray3d
  REAL(KIND=RP)             :: hmx                !! Max crustal depth in km
  REAL(KIND=RP),DIMENSION(4):: sdmn               !! Min standard deviation
  REAL(KIND=RP),DIMENSION(4):: sdmx               !! Max standard deviation
  REAL(KIND=SP)             :: tolerance          !! for ITD
  INTEGER(KIND=IB)          :: maxbumps, MAXG     !! for ITD
  REAL(KIND=SP)             :: norm_ITD
!! 
!! Forward gpell specific
!!
  INTEGER(KIND=IB)      :: ELL_verbose
  REAL(KIND=RP)         :: ELL_prec
  INTEGER(KIND=IB)      :: I_SAMPLING_TYPE_ELL
  INTEGER(KIND=IB)      :: I_SET_STEP_ELL
  REAL(KIND=RP)         :: STEP_SIZE_ELL
  INTEGER(KIND=IB)      :: I_SET_COUNT_ELL
  INTEGER(KIND=IB)      :: COUNT_ELL  
  INTEGER(KIND=IB)      :: I_SET_RANGE_ELL
!!

!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlim
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlim
   INTEGER(KIND=IB)            :: kmin     = 0       ! Min number of layers
   INTEGER(KIND=IB)            :: kmax     = 0       ! Max number of layers
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pk       ! Poisson prior on k
   REAL(KIND=RP)               :: lambda             ! Lambda parameter for Poisson prior on k
   REAL(KIND=RP)               :: hmin               ! Min allowed layer thickness
   REAL(KIND=RP),PARAMETER     :: fact     = 1.00_RP ! factor for rotated space perturbation
   REAL(KIND=RP),PARAMETER     :: factdelay= 1.50_RP ! shrinking factor for delayed rejection (>1.)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpert
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsc
   REAL(KIND=RP)               :: vs01,vs02,vsh1,vsh2,numin,numax
   REAL(KIND=RP)               :: area_bn,area_kskp,area_cr
   REAL(KIND=RP)               :: sigmamin, sigmamax 
!!
!!  Autoregressive model prior variables:
!!
   REAL(KIND=RP)              :: armxH      = 0.2_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP)              :: armxV      = 0.2_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP)              :: armxSWD    = 0.5_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP)              :: armxELL    = 0.5_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimar, maxlimar, maxpertar
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsd, pertarsdsc
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimarSWD, maxlimarSWD, maxpertarSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsdSWD, pertarsdscSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimarELL, maxlimarELL, maxpertarELL
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsdELL, pertarsdscELL
   !REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimarMT, maxlimarMT, maxpertarMT
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsdMT, pertarsdscMT
   
!!
!!  Standard deviation prior variables:
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsd, maxlimsd, maxpertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsd, pertsdsdsc
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsdSWD, maxlimsdSWD, maxpertsdSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdSWD, pertsdsdscSWD
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsdELL, maxlimsdELL, maxpertsdELL
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdELL, pertsdsdscELL
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsdMT, maxlimsdMT, maxpertsdMT
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdMT, pertsdsdscMT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MT Z measurement variance 
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: MTZVAR 
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: sdSWD 
   
!   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE:: fr,fstep   ! Total frequency array
!   REAL(KIND=RP)               :: z_t
!   REAL(KIND=RP)               :: cw
!   REAL(KIND=RP)               :: rw
   CHARACTER(len=64) :: filebase
   INTEGER(KIND=IB)  :: filebaselen
   CHARACTER(LEN=100) :: infileV
   CHARACTER(LEN=100) :: infileR
   CHARACTER(LEN=100) :: infileT
   CHARACTER(LEN=100) :: infileSWD
   CHARACTER(LEN=100) :: infile_sdSWD
   CHARACTER(LEN=100) :: infileELL
   CHARACTER(LEN=100) :: infileMT
   CHARACTER(LEN=100) :: infileMT_ZVAR  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CHARACTER(LEN=100) :: infileref
   CHARACTER(LEN=100) :: infileCdi
   CHARACTER(LEN=100) :: infileCdiSWD
   CHARACTER(LEN=100) :: infileCdiELL
   CHARACTER(LEN=100) :: infileCdiMT
   CHARACTER(LEN=100) :: repfileR
   CHARACTER(LEN=100) :: repfileT
   CHARACTER(LEN=100) :: repfileV
   CHARACTER(LEN=100) :: repfileS
   CHARACTER(LEN=100) :: repfileRF
   CHARACTER(LEN=100) :: repfileSWD
   CHARACTER(LEN=100) :: repfileSWDar
   CHARACTER(LEN=100) :: repfileELL
   CHARACTER(LEN=100) :: repfileELLar
   CHARACTER(LEN=100) :: repfileMT
   CHARACTER(LEN=100) :: repfileMTar
   CHARACTER(LEN=100) :: parfile
   CHARACTER(LEN=64) :: logfile
   CHARACTER(LEN=64) :: seedfile
   CHARACTER(len=64) :: mapfile
   CHARACTER(len=64) :: covfile
   CHARACTER(LEN=64)  :: obsfile
   CHARACTER(LEN=64)  :: arfile
   CHARACTER(LEN=64)  :: predfile
   CHARACTER(LEN=64)  :: obsfileSWD
   CHARACTER(LEN=64)  :: arfileSWD
   CHARACTER(LEN=64)  :: predfileSWD
   CHARACTER(LEN=64)  :: obsfileELL
   CHARACTER(LEN=64)  :: arfileELL
   CHARACTER(LEN=64)  :: predfileELL
   CHARACTER(LEN=64)  :: obsfileMT
   CHARACTER(LEN=64)  :: arfileMT
   CHARACTER(LEN=64)  :: predfileMT
   CHARACTER(len=64) :: sdfile
   CHARACTER(len=64) :: samplefile
   CHARACTER(len=64) :: stepsizefile
   CHARACTER(len=64) :: lincovfile
   CHARACTER(LEN=64)  :: modname = 'sample.geom'

!!
!! Velocity reference model
!!
  INTEGER(KIND=IB)                        :: NVELREF, NPREM
  REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: vel_ref, vel_prem
!!
!! Iterative covariance estimate  parameters
!!
   CHARACTER(len=64)                            :: samplefile_covIter
   CHARACTER(len=64)                            :: samplefile_res_covIter
   CHARACTER(LEN=100)                           :: covparfile
   CHARACTER(LEN=100)                           :: infileCd
   CHARACTER(LEN=100)                           :: infileCdSWD
   CHARACTER(LEN=100)                           :: infileCdELL
   CHARACTER(LEN=100)                           :: infileCdMT
   INTEGER(KIND=IB)                             :: usample_covIter, usample_res_covIter, units, unit2
   INTEGER(KIND=IB)                             :: ICOVest                          
   INTEGER(KIND=IB)                             :: icovIter!!, covIter_nsamples
   INTEGER(KIND=IB)                             :: covIter_zero_nsamples, covITER_period, MAXcovIter
   INTEGER(KIND=IB)                             :: ICOV_iterUpdate, ICOV_iterUpdate_RV, ICOV_iterUpdate_SWD, ICOV_iterUpdate_ELL, ICOV_iterUpdate_MT  
   INTEGER(KIND=IB)                             :: cov_converged
   INTEGER(KIND=IB)                             :: ISD_RV_covIter, ISD_SWD_covIter, ISD_ELL_covIter, ISD_MT_covIter
   INTEGER(KIND=IB)                             :: CHAINTHIN_COVest_period_zeroIter, CHAINTHIN_COVest_period_nonzeroIter                               
   INTEGER(KIND=IB)                             :: NRF2, NMODE2, NMODE_ELL2, NMT2, ncount3                                  
   REAL(KIND=RP),DIMENSION(4)                   :: sdmn_covIter               
   REAL(KIND=RP),DIMENSION(4)                   :: sdmx_covIter              
   REAL(KIND=RP),DIMENSION(4)                   :: sdpar_covIter      
   INTEGER(KIND=IB)                             :: NKEEP2, NKEEP_covIter, NKEEP3, NKEEP_covIter_res        
   INTEGER(KIND=IB)                             :: iSAVEsample_covIter, iSAVEsample_only_zeroIter
   INTEGER(KIND=IB)                             :: iMAP_calc           
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: Cd, Cd_old, Cdi_old      
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: CdSWD, CdSWD_old, CdiSWD_old      
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: CdELL, CdELL_old, CdiELL_old      
   COMPLEX(KIND=RP),DIMENSION(:,:),ALLOCATABLE  :: CdMT, CdMT_old, CdiMT_old      
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE     :: sampleDres, sample2      
   INTEGER(KIND=IB)                             :: NsampleDres
   INTEGER(KIND=IB)                             :: iconverge_criterion, iconverge_criterion_RV, iconverge_criterion_SWD, iconverge_criterion_ELL, iconverge_criterion_MT 
   REAL(KIND=RP)                                :: converge_threshold_RV, converge_threshold_SWD, converge_threshold_ELL, converge_threshold_MT 
   REAL(KIND=RP)                                :: covIter_errRV, covIter_errSWD, covIter_errELL, covIter_errMT 
   INTEGER(KIND=IB), DIMENSION(4)               :: cov_converged_datasets
   INTEGER(KIND=IB), DIMENSION(4)               :: ICOViter_datasets     
   INTEGER(KIND=IB)                             :: nfrac_RV, MAX_NAVE_RV, inonstat_RV, iunbiased_RV, imr_RV
   REAL(KIND=RP)                                :: damp_power_RV
   INTEGER(KIND=IB)                             :: nfrac_SWD, MAX_NAVE_SWD, inonstat_SWD, iunbiased_SWD, imr_SWD
   REAL(KIND=RP)                                :: damp_power_SWD
   INTEGER(KIND=IB)                             :: nfrac_ELL, MAX_NAVE_ELL, inonstat_ELL, iunbiased_ELL, imr_ELL
   REAL(KIND=RP)                                :: damp_power_ELL
   INTEGER(KIND=IB)                             :: nfrac_MT, MAX_NAVE_MT, inonstat_MT, iunbiased_MT, imr_MT
   REAL(KIND=RP)                                :: damp_power_MT
!!
!! Parallel Tempering parameters
!!
  INTEGER(KIND=IB)                         :: NPTCHAINS1            !! # chains T=1
  REAL(KIND=RP)                            :: dTlog                 ! Temperature increment
  INTEGER(KIND=IB)                         :: NT                    !! # tempering levels (temperatures)
  INTEGER(KIND=IB)                         :: NPTCHAINS             !! # parallel tempering chains
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: NCHAINT
  INTEGER(KIND=IB)                         :: ncswap     = 0_IB     !! Temp swap accepted counter
  INTEGER(KIND=IB)                         :: ncswapprop = 0_IB     !! Temp swap proposed counter
  REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: beta_pt               !! Temperature array parallel tempering
  INTEGER(KIND=IB)                         :: ibirth  = 0,ideath  = 0
  INTEGER(KIND=IB)                         :: ibirths = 0,ideaths = 0

!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX
   INTEGER(KIND=IB)           :: NFPMX2
   INTEGER(KIND=IB)           :: ioutside   = 0
   INTEGER(KIND=IB)           :: ireject    = 0, iaccept = 0, iaccept_delay = 0, ireject_delay = 0
   INTEGER(KIND=IB)           :: i_bd,i_bds     ! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_sdpert = 0   ! if sigma is perturbed, don't compute forward model
   INTEGER(KIND=IB)           :: ishearfail = 0 ! if sigma is perturbed, don't compute forward model
   INTEGER(KIND=IB)           :: i_ref_nlay = 0 ! if sigma is perturbed, don't compute forward model

!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch slaves
   INTEGER(KIND=IB)       :: iconv2   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iconv3   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iarfail  = 0       ! Tracks failure of AR model when predicted AR series too large

!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB),PARAMETER    :: NCHAIN     = 1E9_IB  ! # iterations (max # MCMC steps)
   INTEGER(KIND=IB)              :: ICHAINTHIN = 1E0_IB  ! Chain thinning interval
   INTEGER(KIND=IB)              :: NKEEP      = 1E1_IB  ! Number models to keep before writing
   INTEGER(KIND=IB),PARAMETER    :: NAP        = 10      ! Misc parameters in sample (for bookeeping)
   INTEGER(KIND=IB),PARAMETER    :: NDM        = 100     ! No. steps in lin rot est
   INTEGER(KIND=IB)          :: TCHCKPT              !! No. seconds (integer value) between checkpoints
   INTEGER(KIND=IB)          :: icheckpoint          !! No. of checkpoints to data (read from checkpoint/status.txt)

   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)   :: sdevm  ! Std dev for perturbations
!!
!!  Structures for objects and data 
!!
  INTEGER :: imcmc1 = 1   !! Counter for models at T=1 (needs to survive checkpointing!)
  INTEGER :: imcmc2 = 1   !! Counter for mcmc steps to scale diminishing adaptation (needs to survive checkpointing!)
  INTEGER :: NFIELD = 62  !! The number of fields in objstruc
  INTEGER :: objtype1     !! Name of objtype for MPI sending
  INTEGER :: objtype2     !! Name of objtype for MPI sending
  INTEGER :: objtype3     !! Name of objtype for MPI sending
   TYPE :: objstruc
      INTEGER(KIND=IB)                        :: k          ! No. nodes
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: voro       ! 1D voronoi nodes
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:):: voroidx    ! 1D voronoi nodes
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     ! Forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: hiface
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: ziface
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparR        !! Std dev Radial comp
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparV        !! Std dev Vertical comp
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparT        !! Std dev Vertical comp
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparSWD      !! Std dev SWD data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparELL      !! Std dev ELL data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparMT       !! Std dev MT data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveH        !! Std dev Radial comp
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveV        !! Std dev Vertical comp
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveT        !! Std dev Vertical comp
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveSWD      !! Std dev SWD data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveELL      !! Std dev ELL data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveMT       !! Std dev MT data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arpar         !! AR model forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparSWD      !! AR model forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparELL      !! AR model forward parameters
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxar      !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarSWD      !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarELL      !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: gvoroidx   !! Index of live parameters on birth/death node
      INTEGER(KIND=IB)                        :: nunique     !! Number of unique interfaces
      INTEGER(KIND=IB)                        :: NFP         !! Number forward parameters
      REAL(KIND=RP)                           :: beta
      REAL(KIND=RP)                           :: logL        !! log likelihood
      REAL(KIND=RP)                           :: logPr       !! log Prior probability ratio
      REAL(KIND=RP)                           :: tcmp
      INTEGER(KIND=IB)                        :: ireject_bd = 0
      INTEGER(KIND=IB)                        :: iaccept_bd = 0
      INTEGER(KIND=IB)                        :: ireject_bds = 0
      INTEGER(KIND=IB)                        :: iaccept_bds = 0
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsR       !! Observed data for one logL eval.
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredR      !! Predicted data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsV       !! Observed data for one logL eval.
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredV      !! Predicted data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsT       !! Observed data for one logL eval.
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredT      !! Predicted data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: S           !! Predicted data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresR       !! Data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresV       !! Data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresT       !! Data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarR        !! Autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarV        !! Autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarT        !! Autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsSWD     !! Observed data SWD
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredSWD    !! Predicted SWD data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresSWD     !! SWD data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarSWD      !! SWD autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: periods     !! SWD autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsELL     !! Observed data ELL
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredELL    !! Predicted ELL data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresELL     !! ELL data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarELL      !! ELL autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: periods_ELL !! ELL autoregressive model predicted data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)  :: DobsMT      !! Observed data MT
      !REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)  :: MTZVAR      !! Observed data MT impedance (Z) relative variances !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)  :: DpredMT     !! Predicted MT data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)  :: DresMT      !! MT data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)  :: freqMT      !! MT data frequancy
   END TYPE objstruc
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE   :: taper_dpred    !! Taper for pred data

!!!!!!!!!!!!!!Pejman: variables for making the training dataset for the nueral network
REAL(KIND=RP), ALLOCATABLE, DIMENSION(:,:)   :: cur_modMT
REAL(KIND=RP), ALLOCATABLE, DIMENSION(:,:)   :: cur_modSeismic
!!!!!!!!!!!!!!

!!
!! Structure for covariance matrices (only applies for ICOV >= 2)
!!
!  TYPE :: covstruc
!    REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: Cdi   ! Inverse covariance matrix
!  END TYPE covstruc
!!
!! File units for IO
!!

!!
!! Covariance matrices (only applies for ICOV >= 2)
!!
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: Cdi     ! Inverse covariance matrix 
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: CdiSWD  ! Inverse covariance matrix SWD data
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: CdiELL  ! Inverse covariance matrix ELL data
COMPLEX(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: CdiMT  ! Inverse covariance matrix MT data
!!

  INTEGER(KIND=IB) :: usample, ustep, ulog
  INTEGER(KIND=IB) :: reclen

   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE      :: icount

!!
!!  Buffer for sdave for AR discrimination in CHECKBOUNDS_AR
!!
   INTEGER(KIND=IB), PARAMETER               :: NBUF = 100
   REAL(KIND=RP),DIMENSION(:,:,:),ALLOCATABLE:: sdbuf
!!
!!  Global variables
!!
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE   :: sample                    ! Posterior sample
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: tmpmap                    ! temporary for reading map
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: buf_save_snd,buf_save_rcv ! Buffers for MPI sending
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE  :: buffer1                   !
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: buffer2,buffer3           !

!!
!!  MPI global variables
!!
   INTEGER(KIND=IB)            :: rank,NTHREAD,ierr
   INTEGER(KIND=IB)            :: ncount1,ncount2
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB
   INTEGER                     :: to,from,COMM
   INTEGER                     :: status(MPI_STATUS_SIZE)
   INTEGER(KIND=IB)            :: isize1,isize2,isize3
   INTERFACE
      FUNCTION RANDPERM(num)
         USE data_type, ONLY : IB
         IMPLICIT NONE
         INTEGER(KIND=IB), INTENT(IN) :: num
         INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
      END FUNCTION RANDPERM
   END INTERFACE
   REAL(KIND=RP) :: tsave1, tsave2              ! Overall time 

!!
!!  FFTW stuff
!!
   INTEGER*8 :: planR,planC

  CONTAINS
  !==============================================================================
  integer function newunit(unit)
  !==============================================================================
    integer, intent(out), optional :: unit
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
      inquire(unit=lun,opened=opened)
      if (.not. opened) then
        newunit=lun
        exit
      end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

END MODULE RJMCMC_COM
!=======================================================================
