!=======================================================================
MODULE RAY3D_COM
  USE RAY3D_PAR
  USE DATA_TYPE
  IMPLICIT NONE
  INTEGER :: nlay

  INTEGER:: nsamp,misfit_type
  REAL::    mod_min(maxlay,7)
  REAL::    mod_max(maxlay,7)
  LOGICAL:: mod_flag(maxlay,7)
  REAL::    h(maxlayd),rho(maxlayd),alpha(maxlayd),vpvs(maxlayd),beta(maxlayd), &
            strike(maxlayd),dip(maxlayd)
  REAL:: baz(maxtr),slow(maxtr)
  REAL:: data_rf(2,maxsamp,maxtr),weight(maxtr)
  REAL(KIND=SP)          :: shift       = 5._RP
!  REAL(KIND=SP)          :: width       = sampling_dt/4._RP
  REAL(KIND=SP)          :: width       = -1._RP
  REAL(KIND=SP)          :: dt,time_window,t_af
! Common block for raysum calling routines

! Common block for rays3d calling routines
!        common /ray3d_com/nlay,mod_min,mod_flag,baz,slow,
!     *     ntr,data_rf,dt,shift,t_af,nsamp,weight,width,
!     *     w1,w2,c,time_window,misfit_type
!


!	Array sizes for neighbourhood algorithm routines
!
!	PARAMETER		MEANING
!                                       The NA routines use the
!                                       following include file to define
!                                       all parameters controlling
!                                       memory required by the arrays.
!
!                                       The following is a description
!                                       of each parameter in na_param.in
!                                       that you are allowed to change to suit
!                                       your application. They are
!                                       checked against input parameters
!                                       read in from file na.in and an 
!                                       error is reported if any values
!                                       are too small.
!                                       
!       PARAMETER               MEANING
!
!       nd_max                  Maximum number of variables, i.e. the 
!                               dimension of the parameter space.
!
!       nit_max                 Maximum number of iterations
!
!       nsample_max             Maximum number of models generated
!                               in each sample
!
!       nh_max                  Maximum size of header of NAD output file
!                               (leave unchanged if you choose not to add
!                                any material to the NAD header. See manual)
!
!	maxseq			Maximum number of random sequences generated
!				by quasi random number generator. This
!				value can be set to 1 if the quasi
!				random number generator is not used.
!
!-----------------------------------------------------------------------
!
!				The following parameters are fixed
!				and should not be changed.
!
!	nmod_max		Maximum number of models to be generated
!			        (determined by nsample_max and nit_max)
!
!	nsleep_max		Maximum number of samples skipped over 
!				(Currently fixed at 1. Do not change)
!
!-----------------------------------------------------------------------
!
       INTEGER, parameter:: nsample_max=100, nit_max=5000
!sed    parameter (nd_max=27  )
       INTEGER, parameter:: nd_max=50,nh_max=1000,nsleep_max=1
       INTEGER, parameter:: maxseq=nd_max*nsample_max
!       parameter (maxseq=50  )
       INTEGER, parameter:: nmod_max = nsample_max*(nit_max+1)
END MODULE RAY3D_COM
!=======================================================================
