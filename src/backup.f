!!=======================================================================

SUBROUTINE LOGLHOOD(obj,ipred)
!!=======================================================================
!! This is LOGLHOOD function that combines various data sets for 
!! joint inversion. 
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INCLUDE 'raysum/params.h'
INTEGER(KIND=IB):: ipred
TYPE (objstruc) :: obj
REAL(KIND=RP)   :: logL_RV, logL_SWD
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod

!!
!! Generate input vector for SWD
!!
!! curmod = (/ thick, rho, alph, beta, %P, %S, tr, pl, st, di /)
CALL MAKE_CURMOD(obj,curmod)

!! log likelihood RV data
IF(I_RV == 1)THEN
  !! Radial and vertical component inversion
  CALL LOGLHOOD_RV(obj,curmod,ipred,logL_RV)
ELSEIF(I_RV == -1)THEN
  !! Receiver function inversion
  CALL LOGLHOOD_RF(obj,curmod,ipred,logL_RV)
ELSE
  logL_RV = 0._RP
ENDIF
!! log likelihood SWD data
IF(I_SWD == 1)THEN
  CALL LOGLHOOD_SWD(obj,curmod,ipred,logL_SWD)
ELSE
  logL_SWD = 0._RP
ENDIF

!!
!! Joint likelihood (assumes independent errors on the vasious data sets)
!!
obj%logL = logL_RV + logL_SWD
!IF(obj%logL == 0.)THEN
!  CALL PRINTPAR(obj)
!  STOP
!ENDIF

RETURN
END SUBROUTINE LOGLHOOD
!!=======================================================================

SUBROUTINE LOGLHOOD_RV(obj,curmod,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted R and V components
!!
USE RJMCMC_COM
USE MPI
USE ieee_arithmetic
!USE NR
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,ipred,iaz
INTEGER(KIND=IB)                 :: ierr_neg,ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpH,EtmpV,EtmpT
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
!! Forward specific variables
REAL(KIND = SP)                  :: aa(3,3,3,3,NLMX),ar_list(3,3,3,3,NLMX,2),rot(3,3,NLMX)
REAL(KIND = SP)                  :: amp_in,tt(maxph,maxtr),amp(3,maxph,maxtr)
LOGICAL :: bailout,errflag
CHARACTER(len=64) :: phname
LOGICAL :: ISNAN
REAL(KIND=RP) :: DpredHG(NRF1,NTIME),DpredVG(NRF1,NTIME),DpredTG(NRF1,NTIME)
REAL(KIND=RP) :: DpredH(NTIME2),DpredV(NTIME2),DpredT(NTIME2)
REAL(KIND=RP) :: HHG(NHHG),VVG(NHHG),HGHG(NHGHG),VGVG(NHGHG),TTG(NHHG),TGTG(NHGHG)
REAL(KIND=RP) :: S(NSRC)

!!
!!  Build a_ijkl tensors
!!
call buildmodel(aa,ar_list,rot,curmod(:,1),curmod(:,2),&
     curmod(:,3),curmod(:,4),isoflag,curmod(:,5),curmod(:,6),&
     curmod(:,7),curmod(:,8),curmod(:,9),curmod(:,10),obj%nunique+1)

!!
!! Compute multiples
!!
width = 0.5
numph = 0
iphase = 1
!if (mults .eq. 0) then
!  call ph_direct(phaselist,nseg,numph,obj%nunique+1,iphase)
!end if
!!! This is irrelevant, looks at lowermost interface only...
!!if (mults .eq. 1) then
!!  ilay = obj%nunique
!!  call ph_direct_rf(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
!!  call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
!!end if
!if (mults .eq. 2) then
!  do ilay=1,obj%nunique
!    call ph_direct_rf(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
!  end do
!  do ilay=1,obj%nunique
!    call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
!  end do
!end if

 numph=0
  call ph_direct(phaselist,nseg,numph,obj%nunique+1,iphase)
  do ilay=1,obj%nunique
    call ph_fsmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
  end do


!!
!! Get arrival times and amplitudes
!!
amp_in = 1._RP
ierr_neg = 0
call get_arrivals(tt,amp,curmod(:,1),curmod(:,2),isoflag,&
     curmod(:,9),curmod(:,10),aa,ar_list,rot,baz2,slow2,sta_dx,&
     sta_dy,phaselist,ntr,nseg,numph,obj%nunique+1,amp_in,errflag,ierr_neg)
IF(errflag) raysumfail = raysumfail + 1
!IF(ierr_neg == 1) raysumfail = raysumfail + 1
IF(ierr_neg == 1)THEN
  PRINT*,'ERROR'
ENDIF

!! Normalize amplitudes by direct P
call norm_arrivals(amp,baz2,slow2,curmod(1,3),curmod(1,4),&
                   curmod(1,2),ntr,numph,1,1)
!!
!! Make traces (in NEZ) and rotate into RTV (out_rot = 1)
!!
call make_traces(tt,amp,ntr,numph,NTIME,sampling_dt,width,align,shift,synth_cart)
!if (out_rot .eq. 1) then
call rot_traces(synth_cart,baz2,ntr,NTIME,synth_ph)
!end if
!if (out_rot .eq. 2) then
!  call fs_traces(synth_cart,baz2,slow2,curmod(:,3),curmod(:,4),&
!       curmod(:,2),ntr,NTIME,synth_ph)
!end if

!! GF predictions:
!! GF predictions:
DO iaz = 1,NRF1
  DpredHG(iaz,1:NTIME)  = REAL(synth_ph(1,1:NTIME,iaz),RP)
  DpredVG(iaz,1:NTIME)  = REAL(synth_ph(3,1:NTIME,iaz),RP)
  DpredTG(iaz,1:NTIME)  = REAL(synth_ph(2,1:NTIME,iaz),RP)
ENDDO
IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='DpredRG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=51,FILE='DpredVG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=52,FILE='DpredTG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  DO iaz = 1,NRF1
    WRITE(50,201) DpredHG(iaz,:)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(51,201) DpredVG(iaz,:)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(52,201) DpredTG(iaz,:)
  ENDDO
  CLOSE(50)
  CLOSE(51)
  CLOSE(52)
  201 FORMAT(1024ES18.10)
ENDIF
DO iaz = 1,NRF1
  !! Time domain convolution for source equation terms:
  !! SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
  CALL CONVT(obj%DobsH(iaz,:),DpredHG(iaz,:),HHG,NTIME2,NTIME,NHHG)
  CALL CONVT(obj%DobsV(iaz,:),DpredVG(iaz,:),VVG,NTIME2,NTIME,NHHG)
  CALL CONVT(DpredHG(iaz,:),DpredHG(iaz,:),HGHG,NTIME,NTIME,NHGHG)
  CALL CONVT(DpredVG(iaz,:),DpredVG(iaz,:),VGVG,NTIME,NTIME,NHGHG)
  IF(I_T == 0)THEN
    CALL CONVT(obj%DobsT(iaz,:),DpredTG(iaz,:),TTG,NTIME2,NTIME,NHHG)
    CALL CONVT(DpredTG(iaz,:),DpredTG(iaz,:),TGTG,NTIME,NTIME,NHGHG)
  ENDIF

  !! Estimate source by time domain dconvolution 
  !! (accounting for different std dev on H and V):
  IF(ICOV == 1)THEN
    IF(I_T == 0)THEN
      CALL DCONVT(HHG/obj%sdparH(iaz)**2+VVG/obj%sdparV(iaz)**2,HGHG/obj%sdparH(iaz)**2+VGVG/obj%sdparV(iaz)**2,S,NHHG,NHGHG,NSRC)
    ELSEIF(I_T == 1)THEN
      CALL DCONVT(HHG/obj%sdparH(iaz)**2+VVG/obj%sdparV(iaz)**2+TTG/obj%sdparT(iaz)**2, &
                  HGHG/obj%sdparH(iaz)**2+VGVG/obj%sdparV(iaz)**2+TGTG/obj%sdparT(iaz)**2,S,NHHG,NHGHG,NSRC)
    ENDIF
  ELSE
    PRINT*,'WARNING ICOV problem'
    CALL DCONVT(HHG+VVG,HGHG+VGVG,S,NHHG,NHGHG,NSRC)
  ENDIF
  !! Data predictions via time domain convolution:
  CALL CONVT(S,DpredHG(iaz,:),DpredH,NSRC,NTIME,NTIME2)
  CALL CONVT(S,DpredVG(iaz,:),DpredV,NSRC,NTIME,NTIME2)
  IF(I_T == 1)CALL CONVT(S,DpredTG(iaz,:),DpredT,NSRC,NTIME,NTIME2)

  obj%S(iaz,1:NSRC) = S(1:NSRC)
  obj%DpredH(iaz,1:NTIME) = DpredH(1:NTIME)
  obj%DpredV(iaz,1:NTIME) = DpredV(1:NTIME)
  obj%DresH(iaz,1:NTIME) = obj%DobsH(iaz,1:NTIME)-obj%DpredH(iaz,1:NTIME)
  obj%DresV(iaz,1:NTIME) = obj%DobsV(iaz,1:NTIME)-obj%DpredV(iaz,1:NTIME)
  IF(I_T == 1)THEN
    obj%DpredT(iaz,1:NTIME) = DpredT(1:NTIME)
    obj%DresT(iaz,1:NTIME) = obj%DobsT(iaz,1:NTIME)-obj%DpredT(iaz,1:NTIME)
  ENDIF
ENDDO

ibadlogL = 0
!IF(IAR == 1)THEN
!   obj%Dar  = 0._RP
!   !!
!   !!  Compute autoregressive model
!   !!
!   CALL ARPRED(obj,1,1,NTIME)
!   !! Recompute predicted data as ith autoregressive model
!   obj%Dres = obj%Dres-obj%Dar
!
!   !! Check if predicted AR model data are outside max allowed bounds
!   CALL CHECKBOUNDS_ARMX(obj,ibadlogL)
!ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpH = 0._RP
  EtmpV = 0._RP
  EtmpT = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpH(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresH(iaz,:)**2._RP)/(2._RP*obj%sdparH(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparH(iaz)))
      EtmpV(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresV(iaz,:)**2._RP)/(2._RP*obj%sdparV(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparV(iaz)))
      IF(I_T == 1)THEN
        EtmpT(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                    -(SUM(obj%DresT(iaz,:)**2._RP)/(2._RP*obj%sdparT(iaz)**2._RP)&
                    +REAL(NTIME,RP)*LOG(obj%sdparT(iaz)))
      ENDIF
    ENDDO
  ENDIF
  logL = SUM(EtmpH) + SUM(EtmpV)
  IF(I_T == 1)logL = logL + SUM(EtmpT)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
  logL = -HUGE(1._RP)
  ibadlogL = 0
ENDIF

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_RV
!=======================================================================

SUBROUTINE LOGLHOOD_RF(obj,curmod,ipred,logL)
!=======================================================================
USE RJMCMC_COM
USE MPI
USE ieee_arithmetic
!USE NR
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,iparcur,ipred,iaz
INTEGER(KIND=IB)                 :: ierr_neg,ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpH
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
!! Forward specific variables
REAL(KIND = SP)                  :: aa(3,3,3,3,NLMX),ar_list(3,3,3,3,NLMX,2),rot(3,3,NLMX)
REAL(KIND = SP)                  :: amp_in,tt(maxph,maxtr),amp(3,maxph,maxtr)
LOGICAL :: bailout,errflag
CHARACTER(len=64) :: phname
LOGICAL :: ISNAN

!!
!!  Compute predicted RF
!!
!!  Build a_ijkl tensors
!!
call buildmodel(aa,ar_list,rot,curmod(:,1),curmod(:,2),&
     curmod(:,3),curmod(:,4),isoflag,curmod(:,5),curmod(:,6),&
     curmod(:,7),curmod(:,8),curmod(:,9),curmod(:,10),obj%nunique+1)

numph = 0
iphase = 1
if (mults .eq. 0) then
  call ph_direct(phaselist,nseg,numph,obj%nunique+1,iphase)
end if
!! This is irrelevant, looks at lowermost interface only...
!if (mults .eq. 1) then
!  ilay = obj%nunique
!  call ph_direct_rf(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
!  call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
!end if
if (mults .eq. 2) then
  do ilay=1,obj%nunique
    call ph_direct_rf(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
  end do
  do ilay=1,obj%nunique
    call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
  end do
end if
!!
!! Get arrival times and amplitudes
!!
amp_in = 1._RP
ierr_neg = 0
call get_arrivals(tt,amp,curmod(:,1),curmod(:,2),isoflag,&
     curmod(:,9),curmod(:,10),aa,ar_list,rot,baz2,slow2,sta_dx,&
     sta_dy,phaselist,ntr,nseg,numph,obj%nunique+1,amp_in,errflag,ierr_neg)
IF(errflag) raysumfail = raysumfail + 1

!! Normalize amplitudes by direct P
call norm_arrivals(amp,baz2,slow2,curmod(1,3),curmod(1,4),&
                   curmod(1,2),ntr,numph,1,1)

!! Make traces
call make_traces(tt,amp,ntr,numph,NTIME,sampling_dt,width,align,shift,synth_cart)
!if (out_rot .eq. 1) then
call rot_traces(synth_cart,baz2,ntr,NTIME,synth_ph)
!end if
!if (out_rot .eq. 2) then
!  call fs_traces(synth_cart,baz2,slow2,curmod(1,3),curmod(1,4),&
!       curmod(1,2),ntr,NTIME,synth_ph)
!end if

!! GF predictions:
obj%DpredH(1,1:NTIME)  = REAL(synth_ph(1,1:NTIME,1),RP)!/MAXVAL(synth_ph(1,1:NTIME,1))
obj%DresH(1,1:NTIME) = obj%DobsH(1,1:NTIME)-obj%DpredH(1,1:NTIME)
IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='DpredRF.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  WRITE(50,201) (obj%DpredH(1,ifr),ifr=1,NTIME)
  CLOSE(50)
  201 FORMAT(ES18.10)
ENDIF

ibadlogL = 0
IF(IAR == 1)THEN
   obj%DarH  = 0._RP
   !!
   !!  Compute autoregressive model
   !!
   CALL ARPRED_RF(obj,1,1,NTIME)
   !! Recompute predicted data as ith autoregressive model
   obj%DresH = obj%DresH-obj%DarH

   !! Check if predicted AR model data are outside max allowed bounds
   CALL CHECKBOUNDS_ARMXRF(obj)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpH = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpH(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresH(iaz,:)**2._RP)/(2._RP*obj%sdparH(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparH(iaz)))
    ENDDO
  ENDIF
  logL = SUM(EtmpH)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
  logL = -HUGE(1._RP)
  ibadlogL = 0
ENDIF

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_RF
!!=======================================================================

SUBROUTINE LOGLHOOD_SWD(obj,curmod,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted SWD data and compute logL_SWD
!!
USE RJMCMC_COM
USE MPI
USE ieee_arithmetic
!USE NR
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                      :: ipred,ierr_swd,imod,ibadlogL,ilay
TYPE (objstruc)                       :: obj
REAL(KIND=SP),DIMENSION(maxlay,10)    :: curmod
REAL(KIND=SP),DIMENSION(maxlay+NPREM,10):: curmod2
REAL(KIND=SP),DIMENSION(NDAT_SWD)     :: periods,DpredSWD
REAL(KIND=RP),DIMENSION(NMODE)        :: EtmpSWD
REAL(KIND=RP)                         :: logL,factvs,factvpvs
REAL(KIND=RP)                         :: tstart, tend, tcmp   ! Overall time 
LOGICAL :: ISNAN

!!
!!  curmod = thick  rho  alph  beta  %P  %S  tr  pl  st  dip 
!!
factvs   = obj%voro(obj%k,2)  !! Vs is in km/s here
factvpvs = obj%voro(obj%k,3)
curmod2  = 0.
curmod2(1:obj%nunique+1,:) = curmod(1:obj%nunique+1,:)
!! last layer thickness
curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%voro(obj%k,1)*1000.

curmod2(obj%nunique+2:obj%nunique+NPREM,1)   = (vel_prem(1,2:NPREM)-vel_prem(1,1:NPREM-1)) * 1000. !! thickness in km
curmod2(obj%nunique+1+NPREM,1)               = 0.                                                !! HS thickness is 0 
curmod2(obj%nunique+2:obj%nunique+1+NPREM,2) = vel_prem(4,1:NPREM) * 1000.   !! Density
curmod2(obj%nunique+2:obj%nunique+1+NPREM,4) = (vel_prem(2,1:NPREM) + factvs) * 1000. !! Vs
curmod2(obj%nunique+2:obj%nunique+1+NPREM,3) = (vel_prem(2,1:NPREM)+factvs)*(vel_prem(3,1:NPREM)+factvpvs)*1000.  !! Vp
!IF(IMAP == 1)THEN
!  WRITE(*,*) 'curmod2 (including PREM)'
!  DO ilay=1,obj%nunique+1+NPREM+1
!    WRITE(*,206)ilay,curmod2(ilay,1:10)
!  ENDDO
!  206   FORMAT(I3,10F12.4)
!ENDIF

periods = REAL(obj%periods(1,1:NDAT_SWD),SP)
!!
!!  Need to append PREM perturbed by half-space perturbation here to 
!!  ensure that long period SWD can be properly modelled. 
!!
CALL dispersion(obj%nunique+1+NPREM,curmod2(1:obj%nunique+1+NPREM,2)/1000., & 
     curmod2(1:obj%nunique+1+NPREM,3)/1000.,curmod2(1:obj%nunique+1+NPREM,4)/1000.,&
     curmod2(1:obj%nunique+1+NPREM,1)/1000.,DpredSWD,&
     periods,NDAT_SWD,ierr_swd)
!CALL dispersion(obj%nunique+1,curmod(1:obj%nunique+1,2)/1000., & 
!     curmod(1:obj%nunique+1,3)/1000.,curmod(1:obj%nunique+1,4)/1000.,&
!     curmod(1:obj%nunique+1,1)/1000.,DpredSWD,&
!     periods,NDAT_SWD,ierr_swd)

IF(ierr_swd /= 0)THEN
  !write(*,*)'THIS MODEL IS WEIRD, Cannot compute dispersion'
  !write(*,*)'ier=',ierr_swd,'rank',rank
  !write(*,*)'ier < 0 ; INPUT ERROR'                                 
  !write(*,*)'ier      = 0 ; NO ERROR '                              
  !write(*,*)'ier       = 1 ; SLOW CONVERGENCE'                      
  !write(*,*)'ier        = 2 ; ROOT NOT FOUND'
  !write(*,*)'-----------------------------------------------'
  !stop
  !CALL PRINTPAR(obj)
  logL = -HUGE(1._RP)
  RETURN
ENDIF

obj%DpredSWD(1,1:NDAT_SWD) = REAL(DpredSWD,RP)
obj%DresSWD(1,1:NDAT_SWD) = obj%DobsSWD(1,1:NDAT_SWD)-obj%DpredSWD(1,1:NDAT_SWD)

ibadlogL = 0
IF(IAR == 1)THEN
  obj%DarSWD  = 0._RP
  !!
  !!  Compute autoregressive model
  !!
  CALL ARPRED_SWD(obj,1,1,NDAT_SWD)
  !! Recompute predicted data as ith autoregressive model
  obj%DresSWD = obj%DresSWD-obj%DarSWD

  !! Check if predicted AR model data are outside max allowed bounds
  CALL CHECKBOUNDS_ARMXSWD(obj,ibadlogL)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpSWD = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one per mode)
    !!
    DO imod = 1,NMODE
      EtmpSWD(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_SWD,RP)/2._RP)) &
                    -(SUM(obj%DresSWD(imod,:)**2._RP)/(2._RP*obj%sdparSWD(imod)**2._RP)&
                    +REAL(NDAT_SWD,RP)*LOG(obj%sdparSWD(imod)))
    ENDDO
  ENDIF
  logL = SUM(EtmpSWD)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
   logL = -HUGE(1._RP)
   ibadlogL = 0
ENDIF

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_SWD
!!=======================================================================

SUBROUTINE MAKE_CURMOD(obj,curmod)
!!=======================================================================
!!
!! Build curmod from model
!!
USE RJMCMC_COM
USE MPI
USE ieee_arithmetic
!USE NR
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ipar,ilay,id,iparcur
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10):: curmod

curmod = 0._RP

id=0
DO ilay=1,obj%nunique+1
  !curmod(ilay,:) = curmod_glob
  DO ipar=1,NPL
    id=id+1
    iparcur = idxpar(ipar)
    IF(ilay <= obj%nunique)THEN
      IF(iparcur == 1)THEN
        curmod(ilay,iparcur)    = REAL(obj%hiface(ilay),SP)  !! Model works with layer thickness!
      ELSE
        curmod(ilay,iparcur)    = REAL(obj%par(id),SP)
      ENDIF
    ELSE
      IF(iparcur == 1)THEN
        curmod(ilay,iparcur) = 0._SP
        id=id-1
      ELSE
        curmod(ilay,iparcur)    = REAL(obj%par(id),SP)
      ENDIF
    ENDIF
  ENDDO
ENDDO
IF(I_VPVS == 1)THEN
!  IF(I_SMP_VS == 0)THEN
!    !! for sampling Vp:
!    curmod(1:obj%nunique+1,1:3) = curmod(1:obj%nunique+1,1:3)*1000._SP
!    curmod(1:obj%nunique+1,4)   = curmod(1:obj%nunique+1,3)/curmod(1:obj%nunique+1,4)
!  ELSE
    !! for sampling Vs:
    curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
    curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
    curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,3)*curmod(1:obj%nunique+1,4)
!  ENDIF
ELSEIF(I_VPVS == -1)THEN
  IF(I_SMP_VS == 0)THEN
    !! for sampling Vp:
    curmod(1:obj%nunique+1,1:3) = curmod(1:obj%nunique+1,1:3)*1000._SP
    curmod(1:obj%nunique+1,4)   = curmod(1:obj%nunique+1,3)/VPVS
  ELSE
    !! for sampling Vs:
    curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
    curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
    curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,3)*VPVS
  ENDIF
ELSE
  curmod(1:obj%nunique+1,1:4) = curmod(1:obj%nunique+1,1:4)*1000._SP
ENDIF

!! Birch Law for analytical density relationship
!curmod(1:obj%nunique+1,2)   = 1000.*0.77+0.32*curmod(:,3)
!! What Thomas uses:
curmod(1:obj%nunique+1,2) = (2.35+0.036*((curmod(1:obj%nunique+1,3)/1000.)-3.0)**2.)*1000.
IF(IMAP == 1)THEN
  DO ilay=1,obj%nunique+1
    WRITE(*,206)ilay,curmod(ilay,1:10)
  ENDDO
  206   FORMAT(I3,10F12.4)
ENDIF
RETURN
END SUBROUTINE MAKE_CURMOD
!!=======================================================================

SUBROUTINE INTERPLAYER_novar(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! This does not allow for variable layer complexity.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE MPI
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP)                    :: vref,vpvsref

!  obj%k       = 3
!  obj%NFP     = (obj%k * NPL) + (NPL-1)
!  obj%voro    = 0._RP
!  obj%voroidx = 0
!  !! True parameters for simulation:
!  obj%voro(1,:) = (/  0.0_RP, 6.0_RP, 1.8_RP,  0.0_RP/)
!  obj%voro(2,:) = (/ 10.0_RP, 7.0_RP, 0.0_RP,  5.0_RP/)
!  obj%voro(3,:) = (/ 30.0_RP, 8.0_RP, 0.0_RP, 30.0_RP/)
!  obj%voroidx(1,:) = (/ 1, 1, 1, 0/)
!  obj%voroidx(2,:) = (/ 1, 1, 0, 1/)
!  obj%voroidx(3,:) = (/ 1, 1, 0, 1/)

!PRINT*,''
!PRINT*,''
!PRINT*,'Before:',rank
!DO ilay = 1,obj%k
!  WRITE(*,203)ilay,obj%voro(ilay,:)
!ENDDO
!PRINT*,''
!PRINT*,''

!! Sort node according to increasing depth:
obj%nunique = obj%k-1
CALL QSORTC2D(obj%voro(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%ziface = 0._RP
obj%ziface(1:obj%k-1) = obj%voro(2:obj%k,1)

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%k-1
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%k,1) = obj%ziface(1:obj%k)
!partmp(1:obj%k,2:NPL) = obj%voro(1:obj%k,2:NPL)

!!
!! Apply reference profile:
!!
IF(I_VREF == 1)THEN
  DO ivo = 1,obj%k
    CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
    partmp(ivo,2) = vref + obj%voro(ivo,2)
    partmp(ivo,3) = vpvsref + obj%voro(ivo,3)
    !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
  ENDDO
ENDIF

obj%par = 0._RP
DO ilay = 1,obj%k-1
  obj%par((ilay-1)*NPL+1:ilay*NPL) = partmp(ilay,:)
ENDDO
obj%par((obj%k-1)*NPL+1:(obj%k*NPL)-1) = partmp(ilay,2:)

203 FORMAT(I4,20F8.2)

END SUBROUTINE INTERPLAYER_novar
!!=======================================================================

SUBROUTINE INTERPLAYER(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! I.e., some layers are duplicates when nodes are not populated:
!!  Parameter 1:  Parameter 2:     Parameter 3:
!!   --o------      ---o---       ------o--------
!!     |               |                |
!!   ------o--      -------       --o------------
!!         |           |            |
!!         |           |            |
!!         |           |            |
!!   ---------      -------       ------------o--
!!         |           |                      |
!!   ---o-----      -------       ---------------
!!      |              |                      |
!!      |              |                      |
!!   ---------      -------       ---------------
!!
!! The layer node always defines the volume partition below the 
!! node position. The node position defines the interface.
!! The first layer is always fixed at 0 and all nodes are populated.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE MPI
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(NPL-1)  :: niface
INTEGER(KIND=RP),DIMENSION(NPL*NLMX):: ifaceidx
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP),DIMENSION(NLMX,NPL-1)  :: voroh
REAL(KIND=RP),DIMENSION(NLMX+1,NPL-1,2):: voroz

REAL(KIND=RP),DIMENSION(NLMX,NPL)::  voro(NLMX,NPL),voroidx(NLMX,NPL)
REAL(KIND=RP)                    :: vref,vpvsref

!  obj%k       = 3
!  obj%NFP     = (obj%k * NPL) + (NPL-1)
!  obj%voro    = 0._RP
!  obj%voroidx = 0
!  !! True parameters for simulation:
!  obj%voro(1,:) = (/  0.0_RP, 6.0_RP, 1.8_RP,  0.0_RP/)
!  obj%voro(2,:) = (/ 10.0_RP, 7.0_RP, 0.0_RP,  5.0_RP/)
!  obj%voro(3,:) = (/ 30.0_RP, 8.0_RP, 0.0_RP, 30.0_RP/)
!  obj%voroidx(1,:) = (/ 1, 1, 1, 0/)
!  obj%voroidx(2,:) = (/ 1, 1, 0, 1/)
!  obj%voroidx(3,:) = (/ 1, 1, 0, 1/)

!PRINT*,''
!PRINT*,''
!PRINT*,'Before:',rank
!DO ilay = 1,obj%k
!  WRITE(*,203)ilay,obj%voro(ilay,:)
!ENDDO
!PRINT*,''
!PRINT*,''

!obj%voroidx(1,:) = (/ 1, 1, 1/)
!obj%voroidx(2,:) = (/ 1, 1, 0/)
!obj%voroidx(3,:) = (/ 1, 1, 0/)

IF(IDIP == 1)THEN
  NPL_tmp = NPL-1
ELSE
  NPL_tmp = NPL
ENDIF
CALL QSORTC2D(obj%voro(1:obj%k,:),obj%voroidx(1:obj%k,:))

!!
!! Use local variable for voro and voroidx to allow easy change
!! from perturbation value to perturbation+vref
!!
voro    = 0._RP
voroidx = 0
voro    = obj%voro
voroidx = obj%voroidx

!!
!! Apply reference profile:
!!
IF(I_VREF == 1)THEN
  DO ivo = 1,obj%k
    CALL GETREF(vref,vpvsref,voro(ivo,1))
    !WRITE(*,201)ivo,voro(ivo,1),voro(ivo,2),voro(ivo,3)
    IF(voroidx(ivo,2) == 1)voro(ivo,2) = vref + voro(ivo,2)
    IF(voroidx(ivo,3) == 1)voro(ivo,3) = vpvsref + voro(ivo,3)
    !WRITE(*,201)ivo,voro(ivo,1),voro(ivo,2),voro(ivo,3),vref,vpvsref
    !WRITE(*,*)''
  ENDDO
ENDIF
201 FORMAT(I,5F12.6)
!!
!!  Find interfaces for each parameter
!!
niface = 0._RP
DO ipar = 2,NPL
  niface(ipar-1) = SUM(voroidx(:,ipar))-1
ENDDO
obj%ziface = 0._RP
voroz  = 0._RP
voroh  = 0._RP
iface  = 0
DO ipar = 2,NPL
  vorotmp = 0._RP
  itmp = 0
  DO ivo = 1,obj%k
    IF(voroidx(ivo,ipar) == 1)THEN
      itmp = itmp + 1
      vorotmp(itmp,:) = (/voro(ivo,1),voro(ivo,ipar)/)
    ENDIF
  ENDDO
  IF(niface(ipar-1) > 0)THEN
    DO ivo = 1,niface(ipar-1)
      iface = iface + 1
      obj%ziface(iface) = vorotmp(ivo+1,1)
      voroz(ivo,ipar-1,1) = obj%ziface(iface)
      voroz(ivo,ipar-1,2) = vorotmp(ivo,2)
    ENDDO
    voroz(itmp,ipar-1,2) = vorotmp(itmp,2)
  ELSE
    voroz(1,ipar-1,2) = vorotmp(1,2)
  ENDIF
ENDDO

ntot = SUM(niface)
CALL QSORTC1D(obj%ziface(1:ntot))

!!
!!  Find unique interfaces
!!
obj%nunique = obj%k-1
IF(ntot > 1)THEN
  DO ivo = 1,ntot
    itmp = 0
    ifaceidx = 0
    DO ivo2 = ivo+1,ntot
      IF(obj%ziface(ivo) /= obj%ziface(ivo2))THEN
        itmp = itmp + 1
        ifaceidx(itmp) = ivo2
      ENDIF
    ENDDO
    IF(itmp > 0)THEN
      obj%ziface(ivo+1:ivo+itmp) = obj%ziface(ifaceidx(1:itmp))
      obj%ziface(ivo+itmp+1:) = 0._RP
    ENDIF
  ENDDO
ENDIF

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%nunique
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%nunique,1) = obj%ziface(1:obj%nunique)
DO ipar = 2,NPL_tmp
  ivo = 1
  IF(niface(ipar-1) /= 0)THEN
    DO ilay = 1,obj%nunique
      IF(obj%ziface(ilay) <= voroz(ivo,ipar-1,1))THEN
        partmp(ilay,ipar) = voroz(ivo,ipar-1,2)
      ELSE
        ivo = ivo + 1
        partmp(ilay,ipar) = voroz(ivo,ipar-1,2)
      ENDIF
      IF(obj%ziface(ilay) >= voroz(niface(ipar-1),ipar-1,1))THEN
        partmp(ilay+1:obj%nunique+1,ipar) = voroz(ivo+1,ipar-1,2)
        EXIT
      ENDIF
    ENDDO
  ELSE
    partmp(1:obj%nunique+1,ipar) = voroz(1,ipar-1,2)
  ENDIF
ENDDO
partmp(1:obj%nunique,1) = obj%hiface(1:obj%nunique)
IF(IDIP == 1)THEN
  ipar = NPL
  ivo = 1
  partmp(:,ipar) = 0._RP
  IF(SUM(voroidx(:,ipar)) > 0)THEN
    DO ilay = 2,obj%nunique+1
      IF(voroidx(ilay,ipar) == 1)THEN
        partmp(ilay,ipar) = voro(ilay,ipar)
      ENDIF
    ENDDO
  ENDIF
ENDIF

obj%par = 0._RP
DO ilay = 1,obj%nunique
  obj%par((ilay-1)*NPL+1:ilay*NPL) = partmp(ilay,:)
ENDDO
obj%par(obj%nunique*NPL+1:((obj%nunique+1)*NPL)-1) = partmp(ilay,2:)

!IF(obj%voro(1,1)<0._RP)THEN
!PRINT*,''
!PRINT*,'After:',rank
!DO ilay = 1,obj%nunique+1
!  WRITE(*,203)ilay,partmp(ilay,:)
!ENDDO
!CALL PRINTPAR(obj)
!PRINT*,''
!PRINT*,obj%par
!PRINT*,''
!STOP
!PRINT*,'ntot',ntot
!PRINT*,'unique',obj%nunique
!PRINT*,'niface',niface
!PRINT*,'ziface',obj%ziface(1:ntot)
!DO ilay = 1,obj%nunique+1
!  WRITE(*,203)ilay,obj%par((ilay-1)*NPL+1:ilay*NPL)
!ENDDO
!PRINT*,''
!PRINT*,''
!  PRINT*,'INTERP'
!  STOP
!ENDIF
!
!STOP
!200 FORMAT(2I4,20F8.2)
!201 FORMAT(a,20F8.2)
!202 FORMAT(a,20I4)
203 FORMAT(I4,20F8.2)

END SUBROUTINE INTERPLAYER
!=======================================================================

SUBROUTINE GETREF(vref,vpvsref,z)
!!==============================================================================
!!
!! Reads observed data.
!!
USE RJMCMC_COM
IMPLICIT NONE
REAL(KIND=RP)   :: z,vref,vpvsref,grad,dz
INTEGER(KIND=IB):: ipar,iint
IF(z >= vel_ref(1,NVELREF))THEN
  vref = vel_ref(2,NVELREF)
  vpvsref = vel_ref(3,NVELREF)
ELSEIF(z == 0.)THEN
  vref = vel_ref(2,1)
  vpvsref = vel_ref(3,1)
ELSE
  iint = 0
  DO ipar=1,NVELREF
    IF((z-vel_ref(1,ipar)) <= 0.)THEN
      EXIT
    ENDIF
    iint = iint + 1
  ENDDO
  grad = (vel_ref(2,iint+1)-vel_ref(2,iint))/(vel_ref(1,iint+1)-vel_ref(1,iint))
  dz = (z-vel_ref(1,iint))
  vref = vel_ref(2,iint) + dz*grad
  grad = (vel_ref(3,iint+1)-vel_ref(3,iint))/(vel_ref(1,iint+1)-vel_ref(1,iint))
  vpvsref = vel_ref(3,iint) + dz*grad
ENDIF
RETURN
END SUBROUTINE GETREF
!!=======================================================================

SUBROUTINE ARPRED_RF(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxar(ifr) == 1)THEN
   k = 1
   obj%DarH(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   dres1 = obj%DresH(ifr,idata:idatb)

   dar1(1)=0._RP          ! Matlab sets first point to zero...
   DO i=2,idatb-idata+1
      dar1(i) = 0
      IF(k >= i)THEN
         DO j=1,i-1
            dar1(i) = dar1(i) + obj%arpar((ifr-1)+j) * dres1(i-j)
         ENDDO
      ELSE
         DO j=1,k
            dar1(i) = dar1(i) + obj%arpar((ifr-1)+j) * dres1(i-j)
         ENDDO
      ENDIF
   ENDDO
   obj%DarH(ifr,idata:idatb) = dar1
   obj%DarH(ifr,idata) = 0._RP
   obj%DarH(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_RF
!!=======================================================================

SUBROUTINE ARPRED_SWD(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxarSWD(ifr) == 1)THEN
   k = 1
   obj%DarSWD(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   dres1 = obj%DresSWD(ifr,idata:idatb)

   dar1(1)=0._RP          ! Matlab sets first point to zero...
   DO i=2,idatb-idata+1
      dar1(i) = 0
      IF(k >= i)THEN
         DO j=1,i-1
            dar1(i) = dar1(i) + obj%arparSWD((ifr-1)+j) * dres1(i-j)
         ENDDO
      ELSE
         DO j=1,k
            dar1(i) = dar1(i) + obj%arparSWD((ifr-1)+j) * dres1(i-j)
         ENDDO
      ENDIF
   ENDDO
   obj%DarSWD(ifr,idata:idatb) = dar1
   obj%DarSWD(ifr,idata) = 0._RP
   obj%DarSWD(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_SWD
!!=======================================================================

SUBROUTINE CHECKBOUNDS_ARMXRF(obj,ibadlogL)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr,ibadlogL
TYPE(objstruc):: obj

DO ifr = 1,NRF1
   IF(MAXVAL(obj%DarH(ifr,:)) > armxH)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
   IF(MINVAL(obj%DarH(ifr,:)) < -armxH)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
ENDDO

RETURN
END SUBROUTINE CHECKBOUNDS_ARMXRF
!!=======================================================================

SUBROUTINE CHECKBOUNDS_ARMXSWD(obj,ibadlogL)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr,ibadlogL
TYPE(objstruc):: obj

DO ifr = 1,NMODE
   IF(MAXVAL(obj%DarSWD(ifr,:)) > armxSWD)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
   IF(MINVAL(obj%DarSWD(ifr,:)) < -armxSWD)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
ENDDO

RETURN
END SUBROUTINE CHECKBOUNDS_ARMXSWD
!=======================================================================

SUBROUTINE LOGLHOOD2(obj)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj

!!
!!  Compute log likelihood
!!
obj%logL = 1._RP

RETURN
END SUBROUTINE LOGLHOOD2
!!=======================================================================
!!EOF
