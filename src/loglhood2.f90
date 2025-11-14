 ! DO ivo = 1,obj%k
  !  CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
   ! partmp(ivo,2) = vref + obj%voro(ivo,2)
   ! partmp(ivo,3) = vpvsref + obj%voro(ivo,3)
   ! !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
 ! ENDDO 
!!=======================================================================

SUBROUTINE LOGLHOOD(obj,ipred)
!!=======================================================================
!! This is LOGLHOOD function that combines various data sets for 
!! joint inversion. 
USE RJMCMC_COM
IMPLICIT NONE
INCLUDE 'raysum/params.h'
INTEGER(KIND=IB):: ipred
TYPE (objstruc) :: obj
REAL(KIND=RP)   :: logL_RV, logL_SWD, logL_MT
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod

!!
!! Generate input vector for SWD
!!
!! curmod = (/ thick, rho, alph, beta, %P, %S, tr, pl, st, di /)
IF(I_RV == 1 .OR. I_RV == -1 .OR. I_SWD == 1)THEN
  curmod = 0._SP
  IF(iraysum == 1)THEN
    CALL MAKE_CURMOD(obj,curmod)
  ELSE
    !! For ray3d or Shibutani code
    CALL MAKE_CURMOD_ray3d(obj,curmod(:,1:7))
  ENDIF
ENDIF
IF(I_RV == 1)THEN
  !! log likelihood for direct seimogram inversion (radial & vertical data)
  IF(iraysum == 1)THEN
    !! inversion via raysum
    CALL LOGLHOOD_RV(obj,curmod,ipred,logL_RV)
  ELSE
    !! inversion via ray 3D
    CALL LOGLHOOD_RV_ray3d(obj,ipred,logL_RV)
  ENDIF
ELSEIF(I_RV == 2)THEN
  !! log likelihood for direct seimogram inversion (radial & vertical data)
  !! inversion via raysum
  CALL LOGLHOOD_RV_CROSS(obj,curmod,ipred,logL_RV)
ELSEIF(I_RV == -1)THEN
  !! Receiver function inversion
  IF(iraysum == 1)THEN
    !! inversion via raysum
    CALL LOGLHOOD_RF(obj,curmod,ipred,logL_RV)
  ELSEIF(iraysum == 0)THEN
    !! inversion via ray 3D
    CALL LOGLHOOD_RF_ray3d(obj,ipred,logL_RV)
  ELSEIF(iraysum == -1)THEN
    !! inversion via Shibutani implementation
    CALL LOGLHOOD_RF_shibutani(obj,ipred,logL_RV)
  ENDIF
ELSE
  logL_RV = 0._RP
ENDIF
!! log likelihood SWD data
IF(I_SWD == 1)THEN
  CALL LOGLHOOD_SWD(obj,curmod,ipred,logL_SWD)
ELSE
  logL_SWD = 0._RP
ENDIF

!! log likelihood MT data
IF(I_MT == 1)THEN
  CALL LOGLHOOD_MT(obj,ipred,logL_MT)
ELSE
  logL_MT = 0._RP
ENDIF

!!
!! Joint likelihood (assumes independent errors on the vasious data sets)
!!
obj%logL = logL_RV + logL_SWD + logL_MT
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
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,ipred,iaz
INTEGER(KIND=IB)                 :: ierr_neg,ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpR,EtmpV,EtmpT
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
!! Forward specific variables
REAL(KIND = SP)                  :: aa(3,3,3,3,maxlay),ar_list(3,3,3,3,maxlay,2),rot(3,3,maxlay)
REAL(KIND = SP)                  :: amp_in,tt(maxph,maxtr),amp(3,maxph,maxtr)
LOGICAL :: bailout,errflag
CHARACTER(len=64) :: phname
LOGICAL :: ISNAN
REAL(KIND=RP) :: DpredRG(NRF1,NTIME),DpredVG(NRF1,NTIME),DpredTG(NRF1,NTIME)
REAL(KIND=RP) :: DpredR(NTIME2),DpredV(NTIME2),DpredT(NTIME2)
REAL(KIND=RP) :: RRG(NRRG),VVG(NRRG),RGRG(NRGRG),VGVG(NRGRG),TTG(NRRG),TGTG(NRGRG)
REAL(KIND=RP) :: S(NSRC)


!IF(IMAP == 1)THEN
!  PRINT*,'CURMOD IN raysum'
!  DO ilay=1,obj%nunique+1
!    WRITE(*,206)ilay,curmod(ilay,1:10)
!  ENDDO
!  206   FORMAT(I3,10F12.4)
!ENDIF

!!
!!  Build a_ijkl tensors
!!
isoflag(1:obj%nunique+1) = .TRUE.
call buildmodel(aa,ar_list,rot,curmod(:,1),curmod(:,2),&
     curmod(:,3),curmod(:,4),isoflag,curmod(:,5),curmod(:,6),&
     curmod(:,7),curmod(:,8),curmod(:,9),curmod(:,10),obj%nunique+1)

!!
!! Compute multiples
!!
!print*,'width',width
numph = 0
iphase = 1
!! no multiples
if (mults .eq. 0) then
  call ph_direct(phaselist,nseg,numph,obj%nunique+1,iphase)
end if
!! Only multiple of last interface (can be Moho but not necessarily...)
if (mults .eq. 1) then
  ilay = obj%nunique
  call ph_direct_rf(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
  call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
end if
!! Josip's version of multiple (all significant mults)
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
!IF(ierr_neg == 1) raysumfail = raysumfail + 1
IF(ierr_neg == 1)THEN
  PRINT*,'ERROR'
ENDIF
!! Normalize amplitudes by direct P
call norm_arrivals(amp,baz2,slow2,curmod(1,3),curmod(1,4),&
                   curmod(1,2),ntr,numph,1,1)
!IF(IMAP == 1)THEN
!  !! Write out arrivals
!  call writearrivals(iounit1,tt,amp,ntr,numph)
!  close(unit=iounit1)
!ENDIF
!!
!! Make traces (in NEZ) and rotate into RTV (out_rot = 1)
!!
call make_traces(tt,amp,ntr,numph,NTIME,sampling_dt,width2,align,shift2,synth_cart)
!if (out_rot .eq. 1) then
!call writetraces(iounit2,synth_cart,ntr,NTIME,sampling_dt,align,shift2)
!close(unit=iounit2)

call rot_traces(synth_cart,baz2,ntr,NTIME,synth_ph)
!call writetraces(iounit2,synth_ph,ntr,NTIME,sampling_dt,align,shift2)
!close(unit=iounit2)

!end if
!if (out_rot .eq. 2) then
!  call fs_traces(synth_cart,baz2,slow2,curmod(:,3),curmod(:,4),&
!       curmod(:,2),ntr,NTIME,synth_ph)
!end if

DO iaz = 1,NRF1
  DpredRG(iaz,1:NTIME)  = REAL(synth_ph(1,1:NTIME,iaz),RP)
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
    WRITE(50,201) DpredRG(iaz,:)
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
!  open(unit=iounit1,file='sample.ph',status='unknown')
!  call writephases(iounit1,phaselist,nseg,numph)
!  close(unit=iounit1)
!  write(*,*) 'Phases written to sample.ph'
  201 FORMAT(1024ES18.10)
ENDIF

DO iaz = 1,NRF1
  !! Time domain convolution for source equation terms:
  !! SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
  CALL CONVT(obj%DobsR(iaz,:),DpredRG(iaz,:),RRG,NTIME2,NTIME,NRRG)
  CALL CONVT(obj%DobsV(iaz,:),DpredVG(iaz,:),VVG,NTIME2,NTIME,NRRG)
  CALL CONVT(DpredRG(iaz,:),DpredRG(iaz,:),RGRG,NTIME,NTIME,NRGRG)
  CALL CONVT(DpredVG(iaz,:),DpredVG(iaz,:),VGVG,NTIME,NTIME,NRGRG)
  IF(I_T == 1)THEN
    CALL CONVT(obj%DobsT(iaz,:),DpredTG(iaz,:),TTG,NTIME2,NTIME,NRRG)
    CALL CONVT(DpredTG(iaz,:),DpredTG(iaz,:),TGTG,NTIME,NTIME,NRGRG)
  ENDIF

  !! Estimate source by time domain dconvolution 
  !! (accounting for different std dev on H and V):
  IF(ICOV == 1)THEN
    IF(I_T == 0)THEN
      CALL DCONVT(RRG/obj%sdparR(iaz)**2+VVG/obj%sdparV(iaz)**2, &
      RGRG/obj%sdparR(iaz)**2+VGVG/obj%sdparV(iaz)**2,S,NRRG,NRGRG,NSRC)
    ELSEIF(I_T == 1)THEN
      CALL DCONVT(RRG/obj%sdparR(iaz)**2+VVG/obj%sdparV(iaz)**2+TTG/obj%sdparT(iaz)**2, &
      RGRG/obj%sdparR(iaz)**2+VGVG/obj%sdparV(iaz)**2+TGTG/obj%sdparT(iaz)**2,S,NRRG,NRGRG,NSRC)
    ENDIF
  ELSE
    PRINT*,'WARNING ICOV problem'
    CALL DCONVT(RRG+VVG,RGRG+VGVG,S,NRRG,NRGRG,NSRC)
  ENDIF
  !! Data predictions via time domain convolution:
  CALL CONVT(S,DpredRG(iaz,:),DpredR,NSRC,NTIME,NTIME2)
  CALL CONVT(S,DpredVG(iaz,:),DpredV,NSRC,NTIME,NTIME2)
  IF(I_T == 1)CALL CONVT(S,DpredTG(iaz,:),DpredT,NSRC,NTIME,NTIME2)

  obj%S(iaz,1:NSRC) = S(1:NSRC)
  obj%DpredR(iaz,1:NTIME) = DpredR(1:NTIME)
  obj%DpredV(iaz,1:NTIME) = DpredV(1:NTIME)
  obj%DresR(iaz,1:NTIME) = obj%DobsR(iaz,1:NTIME)-obj%DpredR(iaz,1:NTIME)
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
  EtmpR = 0._RP
  EtmpV = 0._RP
  EtmpT = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpR(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresR(iaz,:)**2._RP)/(2._RP*obj%sdparR(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparR(iaz)))
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
  logL = SUM(EtmpR) + SUM(EtmpV)
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
!!=======================================================================

SUBROUTINE LOGLHOOD_RV_CROSS(obj,curmod,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted R and V components
!!
USE RJMCMC_COM
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,ipred,iaz
INTEGER(KIND=IB)                 :: ierr_neg,ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpR,EtmpV,EtmpT
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
!! Forward specific variables
REAL(KIND = SP)                  :: aa(3,3,3,3,maxlay),ar_list(3,3,3,3,maxlay,2),rot(3,3,maxlay)
REAL(KIND = SP)                  :: amp_in,tt(maxph,maxtr),amp(3,maxph,maxtr)
LOGICAL :: bailout,errflag
CHARACTER(len=64) :: phname
LOGICAL :: ISNAN
REAL(KIND=RP) :: DpredRG(NRF1,NTIME),DpredVG(NRF1,NTIME),DpredTG(NRF1,NTIME)
REAL(KIND=RP) :: DpredR(NTIME2),DpredV(NTIME2),DpredT(NTIME2)
REAL(KIND=RP) :: VRG(NRRG),RVG(NRRG)


!IF(IMAP == 1)THEN
!  PRINT*,'CURMOD IN raysum'
!  DO ilay=1,obj%nunique+1
!    WRITE(*,206)ilay,curmod(ilay,1:10)
!  ENDDO
!  206   FORMAT(I3,10F12.4)
!ENDIF

!!
!!  Build a_ijkl tensors
!!
isoflag(1:obj%nunique+1) = .TRUE.
call buildmodel(aa,ar_list,rot,curmod(:,1),curmod(:,2),&
     curmod(:,3),curmod(:,4),isoflag,curmod(:,5),curmod(:,6),&
     curmod(:,7),curmod(:,8),curmod(:,9),curmod(:,10),obj%nunique+1)

!!
!! Compute multiples
!!
!print*,'width',width
numph = 0
iphase = 1
!! no multiples
if (mults .eq. 0) then
  call ph_direct(phaselist,nseg,numph,obj%nunique+1,iphase)
end if
!! Only multiple of last interface (can be Moho but not necessarily...)
if (mults .eq. 1) then
  ilay = obj%nunique
  call ph_direct_rf(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
  call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
end if
!! Josip's version of multiple (all significant mults)
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
!IF(ierr_neg == 1) raysumfail = raysumfail + 1
IF(ierr_neg == 1)THEN
  PRINT*,'ERROR'
ENDIF
!! Normalize amplitudes by direct P
call norm_arrivals(amp,baz2,slow2,curmod(1,3),curmod(1,4),&
                   curmod(1,2),ntr,numph,1,1)
!IF(IMAP == 1)THEN
!  !! Write out arrivals
!  call writearrivals(iounit1,tt,amp,ntr,numph)
!  close(unit=iounit1)
!ENDIF
!!
!! Make traces (in NEZ) and rotate into RTV (out_rot = 1)
!!
call make_traces(tt,amp,ntr,numph,NTIME,sampling_dt,width2,align,shift2,synth_cart)
!if (out_rot .eq. 1) then
!call writetraces(iounit2,synth_cart,ntr,NTIME,sampling_dt,align,shift2)
!close(unit=iounit2)

call rot_traces(synth_cart,baz2,ntr,NTIME,synth_ph)
!call writetraces(iounit2,synth_ph,ntr,NTIME,sampling_dt,align,shift2)
!close(unit=iounit2)

!end if
!if (out_rot .eq. 2) then
!  call fs_traces(synth_cart,baz2,slow2,curmod(:,3),curmod(:,4),&
!       curmod(:,2),ntr,NTIME,synth_ph)
!end if

!! GF predictions:
!PRINT*,'synth',synth_ph(1,1:NTIME,iaz)

DO iaz = 1,NRF1
  DpredRG(iaz,1:NTIME)  = REAL(synth_ph(1,1:NTIME,iaz),RP)
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
    WRITE(50,201) DpredRG(iaz,:)
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
!  open(unit=iounit1,file='sample.ph',status='unknown')
!  call writephases(iounit1,phaselist,nseg,numph)
!  close(unit=iounit1)
!  write(*,*) 'Phases written to sample.ph'
  201 FORMAT(1024ES18.10)
ENDIF

DO iaz = 1,NRF1
  !! Time domain convolution for source equation terms:
  !! SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
  CALL CONVT(DpredVG(iaz,:),obj%DobsR(iaz,:),VRG,NTIME,NTIME2,NRRG)
  CALL CONVT(DpredRG(iaz,:),obj%DobsV(iaz,:),RVG,NTIME,NTIME2,NRRG)
  obj%DresR(iaz,1:NTIME) = 0.
  obj%DresR(iaz,1:NTIME) = VRG(1:NTIME)-RVG(1:NTIME)
  obj%DpredR(iaz,1:NTIME) = VRG(1:NTIME)
  obj%DpredV(iaz,1:NTIME) = RVG(1:NTIME)
ENDDO
IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='vRG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=51,FILE='rVG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  WRITE(50,201) VRG
  WRITE(51,201) RVG
  CLOSE(50)
  CLOSE(51)
!  open(unit=iounit1,file='sample.ph',status='unknown')
!  call writephases(iounit1,phaselist,nseg,numph)
!  close(unit=iounit1)
!  write(*,*) 'Phases written to sample.ph'
ENDIF

ibadlogL = 0

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpR = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpR(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresR(iaz,:)**2._RP)/(2._RP*obj%sdparR(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparR(iaz)))
    ENDDO
  ENDIF
  logL = SUM(EtmpR)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
  logL = -HUGE(1._RP)
  ibadlogL = 0
ENDIF

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_RV_CROSS
!!=======================================================================

SUBROUTINE LOGLHOOD_RV_ray3d(obj,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted R and V components
!!
USE RJMCMC_COM
USE ieee_arithmetic
USE RAY3D_COM
IMPLICIT NONE
!include 'ray3d/na_param.inc'
!include 'ray3d/ray3d_param.inc'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,iparcur,ipred,iaz
INTEGER(KIND=IB)                 :: ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpR,EtmpV,EtmpT
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
INTEGER(KIND=IB):: ifail
LOGICAL:: ISNAN
REAL(KIND=RP):: DpredRG(NRF1,NTIME),DpredTG(NRF1,NTIME),DpredVG(NRF1,NTIME)
REAL(KIND=RP):: DpredR(NTIME2),DpredT(NTIME2),DpredV(NTIME2)
REAL(KIND=RP):: RRG(NRRG),TTG(NRRG),VVG(NRRG),RGRG(NRGRG),TGTG(NRGRG),VGVG(NRGRG)
REAL(KIND=RP):: S(NSRC)

! ray3d specific variables
INTEGER(KIND=IB):: nly,jlay,iqq,npts
REAL(KIND=SP)   :: hsub,hsubinc,bazs,ps
real spike(maxsamp,3),rad(maxsamp),tra(maxsamp),ver(maxsamp), &
     rfrad(maxsamp),rftra(maxsamp)

ifail = 0  ! Added by SED to signal failure of forward model

!! zero traces

iqq = 0
DO iaz = 1,NRF1
! calculate seismograms using ay3d for each trace

  spike  = 0.
  bazs = baz2(iaz)   !! baz is in deg. (for ray3d), baz2 is in radians (for raysum)
  ps   = REAL(slow2(iaz),SP)
  npts = NTIME
  dt = REAL(sampling_dt,SP)
  t_af = ((npts-1)*dt)-(shift2-dt)

  call ray3dn(strike,dip,h,alpha,beta,rho,bazs,ps,dt, &
       t_af+shift2-dt,shift2-dt,obj%nunique+1,npts,spike,iqq,ifail)

  IF (ifail == 1) THEN   ! Added by SED
    obj%logL = -HUGE(1._RP)
    RETURN
  ENDIF
  IF(NTIME /= npts) PRINT*,'ERROR: seismogram record length mismatch',NTIME,npts
  !! GF predictions:
  DpredRG(iaz,1:NTIME)  = spike(1:NTIME,1)   ! Jan Dettmer: these are impulse respones
  DpredTG(iaz,1:NTIME)  = spike(1:NTIME,2)
  DpredVG(iaz,1:NTIME)  = spike(1:NTIME,3)
ENDDO

IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='DpredRG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=51,FILE='DpredVG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=52,FILE='DpredTG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  DO iaz = 1,NRF1
    WRITE(50,201) DpredRG(iaz,1:NTIME)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(51,201) DpredVG(iaz,1:NTIME)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(52,201) DpredTG(iaz,1:NTIME)
  ENDDO
  CLOSE(50)
  CLOSE(51)
  CLOSE(52)
  201 FORMAT(1024ES18.10)
ENDIF

DO iaz = 1,NRF1
  !! Time domain convolution for source equation terms:
  !! SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
  CALL CONVT(obj%DobsR(iaz,:),DpredRG(iaz,:),RRG,NTIME2,NTIME,NRRG)
  CALL CONVT(obj%DobsV(iaz,:),DpredVG(iaz,:),VVG,NTIME2,NTIME,NRRG)
  CALL CONVT(DpredRG(iaz,:),DpredRG(iaz,:),RGRG,NTIME,NTIME,NRGRG)
  CALL CONVT(DpredVG(iaz,:),DpredVG(iaz,:),VGVG,NTIME,NTIME,NRGRG)
  IF(I_T == 1)THEN
    CALL CONVT(obj%DobsT(iaz,:),DpredTG(iaz,:),TTG,NTIME2,NTIME,NRRG)
    CALL CONVT(DpredTG(iaz,:),DpredTG(iaz,:),TGTG,NTIME,NTIME,NRGRG)
  ENDIF

  !! Estimate source by time domain dconvolution 
  !! (accounting for different std dev on H and V):
  IF(ICOV == 1)THEN
    IF(I_T == 0)THEN
      CALL DCONVT(RRG/obj%sdparR(iaz)**2+VVG/obj%sdparV(iaz)**2, &
      RGRG/obj%sdparR(iaz)**2+VGVG/obj%sdparV(iaz)**2,S,NRRG,NRGRG,NSRC)
    ELSEIF(I_T == 1)THEN
      CALL DCONVT(RRG/obj%sdparR(iaz)**2+VVG/obj%sdparV(iaz)**2+TTG/obj%sdparT(iaz)**2, &
      RGRG/obj%sdparR(iaz)**2+VGVG/obj%sdparV(iaz)**2+TGTG/obj%sdparT(iaz)**2,S,NRRG,NRGRG,NSRC)
    ENDIF
  ELSE
    PRINT*,'WARNING ICOV problem'
    CALL DCONVT(RRG+VVG,RGRG+VGVG,S,NRRG,NRGRG,NSRC)
  ENDIF
  !! Data predictions via time domain convolution:
  CALL CONVT(S,DpredRG(iaz,:),DpredR,NSRC,NTIME,NTIME2)
  CALL CONVT(S,DpredVG(iaz,:),DpredV,NSRC,NTIME,NTIME2)
  IF(I_T == 1)CALL CONVT(S,DpredTG(iaz,:),DpredT,NSRC,NTIME,NTIME2)

  obj%S(iaz,1:NSRC) = S(1:NSRC)
  obj%DpredR(iaz,1:NTIME) = DpredR(1:NTIME)
  obj%DpredV(iaz,1:NTIME) = DpredV(1:NTIME)
  obj%DresR(iaz,1:NTIME) = obj%DobsR(iaz,1:NTIME)-obj%DpredR(iaz,1:NTIME)
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
  EtmpR = 0._RP
  EtmpV = 0._RP
  EtmpT = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpR(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresR(iaz,:)**2._RP)/(2._RP*obj%sdparR(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparR(iaz)))
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
  logL = SUM(EtmpR) + SUM(EtmpV)
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
END SUBROUTINE LOGLHOOD_RV_ray3d
!=======================================================================

SUBROUTINE LOGLHOOD_RF(obj,curmod,ipred,logL)
!=======================================================================
USE RJMCMC_COM
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,iparcur,ipred,iaz
INTEGER(KIND=IB)                 :: ierr_neg,ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10) :: curmod
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpR
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
REAL(KIND=RP), DIMENSION(NRF1,NTIME) :: DpredR
!!!!!Pejman: for creating synthtic impulse responces**************************
!REAL(KIND=RP), ALLOCATABLE, DIMENSION(:,:) :: DpredRG
!REAL(KIND=RP), ALLOCATABLE, DIMENSION(:,:) :: DpredVG
!ALLOCATE(DpredRG(NRF1, NTIME))
!ALLOCATE(DpredVG(NRF1, NTIME))
REAL(KIND=RP), DIMENSION(NRF1,NTIME) :: DpredRG
REAL(KIND=RP), DIMENSION(NRF1,NTIME) :: DpredVG
!****************************************************************************

!!
!!  Compute predicted RF
!!
!!  Build a_ijkl tensors
!!
isoflag(1:obj%nunique+1) = .TRUE.
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
  !PRINT*, 'direct_phases'
  !Write(*,*) phaselist
  do ilay=1,obj%nunique
    call ph_rfmults(phaselist,nseg,numph,obj%nunique+1,ilay,iphase)
  end do
  !PRINT*, 'multiples'
  !Write(*,*) phaselist
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
call make_traces(tt,amp,ntr,numph,NTIME,sampling_dt,width2,align,shift2,synth_cart)
!if (out_rot .eq. 1) then
call rot_traces(synth_cart,baz2,ntr,NTIME,synth_ph)
!end if
!if (out_rot .eq. 2) then
!  call fs_traces(synth_cart,baz2,slow2,curmod(1,3),curmod(1,4),&
!       curmod(1,2),ntr,NTIME,synth_ph)
!end if

!! GF predictions:
DpredR(1,1:NTIME)  = REAL(synth_ph(1,1:NTIME,1),RP)/MAXVAL(synth_ph(3,1:NTIME,1))
IF (ITAPER == 1) DpredR(1,1:NTIME) = DpredR(1,1:NTIME) * taper_dpred(1:NTIME)
obj%DpredR(1,1:NTIME)  = DpredR(1,1:NTIME)
obj%DresR(1,1:NTIME) = obj%DobsR(1,1:NTIME)-obj%DpredR(1,1:NTIME)
IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='DpredRF.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  WRITE(50,201) (obj%DpredR(1,ifr),ifr=1,NTIME)
  CLOSE(50)
  201 FORMAT(ES18.10)
ENDIF
!!!!! Pejman: for creating synthetic impulse responces *****************************
IF(IMAP == 1)THEN
  DpredRG(1,1:NTIME) = REAL(synth_ph(1,1:NTIME,1), RP)
  OPEN(UNIT=100,FILE='DpredRG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  WRITE(100,2001) (DpredRG(1,ifr),ifr=1,NTIME)
  CLOSE(100)

  DpredVG(1,1:NTIME) = REAL(synth_ph(3,1:NTIME,1), RP)
  OPEN(UNIT=200,FILE='DpredVG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  WRITE(200,2001) (DpredVG(1,ifr),ifr=1,NTIME)
  CLOSE(200)
  2001 FORMAT(ES18.10)
ENDIF
!************************************************************************************
ibadlogL = 0
IF(IAR == 1)THEN
   obj%DarR  = 0._RP
   !!
   !!  Compute autoregressive model
   !!
   CALL ARPRED_RF(obj,1,1,NTIME)
   !! Recompute predicted data as ith autoregressive model
   obj%DresR = obj%DresR-obj%DarR

   !! Check if predicted AR model data are outside max allowed bounds
   CALL CHECKBOUNDS_ARMXRF(obj,ibadlogL)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpR = 0._RP
  IF(ICOV==0)THEN
  !!
  !! implicitly sample over sigma
  !!
  DO iaz = 1,NRF1
    EtmpR(iaz) = -REAL(NTIME,RP)/2._RP * LOG( SUM(obj%DresR(iaz,:)**2._RP) / REAL(NTIME,RP) )
  ENDDO

  ELSEIF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpR(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresR(iaz,:)**2._RP)/(2._RP*obj%sdparR(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparR(iaz)))
    ENDDO
!  ENDIF
   ELSEIF(ICOV == 2)THEN
   !! Empricial Cov. Mat estimate and magniutde scaling  scaling: Cd = sC'd not Cd = (s^2)C'd
     DO iaz =1,NRF1
       EtmpR(iaz) = -DOT_PRODUCT(MATMUL(obj%DresR(iaz,:),Cdi),obj%DresR(iaz,:))/(2._RP*obj%sdparR(iaz)) &
                    -REAL(NTIME,RP)/2._RP*LOG(obj%sdparR(iaz))
     ENDDO
   ENDIF
 !  ELSEIF(ICOV == 2)THEN
       !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
       !! dev parameter... (see Dettmer et al. 2014 GJI)
   !  EtmpR(iaz) = -DOT_PRODUCT(MATMUL(obj%DresR(iaz,:),Cdi),obj%DresR(iaz,:))/(2._RP*obj%sdparR(iaz)) &
  !                -REAL(NTIME,RP)/2._RP*LOG(obj%sdparR(iaz))
 !    IF(I_T == 1)THEN
     !  EtmpT(iaz) = -DOT_PRODUCT(MATMUL(obj%DresT(iaz,:),Cdi),obj%DresT(iaz,:))/(2._RP*obj%sdparT(iaz)) &
    !                  -REAL(NTIME,RP)/2._RP*LOG(obj%sdparT(iaz))
   !  ENDIF
  ! ENDIF
 ! ENDDO

  logL = SUM(EtmpR)
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

SUBROUTINE LOGLHOOD_RF_ray3d(obj,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted R and V components
!!
USE RJMCMC_COM
USE ieee_arithmetic
USE RAY3D_COM
IMPLICIT NONE
!include 'ray3d/na_param.inc'
!include 'ray3d/ray3d_param.inc'

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,iparcur,ipred,iaz
INTEGER(KIND=IB)                 :: ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpR,EtmpV,EtmpT
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
!! Forward specific variables
!REAL(KIND=SP)                  :: aa(3,3,3,3,NLMX),ar_list(3,3,3,3,NLMX,2),rot(3,3,NLMX)
!REAL(KIND=SP)                  :: amp_in,tt(maxph,maxtr),amp(3,maxph,maxtr)
!LOGICAL :: bailout,errflag
!CHARACTER(len=64) :: phname
INTEGER(KIND=IB):: ifail
LOGICAL:: ISNAN
REAL(KIND=RP):: DpredRG(NRF1,NTIME),DpredTG(NRF1,NTIME),DpredVG(NRF1,NTIME)
REAL(KIND=RP):: DpredR(NTIME2),DpredT(NTIME2),DpredV(NTIME2)
REAL(KIND=RP):: RRG(NRRG),TTG(NRRG),VVG(NRRG),RGRG(NRGRG),TGTG(NRGRG),VGVG(NRGRG)
REAL(KIND=RP):: S(NSRC)

! ray3d specific variables
INTEGER(KIND=IB):: nly,jlay,iqq,npts
REAL(KIND=SP)   :: hsub,hsubinc,bazs,ps
real spike(maxsamp,3),rad(maxsamp),tra(maxsamp),ver(maxsamp), &
     rfrad(maxsamp),rftra(maxsamp)

ifail = 0  ! Added by SED to signal failure of forward model

!! zero traces

iqq = 0
DO iaz = 1,NRF1
! calculate seismograms using ray3d for each trace

  spike  = 0.
  bazs = baz2(iaz)   !! baz is in deg. (for ray3d), baz2 is in radians (for raysum)
  ps   = REAL(slow2(iaz),SP) !! ps is for km/s units here!
  npts = NTIME
  dt = REAL(sampling_dt,SP)
  t_af = ((npts-1)*dt)-(shift2-dt)

  call ray3dn(strike,dip,h,alpha,beta,rho,bazs,ps,dt, &
       t_af+shift2-dt,shift2-dt,obj%nunique+1,npts,spike,iqq,ifail)

  IF (ifail == 1) THEN   ! Added by SED
    obj%logL = -HUGE(1._RP)
    RETURN
  ENDIF
  IF(NTIME /= npts) PRINT*,'ERROR: seismogram record length mismatch',NTIME,npts
  !! GF predictions:
  DpredRG(iaz,1:NTIME)  = spike(1:NTIME,1)   ! Jan Dettmer: these are impulse respones
  DpredTG(iaz,1:NTIME)  = spike(1:NTIME,2)
  DpredVG(iaz,1:NTIME)  = spike(1:NTIME,3)
  call pwaveqn(npts,dt,spike(1:NTIME,1),spike(1:NTIME,3),rfrad,shift2,wl,width2)
  IF(I_T == 1) call pwaveqn(npts,dt,spike(1:NTIME,2),spike(1:NTIME,3),rftra,shift2-dt,wl,width2)

  obj%DpredR(iaz,1:NTIME) = rfrad(1:NTIME)*taper_dpred(1:NTIME)
  !obj%DpredR(iaz,1:NTIME) = rfrad(1:NTIME)/MAXVAL(rfrad(1:NTIME))
  !obj%DpredR(iaz,1:NTIME) = rfrad(1:NTIME)/0.57 !J Smale: NORMALIZATION FOR ITERDECON
  obj%DresR(iaz,1:NTIME) = obj%DobsR(iaz,1:NTIME)-obj%DpredR(iaz,1:NTIME)
  IF(I_T == 1)THEN
    obj%DpredT(iaz,1:NTIME) = rftra(1:NTIME)*taper_dpred(1:NTIME)
    obj%DresT(iaz,1:NTIME) = obj%DobsT(iaz,1:NTIME)-obj%DpredT(iaz,1:NTIME)
  ENDIF

ENDDO

IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='DpredRG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=51,FILE='DpredVG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=52,FILE='DpredTG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=53,FILE='DpredRFR.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  OPEN(UNIT=54,FILE='DpredRFT.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  DO iaz = 1,NRF1
    WRITE(50,201) DpredRG(iaz,1:NTIME)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(51,201) DpredVG(iaz,1:NTIME)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(52,201) DpredTG(iaz,1:NTIME)
  ENDDO
  DO iaz = 1,NRF1
    WRITE(53,201) obj%DpredR(iaz,1:NTIME)
  ENDDO
  IF(I_T == 1)THEN
    DO iaz = 1,NRF1
      WRITE(54,201) obj%DpredT(iaz,1:NTIME)
    ENDDO
  ENDIF
  CLOSE(50)
  CLOSE(51)
  CLOSE(52)
  CLOSE(53)
  CLOSE(54)
  201 FORMAT(1024ES18.10)
ENDIF

ibadlogL = 0
IF(IAR == 1)THEN
   obj%DarR  = 0._RP
   !!
   !!  Compute autoregressive model
   !!
   CALL ARPRED_RF(obj,1,1,NTIME)
   !! Recompute predicted data as ith autoregressive model
   obj%DresR = obj%DresR-obj%DarR

   !! Check if predicted AR model data are outside max allowed bounds
   CALL CHECKBOUNDS_ARMXRF(obj,ibadlogL)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpR = 0._RP
  EtmpT = 0._RP
  DO iaz = 1,NRF1
    IF(ICOV == 0)THEN
    !!
    !! impilicitly sample over sigma
    !!
      EtmpR(iaz) = -REAL(NTIME,RP)/2._RP * LOG( SUM(obj%DresR(iaz,:)**2._RP) / REAL(NTIME,RP) )
      IF(I_T == 1)THEN
        EtmpT(iaz) = -REAL(NTIME,RP)/2._RP * LOG( SUM(obj%DresT(iaz,:)**2._RP) / REAL(NTIME,RP) )
      ENDIF
    ELSEIF(ICOV == 1)THEN
      !!
      !! Sample over sigma (one for all freqs)
      !!
      EtmpR(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                   -(SUM(obj%DresR(iaz,:)**2._RP)/(2._RP*obj%sdparR(iaz)**2._RP)&
                   +REAL(NTIME,RP)*LOG(obj%sdparR(iaz)))
      IF(I_T == 1)THEN
        EtmpT(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                    -(SUM(obj%DresT(iaz,:)**2._RP)/(2._RP*obj%sdparT(iaz)**2._RP)&
                    +REAL(NTIME,RP)*LOG(obj%sdparT(iaz)))
      ENDIF
    ELSEIF(ICOV==2)THEN
      !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
      !! dev parameter... (see Dettmer et al. 2014 GJI)
      EtmpR(iaz) = -DOT_PRODUCT(MATMUL(obj%DresR(iaz,:),Cdi),obj%DresR(iaz,:))/(2._RP*obj%sdparR(iaz)) &
                   -REAL(NTIME,RP)/2._RP*LOG(obj%sdparR(iaz))
      IF(I_T == 1)THEN
        EtmpT(iaz) = -DOT_PRODUCT(MATMUL(obj%DresT(iaz,:),Cdi),obj%DresT(iaz,:))/(2._RP*obj%sdparT(iaz)) &
                     -REAL(NTIME,RP)/2._RP*LOG(obj%sdparT(iaz))
      ENDIF
    ENDIF 
  ENDDO
  logL = SUM(EtmpR)
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
END SUBROUTINE LOGLHOOD_RF_ray3d
!!=======================================================================

SUBROUTINE LOGLHOOD_RF_shibutani(obj,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted R and V components
!!
USE RJMCMC_COM
USE RAY3D_COM
USE ieee_arithmetic
IMPLICIT NONE

INTEGER(KIND=IB)                 :: ifr,ipar,id,ilay,iparcur,ipred,iaz
INTEGER(KIND=IB)                 :: ibadlogL
TYPE (objstruc)                  :: obj
REAL(KIND=RP),DIMENSION(NRF1)    :: EtmpR,EtmpV,EtmpT
REAL(KIND=RP)                    :: logL
INTEGER(KIND=IB),DIMENSION(obj%nunique):: idxh
INTEGER(KIND=IB),DIMENSION(NRF1) :: idxrand
REAL(KIND=RP)                    :: tstart, tend, tcmp   ! Overall time 
!! Forward specific variables
INTEGER(KIND=IB):: ifail
LOGICAL:: ISNAN

! ray3d specific variables
INTEGER(KIND=IB):: nly,jlay,iqq,npts
REAL(KIND=SP)   :: hsub,hsubinc,bazs,ps
real spike(maxsamp,3),rad(maxsamp),tra(maxsamp),ver(maxsamp), &
     rfrad(maxsamp),qa(maxlay), qb(maxlay),din,fs
! Local constants
real, parameter :: rad2 = 0.017453292, v60 = 8.043

!! Attenuations for Shibutani code (hard coded to values Thomas & Rhys used)
qa = 1450.
qb = 600.

!! zero traces
iqq = 0
iaz = 1
!DO iaz = 1,NRF1
! calculate seismograms using ray3d for each trace

!  bazs = baz2(iaz)   !! baz is in deg. (for ray3d), baz2 is in radians (for raysum)
  ps   = REAL(slow2(iaz),SP) !! ps is for km/s units here (converted in read_input)!
  npts = NTIME
  dt = REAL(sampling_dt,SP)
  fs = 1._SP/dt

  !call ray3dn(strike,dip,h,alpha,beta,rho,bazs,ps,dt, &
  !     t_af+shift2-dt,shift2-dt,obj%nunique+1,npts,spike,iqq,ifail)

  !! Thomas/Rhys fixed ray param...
  !ppara = sin(angle*rad)/v60
  din = asin(ps*beta(obj%nunique+1)*vpvs(obj%nunique+1))/rad2

  !! Shibutani code to compute RF for stack of horiz. stratified layers:
  call theo(obj%nunique+1, beta,rho, h, vpvs, &
            qa, qb, fs, din, width2, wl, &
            shift2, npts, rfrad)


!  IF (ifail == 1) THEN   ! Added by SED
!    obj%logL = -HUGE(1._RP)
!    RETURN
!  ENDIF
  IF(NTIME /= npts) PRINT*,'ERROR: seismogram record length mismatch',NTIME,npts
  !! GF predictions:
!  DpredRG(iaz,1:NTIME)  = spike(1:NTIME,1)   ! Jan Dettmer: these are impulse respones
!  DpredTG(iaz,1:NTIME)  = spike(1:NTIME,2)
!  DpredVG(iaz,1:NTIME)  = spike(1:NTIME,3)
!  call pwaveqn(npts,dt,spike(1:NTIME,1),spike(1:NTIME,3),rfrad,shift2,wl,width2)
!  IF(I_T == 1) call pwaveqn(npts,dt,spike(1:NTIME,2),spike(1:NTIME,3),rftra,shift2-dt,wl,width2)

  obj%DpredR(iaz,1:NTIME) = REAL(rfrad(1:NTIME),RP)
  obj%DresR(iaz,1:NTIME) = obj%DobsR(iaz,1:NTIME)-obj%DpredR(iaz,1:NTIME)

!ENDDO

IF(IMAP == 1)THEN
  OPEN(UNIT=50,FILE='DpredRG.txt',FORM='formatted',STATUS='replace', &
  ACTION='WRITE',POSITION='REWIND',RECL=8192)
  WRITE(50,201) rfrad(1:NTIME)
  CLOSE(50)
  201 FORMAT(1024ES18.10)
ENDIF

ibadlogL = 0
IF(IAR == 1)THEN
   obj%DarR  = 0._RP
   !!
   !!  Compute autoregressive model
   !!
   !PRINT*,rank,'before'
   CALL ARPRED_RF(obj,1,1,NTIME)
   !! Recompute predicted data as ith autoregressive model
   obj%DresR = obj%DresR-obj%DarR

   !! Check if predicted AR model data are outside max allowed bounds
   CALL CHECKBOUNDS_ARMXRF(obj,ibadlogL)
   !PRINT*,rank,'after'
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpR = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one for all freqs)
    !!
    DO iaz = 1,NRF1
      EtmpR(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
                  -(SUM(obj%DresR(iaz,:)**2._RP)/(2._RP*obj%sdparR(iaz)**2._RP)&
                  +REAL(NTIME,RP)*LOG(obj%sdparR(iaz)))
!      IF(I_T == 1)THEN
!        EtmpT(iaz) = LOG(1._RP/(2._RP*PI2)**(REAL(NTIME,RP)/2._RP)) &
!                    -(SUM(obj%DresT(iaz,:)**2._RP)/(2._RP*obj%sdparT(iaz)**2._RP)&
!                    +REAL(NTIME,RP)*LOG(obj%sdparT(iaz)))
!      ENDIF
    ENDDO
  ENDIF
  logL = SUM(EtmpR)
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
END SUBROUTINE LOGLHOOD_RF_shibutani
!!=======================================================================

SUBROUTINE LOGLHOOD_SWD(obj,curmod,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted SWD data and compute logL_SWD
!!
USE RJMCMC_COM
USE ieee_arithmetic
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
REAL(KIND=SP), ALLOCATABLE, DIMENSION(:,:) :: curmodtest

IF(IMAP == 1)THEN
  PRINT*,'CURMOD IN SWD'
  DO ilay=1,obj%nunique+1
    WRITE(*,206)ilay,curmod(ilay,1:10)
  ENDDO
  206   FORMAT(I3,10F12.4)
ENDIF

!!
!!  curmod = thick  rho  alph  beta  %P  %S  tr  pl  st  dip 
!!

!! Find lowest valid entry for perturbations
!factvs   = obj%par(obj%nunique*NPL+1)  !! Vs is in km/s here
!factvpvs = obj%par(obj%nunique*NPL+2)
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,2) == 1)THEN
    factvs   = obj%voro(ilay,2)  !! Vs is in km/s here
    EXIT
  ENDIF
ENDDO
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,3) == 1)THEN
    factvpvs = obj%voro(ilay,3)
    EXIT
  ENDIF
ENDDO
!WRITE(*,*) 'vel_prem = ', SHAPE(vel_prem), vel_prem
curmod2  = 0.
curmod2(1:obj%nunique+1,:) = curmod(1:obj%nunique+1,:)
!! last layer thickness
!curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%par((obj%nunique-1)*NPL+1)*1000.
curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%voro(obj%k,1)*1000.

!PRINT*,'obj voro'
!DO ilay=1,obj%k+1
!  PRINT*,obj%voro(ilay,:)
!ENDDO
!PRINT*,'obj par',obj%par(1:obj%NFP+2)
!PRINT*,'fact',factvs,factvpvs,obj%par((obj%nunique-1)*NPL+1)

curmod2(obj%nunique+2:obj%nunique+NPREM,1)   = (vel_prem(1,2:NPREM)-vel_prem(1,1:NPREM-1)) * 1000. !! thickness in km
curmod2(obj%nunique+1+NPREM,1)               = 0.                                                  !! HS thickness is 0 
curmod2(obj%nunique+2:obj%nunique+1+NPREM,2) = vel_prem(4,1:NPREM) * 1000.            !! Density
curmod2(obj%nunique+2:obj%nunique+1+NPREM,4) = (vel_prem(2,1:NPREM) + factvs) * 1000. !! Vs
curmod2(obj%nunique+2:obj%nunique+1+NPREM,3) = (vel_prem(2,1:NPREM)+factvs)*(vel_prem(3,1:NPREM)+factvpvs)*1000. !! Vp
IF(IMAP == 1)THEN
  WRITE(*,*) 'curmod2 (including PREM)'
  DO ilay=1,obj%nunique+1+NPREM+1
    WRITE(*,206)ilay,curmod2(ilay,1:10)
  ENDDO
  !206   FORMAT(I3,10F12.4)
ENDIF

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
!ALLOCATE(curmodtest(maxlay,10))
!curmodtest = 0.0
!curmodtest(1:obj%nunique+1,1) = (/ 40000.0, 60000.0, 100000.0, 0.0  /)
!curmodtest(1:obj%nunique+1,2) = (/ 2593.3594, 3098.5693, 3308.5215, 3399.7598  /)
!curmodtest(1:obj%nunique+1,3) = (/ 5600.0, 7560.0, 8160.0, 8400.0  /)
!curmodtest(1:obj%nunique+1,4) = (/ 3200.0, 4200.0, 4800.0, 4200.0  /)
!CALL dispersion(obj%nunique+1,curmodtest(1:obj%nunique+1,2)/1000., & 
!     curmodtest(1:obj%nunique+1,3)/1000.,curmodtest(1:obj%nunique+1,4)/1000.,&
!     curmodtest(1:obj%nunique+1,1)/1000.,DpredSWD,&
!     periods,NDAT_SWD,ierr_swd)
!IF(IMAP == 1)THEN
!  PRINT*,'CURMODTEST IN SWD'
!  DO ilay=1,obj%nunique+1
!    WRITE(*,296)ilay,curmodtest(ilay,1:10)
!  ENDDO
!  296   FORMAT(I3,10F12.4)
!ENDIF
!CALL dispersion(4,(/ 2593.3594, 3098.5693, 3308.5215, 3399.7598  /)/1000., & 
!     (/ 5600.0, 7560.0, 8160.0, 8400.0  /)/1000.,(/ 3200.0, 4200.0, 4800.0, 4200.0  /)/1000.,&
!     (/ 40000.0, 60000.0, 100000.0, 0.0  /)/1000.,DpredSWD,&
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
  IF(ICOV == 0)THEN
  !!
  !! implicitly sample over sigma
  !!
  DO imod = 1,NMODE
    EtmpSWD(imod) = -REAL(NDAT_SWD,RP)/2._RP * LOG( SUM(obj%DresSWD(imod,:)**2._RP) / REAL(NDAT_SWD,RP) )
  ENDDO
  ELSEIF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one per mode)
    !!
    DO imod = 1,NMODE
      EtmpSWD(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_SWD,RP)/2._RP)) &
                    -(SUM(obj%DresSWD(imod,:)**2._RP)/(2._RP*obj%sdparSWD(imod)**2._RP)&
                    +REAL(NDAT_SWD,RP)*LOG(obj%sdparSWD(imod)))
    ENDDO
  ELSEIF(ICOV == 2)THEN
    !! Empirical cov mat estimate and magnitude scaling (scaling is done via std
    !! dev parameter... (see Dettmer et al. 2014 GJI)  scaling: Cd = sC'd not Cd = (s^2)C'd 
    DO imod = 1,NMODE
      !EtmpSWD(imod) = -DOT_PRODUCT(MATMUL(obj%DresSWD(imod,:),CdiSWD),obj%DresSWD(imod,:)) !TRANSPOSE(obj%DresSWD(imod,:)))
      !EtmpSWD(imod) = EtmpSWD(imod)/(2._RP*obj%sdparSWD(imod)) &
       !    -REAL(NDAT_SWD,RP)/2._RP*LOG(obj%sdparSWD(imod))
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
IF (IMAP == 1) WRITE(*,*) 'logL_SWD = ', logL
RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_SWD

!!=======================================================================

SUBROUTINE LOGLHOOD_MT(obj,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted MT data and compute logL_MT
!!
USE RJMCMC_COM
USE ieee_arithmetic

IMPLICIT NONE
INTERFACE
   SUBROUTINE MT1DforwardB3(ZMT, appRes, phase, sigma, h, angFre, Nl)
     USE data_type !module contains integer constants for sigle and double precision (SP and RP)
     IMPLICIT NONE
     REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: sigma    !vector of layers conductivities
     REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: h        !vector of layers thicknesses
     REAL(KIND = RP), DIMENSION(:),   INTENT(IN)  :: angFre   !vector of ANGULAR frequencies
     INTEGER,                         INTENT(IN)  :: Nl
     COMPLEX(KIND = RP), DIMENSION(:),   INTENT(OUT) :: ZMT      !vector of impedances
     REAL(KIND = RP),    DIMENSION(:),   INTENT(OUT) :: appRes   !vector of apparent resistivities
     REAL(KIND = RP),    DIMENSION(:),   INTENT(OUT) :: phase    !vector of phases in radians
   END SUBROUTINE MT1DforwardB3
END INTERFACE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                      :: ipred,ierr_swd,imod,ibadlogL,ilay
TYPE (objstruc)                       :: obj
REAL(KIND=RP),DIMENSION(maxlay)       :: sigma, log10sigma
!REAL(KIND=SP),DIMENSION(maxlay+NPREM,10):: curmod2
COMPLEX(KIND=RP),DIMENSION(NDAT_MT)   :: ZMT
REAL(KIND=RP),DIMENSION(NDAT_MT)      :: appRes, phase, freqMT
REAL(KIND=RP),DIMENSION(2*NDAT_MT)    :: DpredMT
!REAL(KIND=RP),DIMENSION(NMODE)        :: EtmpSWD
REAL(KIND=RP)                         :: logL,factvs,factvpvs
REAL(KIND=RP)                         :: tstart, tend, tcmp   ! Overall time 
LOGICAL :: ISNAN
COMPLEX(KIND=RP), DIMENSION(NDAT_MT)  :: Complex_res
!REAL(KIND = RP),PARAMETER             :: PI = DACOS(-1.0D0)

IF(IMAP == 1) WRITE(*,*) 'maxlay = ', maxlay

IF(IMAP == 1)THEN
  PRINT*, 'VORO:'
  DO ilay=1, obj%nunique+1
    WRITE(*,*) ilay, obj%voro(ilay,:)
  END DO
END IF

IF(IMAP == 1)THEN
  PRINT*, 'PAR:'
  WRITE(*,*) obj%par(1:NPL*obj%k)
END IF

!! define isig, inducator of sigma in par

DO ilay=1,obj%nunique
  !isig = (ilay-1)*NPL + NPL
  log10sigma(ilay) = obj%par((ilay-1)*NPL + NPL)
  sigma(ilay) = 10.0**obj%par((ilay-1)*NPL + NPL)
ENDDO
log10sigma(obj%nunique+1) = obj%par(obj%nunique*NPL+NPL-1)
sigma(obj%nunique+1) = 10.0**obj%par(obj%nunique*NPL+NPL-1)

IF(IMAP == 1)THEN
  PRINT*,'CONDUCTIVITY and h IN MT'
  DO ilay=1,obj%nunique+1
    WRITE(*,206)ilay,sigma(ilay), 1000.0*obj%hiface(ilay)
  ENDDO
  206   FORMAT(I3,2F16.4)
ENDIF

!!
!!  curmod = thick  rho  alph  beta  %P  %S  tr  pl  st  dip 
!!
!WRITE(*,*) 'PI2 = ', PI2
!PRINT*, 'Before'
!WRITE(*,*) 'appRes = ', appRes
!WRITE(*,*) 'phase = ', phase
!WRITE(*,*) 'sigma = ', sigma(1:obj%nunique+1)
!WRITE(*,*) 'h = ', 1000.0*obj%hiface(1:obj%nunique+1)
!WRITE(*,*) 'fre = ', obj%freqMT
!WRITE(*,*) 'angfre = ', 2.0*PI*obj%freqMT
!WRITE(*,*) 'k = ', obj%nunique+1
! call MT forward code
freqMT = REAL(obj%freqMT, RP)
CALL MT1DforwardB3(ZMT, appRes, phase, sigma(1:obj%nunique+1), 1000.0_RP*obj%hiface(1:obj%nunique+1), 2.0_RP*PI2*freqMT, obj%nunique+1)
!PRINT*, 'After'
!PRINT*, 'appRes'
!WRITE(*,*) appRes
!PRINT*, 'phase'
!WRITE(*,*) phase
IF (I_ZMT == 1) THEN
  DpredMT(1:NDAT_MT) = REAL(ZMT)
  DpredMT(NDAT_MT+1:2*NDAT_MT) = AIMAG(ZMT)
ELSE
  DpredMT(1:NDAT_MT) = appRes
  DpredMT(NDAT_MT+1:2*NDAT_MT) = phase
END IF

obj%DpredMT(1:2*NDAT_MT) = DpredMT
obj%DresMT(1:2*NDAT_MT) = obj%DobsMT(1:2*NDAT_MT)-obj%DpredMT(1:2*NDAT_MT)

!!!!!!!!!!!!!!!!!!!!!! Pejman: training dataset for training the surrogate
!IF (IMAP==1) THEN
!   cur_modMT = 0.0_RP
!   cur_modMT(1:obj%nunique, 1) = obj%hiface(1:obj%nunique) !! hiface is in km
!   cur_modMT(1:obj%nunique+1, 2) = log10sigma(1:obj%nunique+1)
 
!   PRINT*,'cur_modMT'
!   DO ilay=1,obj%nunique+1
!     WRITE(*,206)ilay,10.0_RP**cur_modMT(ilay, 2), 1000.0_RP*cur_modMT(ilay,1)
!   ENDDO

!END IF
!!!!!!!!!!!!!!!!!!!!!

ibadlogL = 0
!!Add it later, can be off by IAR=0
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
  
  IF(ICOV == 0)THEN
    !! implicitly sample over sigma
    IF(I_ZMT == 1)THEN
      Complex_res = CMPLX(obj%DresMT(1:NDAT_MT), obj%DresMT(NDAT_MT+1:2*NDAT_MT), RP)
      logL = -REAL(NDAT_MT,RP)* LOG( SUM(CONJG(Complex_res)*Complex_res) / REAL(NDAT_MT,RP) )
    ELSE
      logL = -REAL(NDAT_MT,RP)/2._RP * LOG( SUM(obj%DresMT**2._RP) / REAL(NDAT_MT,RP) )
    ENDIF

  ELSEIF(ICOV == 1)THEN
    IF(IMAGSCALE == 0)THEN
      !!
      !!
      !! 
      IF (I_ZMT == 1) THEN       !!!Cd = (s^2)C'd
        
        Complex_res = CMPLX(obj%DresMT(1:NDAT_MT), obj%DresMT(NDAT_MT+1:2*NDAT_MT), RP)
        logL = -REAL(NDAT_MT,RP)*LOG(PI2*obj%sdparMT(1)**2.0_RP) &
               -SUM(LOG(MTZVAR)) &
               -SUM(CONJG(Complex_res)*Complex_res/MTZVAR)/(obj%sdparMT(1)**2.0_RP)

      ELSE

      logL = LOG(1._RP/(2._RP*PI2)**(REAL(2*NDAT_MT,RP)/2._RP)) &
           -(SUM(obj%DresMT**2._RP)/(2._RP*obj%sdparMT(1)**2._RP)&
           +REAL(2*NDAT_MT,RP)*LOG(obj%sdparMT(1)))

      END IF !I_ZMT

    ELSE


         logL = -(REAL(2*NDAT_MT,RP)/2._RP)*LOG(2._RP*PI2) &
        -SUM(obj%DresMT**2._RP/ &
        (2._RP*(obj%DobsMT*obj%sdparMT(1))**2._RP)) &
        -REAL(2*NDAT_MT,RP)*LOG(obj%sdparMT(1)) - &
        SUM(LOG(ABS(obj%DobsMT)))

    END IF !IMAGSCALE

  ELSEIF(ICOV == 2)THEN

    IF(I_ZMT == 1)THEN !!!!!!! check for type (real to complex) conversions. logL must be real
!      Complex_res = CMPLX(obj%DresMT(1:NDAT_MT), obj%DresMT(NDAT_MT+1:2*NDAT_MT), RP)
!      logL = -DOT_PRODUCT(MATMUL(obj%Complex_res,CdiMT),Complex_res) / (obj%sdparMT(1)) &
!             -REAL(NDAT_MT,RP)/2._RP * LOG(obj%sdparMT(1))

     Complex_res = CMPLX(obj%DresMT(1:NDAT_MT), obj%DresMT(NDAT_MT+1:2*NDAT_MT), RP)
     logL = -REAL(NDAT_MT,RP)*LOG(PI2*obj%sdparMT(1)**2.0_RP) &
            -SUM(LOG(MTZVAR)) &
            -SUM(CONJG(Complex_res)*Complex_res/MTZVAR)/(obj%sdparMT(1)**2.0_RP)

    ELSE
!    logL = -DOT_PRODUCT(MATMUL(obj%DresMT,CdiMT),obj%DresMT) / (2._RP*obj%sdparMT(1)) &
!           -REAL(NDAT_MT,RP) * LOG(obj%sdparMT(1))

     logL = LOG(1._RP/(2._RP*PI2)**(REAL(2*NDAT_MT,RP)/2._RP)) &
            -(SUM(obj%DresMT**2._RP)/(2._RP*obj%sdparMT(1)**2._RP)&
            +REAL(2*NDAT_MT,RP)*LOG(obj%sdparMT(1)))

    ENDIF  !!I_ZMT

  ENDIF !!ICOV

  !WRITE(*,*) 'MTDATAstd = ', obj%sdparMT(1)
  !WRITE(*,*) 'ieee_is_nan(logL) = ', ieee_is_nan(logL)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
   logL = -HUGE(1._RP)
   ibadlogL = 0
ENDIF
IF (IMAP == 1) WRITE(*,*) 'logL_MT = ', logL
RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_MT

!!=======================================================================



!!=======================================================================

SUBROUTINE MAKE_CURMOD(obj,curmod)
!!=======================================================================
!!
!! Build curmod from model
!!
USE RJMCMC_COM
USE ieee_arithmetic
IMPLICIT NONE
INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ipar,ilay,id,iparcur, NPL2, ijump
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,10):: curmod

NPL2 = NPL
IF(I_MT == 1) NPL2 = NPL -1

ijump = 0
IF(I_MT == 1) ijump = ijump + 1 !! to jump over sigma in par

curmod = 0._RP
!! curmod = (/ thick, rho, alph, beta, %P, %S, tr, pl, st, di /)
id=0
DO ilay=1,obj%nunique+1
  !curmod(ilay,:) = curmod_glob
  IF(I_MT == 1 .AND. ilay > 1) id=id+ijump !! jump over sigma in par 
  DO ipar=1,NPL2
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
  !! for sampling Vs:
  curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
  curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
  curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,3)*curmod(1:obj%nunique+1,4)
ELSEIF(I_VPVS == -1)THEN
  !! for sampling Vs:
  curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
  curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
  curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,3)*1.75_SP
ELSE
  curmod(1:obj%nunique+1,1:4) = curmod(1:obj%nunique+1,1:4)*1000._SP
ENDIF

!! Birch Law for analytical density relationship
!curmod(1:obj%nunique+1,2)   = 1000.*0.77+0.32*curmod(:,3)
!! What Thomas uses:
curmod(1:obj%nunique+1,2) = (2.35+0.036*((curmod(1:obj%nunique+1,3)/1000.)-3.0)**2.)*1000.

!!
!! Convert angles to radians (trend, plunge, strike, dip)
curmod(1:obj%nunique+1,7)  = curmod(1:obj%nunique+1,7)/180._RP * PI
curmod(1:obj%nunique+1,8)  = curmod(1:obj%nunique+1,8)/180._RP * PI
curmod(1:obj%nunique+1,9)  = curmod(1:obj%nunique+1,9)/180._RP * PI
curmod(1:obj%nunique+1,10) = curmod(1:obj%nunique+1,10)/180._RP * PI

!IF(IMAP == 1)THEN
!  DO ilay=1,obj%nunique+1
!    WRITE(*,206)ilay,curmod(ilay,1:10)
!  ENDDO
!  206   FORMAT(I3,10F12.4)
!ENDIF
RETURN
END SUBROUTINE MAKE_CURMOD
!!=======================================================================

SUBROUTINE MAKE_CURMOD_ray3d(obj,curmod)
!!=======================================================================
!!
!! Build curmod from model
!!
USE RJMCMC_COM
USE RAY3D_COM
USE ieee_arithmetic
IMPLICIT NONE
!INCLUDE 'raysum/params.h'

INTEGER(KIND=IB)                 :: ipar,ilay,id,iparcur,ifail,nly, NPL2, ijump
TYPE (objstruc)                  :: obj
REAL(KIND=SP),DIMENSION(maxlay,7):: curmod

!real    :: mod_min(maxlay,7)
!logical :: mod_flag(maxlay,7)
!real    :: h(maxlay),rho(maxlayd),alpha(maxlayd),beta(maxlayd)
!real    :: strike(maxlayd),dip(maxlayd)

!! Scratch variables
real spike(maxsamp,3),rad(maxsamp),tra(maxsamp),ver(maxsamp)
real rfrad(maxsamp),rftra(maxsamp)
real dt_wind(2,maxsamp,maxtr),rf_wind(2,maxsamp,maxtr)

NPL2 = NPL
IF(I_MT == 1) NPL2 = NPL -1

ijump = 0
IF(I_MT == 1) ijump = ijump + 1

ifail  = 0  ! Added to ray3d by SED to signal failure of forward model
curmod = 0._RP

alpha  = 0.
beta   = 0.
vpvs   = 0.
rho    = 0.
h      = 0.
strike = 0.
dip    = 0.

!! Build velocity model
!! First, restore from model parameters to RF parameters
!do ilay=1,nlay
!  do ipar=1,7
!    if (mod_flag(ilay,ipar)) then
!      id=id+1
!      curmod(ilay,ipar)=model(id)
!    else
!      curmod(ilay,ipar)=mod_min(ilay,ipar)
!    endif
!  enddo
!enddo
!! curmod = (/ thick, rho, alph, beta, %P, %S, tr, pl, st, di /)
id = 0
DO ilay=1,obj%nunique+1
  !curmod(ilay,:) = curmod_glob
  IF(I_MT == 1 .AND. ilay > 1) id=id+ijump
  DO ipar=1,NPL2
    id=id+1
    iparcur = idxpar(ipar)
    IF(ilay <= obj%nunique)THEN
      IF(iparcur == 1)THEN
        curmod(ilay,iparcur) = REAL(obj%hiface(ilay),SP)  !! Model works with layer thickness!
      ELSE
        curmod(ilay,iparcur) = REAL(obj%par(id),SP)
      ENDIF
    ELSE
      IF(iparcur == 1)THEN
        curmod(ilay,iparcur) = 0._SP
        id=id-1
      ELSE
        curmod(ilay,iparcur) = REAL(obj%par(id),SP)
      ENDIF
    ENDIF
  ENDDO
ENDDO

!! Birch Law for analytical density relationship
!curmod(1:obj%nunique+1,2)   = 1000.*0.77+0.32*curmod(:,3)
!! What Thomas uses:
curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,4)*curmod(1:obj%nunique+1,3)
curmod(1:obj%nunique+1,2) = (2.35+0.036*((curmod(1:obj%nunique+1,3))-3.0)**2.)

!!
!! Convert curmod to ray3d parameters
!!
h = 0.; alpha = 0.; beta = 0.; rho = 0.; strike = 0.; dip = 0.
nly=0
beta(1:obj%nunique+1)  = curmod(1:obj%nunique+1,4)
alpha(1:obj%nunique+1) = curmod(1:obj%nunique+1,3)
rho(1:obj%nunique+1)   = curmod(1:obj%nunique+1,2)
h(1:obj%nunique)       = curmod(1:obj%nunique,1)
strike(1:obj%nunique)  = curmod(2:obj%nunique+1,6)
dip(1:obj%nunique)     = curmod(2:obj%nunique+1,7)
vpvs(1:obj%nunique+1)  = alpha(1:obj%nunique+1)/beta(1:obj%nunique+1)

!!
!! Convert to m/s units for SWD...
!!
curmod(1:obj%nunique+1,1) = curmod(1:obj%nunique+1,1)*1000._SP
curmod(1:obj%nunique+1,2) = curmod(1:obj%nunique+1,2)*1000._SP
curmod(1:obj%nunique+1,3) = curmod(1:obj%nunique+1,3)*1000._SP
curmod(1:obj%nunique+1,4) = curmod(1:obj%nunique+1,4)*1000._SP
!! print velocity model on screen
IF(IMAP == 1)THEN
  write(6,*)'nly=',obj%nunique+1
  write(6,*)'Velocity Model:'
  write(6,*)'layer  Vp        Vs    rho  thickness strike   dip'
  do ilay=1,obj%nunique+1
    write(*,'(i3,1x,7f8.2)') ilay,alpha(ilay),beta(ilay),rho(ilay),& 
                             h(ilay),strike(ilay),dip(ilay)
  enddo
ENDIF

RETURN
END SUBROUTINE MAKE_CURMOD_ray3d
!!=======================================================================

SUBROUTINE INTERPLAYER_novar(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! This does not allow for variable layer complexity.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot,ilim
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP)                    :: vref,vpvsref
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort

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

tmpsort = 0._DRP
tmpsort = REAL(obj%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%voro(1:obj%k,:) = REAL(tmpsort,RP)

obj%ziface = 0._RP
obj%ziface(1:obj%k-1) = obj%voro(2:obj%k,1)

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%k-1
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%k,1) = obj%ziface(1:obj%k)
partmp(1:obj%k,2:NPL) = obj%voro(1:obj%k,2:NPL)

IF(I_VREF == 1 .AND. (I_RV == -1 .OR. I_RV == 1 .OR. I_SWD == 1) )THEN
  DO ivo = 1,obj%k
    CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
    partmp(ivo,2) = vref + partmp(ivo,2)      !! if only I_MT == 1, this produces wrong result
    partmp(ivo,3) = vpvsref + partmp(ivo,3)   !! if only I_MT== 1, partmp(ivo, 3) produces runtime error: no third column
    !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
  ENDDO 
END IF
!!
!! Apply reference profile:
!!
!IF(I_VREF == 0)THEN
 ! partmp(1:obj%k,2:NPL) = obj%voro(1:obj%k,2:NPL)
!ELSE
 ! DO ivo = 1,obj%k
  !  CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
   ! partmp(ivo,2) = vref + obj%voro(ivo,2)
   ! partmp(ivo,3) = vpvsref + obj%voro(ivo,3)
   ! !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
 ! ENDDO 
!ENDIF

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
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort
REAL(KIND=DRP),DIMENSION(NLMX*NPL):: tmpsort2

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

IF((I_RV==-1 .OR. I_RV==1 .OR. I_SWD==1) .AND. IDIP == 1)THEN
  NPL_tmp = NPL-1
ELSE
  NPL_tmp = NPL
ENDIF

tmpsort = 0._DRP
tmpsort = REAL(obj%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%voro(1:obj%k,:) = REAL(tmpsort,RP)

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
IF((I_RV==-1 .OR. I_RV==1 .OR. I_SWD==1) .AND. I_VREF == 1)THEN
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
tmpsort2 = 0._DRP
tmpsort2 = REAL(obj%ziface(1:ntot),DRP)
CALL QSORTC1D(tmpsort2(1:ntot))
obj%ziface(1:ntot) = REAL(tmpsort2,RP)

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
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxar(ifr) == 1)THEN
   k = 1
   obj%DarR(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   dres1 = obj%DresR(ifr,idata:idatb)

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
   obj%DarR(ifr,idata:idatb) = dar1
   obj%DarR(ifr,idata) = 0._RP
   obj%DarR(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_RF
!!=======================================================================

SUBROUTINE ARPRED_SWD(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
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
   IF(MAXVAL(obj%DarR(ifr,:)) > armxH)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
   IF(MINVAL(obj%DarR(ifr,:)) < -armxH)THEN
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
