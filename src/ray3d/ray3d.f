!===================================================================================== 
      SUBROUTINE user_init(ntr)
!=====================================================================================
! Initialize ray3d forward modelling
        USE ray3d_com
        IMPLICIT NONE
        !include 'ray3d_param.inc'
        !include 'na_param.inc'
        
        integer ntr,ios,nlay2,i,k,j,nd,ilay,iparam
        real c
        !real ranges(2,nd_max)
        !real mod_min(maxlay,7),mod_max(maxlay,7)
        !logical mod_flag(maxlay,7)
        !real baz(maxtr),slow(maxtr)
        !real data_rf(2,maxsamp,maxtr),weight(maxtr)

! Common block for rays3d calling routines
!        common /ray3d_com/nlay,mod_min,mod_flag,baz,slow,
!     *     ntr,data_rf,dt,shift,t_af,nsamp,weight,width,
!     *     w1,w2,c,time_window,misfit_type

! Initiation
      mod_max = 0.
      mod_min = 0. 
     
! Read in limits of model parameters from model.max and model.min
!        open(iounit1,file='model.max',status='old',iostat=ios)
!        if (ios.ne.0) then
!          write(*,*) 'Open file error: model.max. Program stops.'
!          stop
!        endif
!        open(iounit2,file='model.min',status='old',iostat=ios)
!        if (ios.ne.0) then
!          write(*,*) 'Open file error: model.min. Program stops.'
!        endif
!        read(iounit1,*) nlay
!        read(iounit2,*) nlay2
!        if (nlay.ne.nlay2) then
!          write(*,*) 'ERROR: number of model parameters does not match'
!          write(*,*) 'Check files model.max and model.min'
!          write(*,*) 'Program stops.'
!        else
!          do i=1,nlay
!            read(iounit1,*)k,(mod_max(i,j),j=1,7)
!            read(iounit2,*)k,(mod_min(i,j),j=1,7)
!          enddo
!        endif

! Shuffle parameters into format NA wants -- use modflag to
! keep track of which parameters are used.
!        nd=0
!        do ilay=1,nlay
!          do iparam=1,7
!            mod_flag(ilay,iparam) = (mod_min(ilay,iparam) .ne.
!     *                               mod_max(ilay,iparam))
!            if (mod_flag(ilay,iparam)) then
!              nd=nd+1
!              if (nd .gt. nd_max) then
!                write (*,*) 'ERROR -- too many variable parameters!'
!                stop
!              else
!                ranges(1,nd)=mod_min(ilay,iparam)
!                ranges(2,nd)=mod_max(ilay,iparam)
!              end if              
!            end if
!          end do
!        end do   
        
        write (*,*) 'Reading data...'
        open(unit=iounit1,file='rfdata.in',status='old',iostat=ios)
!        if (ios.eq.0) then
          read(iounit1,*) ntr,nsamp,misfit_type,c,width,shift,t_af,
     +                    dt,time_window

          !ntr = 1  ! *******************************

!         write(*,*)  'number of traces:',ntr,' number of samples=',nsamp
!         !write(*,20) ' Water leveling constant:',c
!         !write(*,20) ' Gaussian width         :',width
!         !write(*,20) ' Time before 1st arrival:',shift
!         write(*,20) ' Time after 1st arrival :',t_af
!         !write(*,20) ' Sampling interval      :',dt
!         write(*,20) ' Inversion time window  :',time_window
!         !write(*,20) ' R-/T-comp weighting    :',w1,w2
!         write(*,21) ' Misfit calculation type:',misfit_type
!20       format(a25,2f8.3)
!21       format(a25,i3)
!         do i=1,ntr
!            read(iounit1,*) baz(i),slow(i),weight(i)
!            write(*,10) i,baz(i),slow(i),weight(i)
!10          format(' trace ',I4,': baz=',f10.4,' slowness=',f9.5,
!     +        ' Weight: ',f4.2)
!            do j=1,nsamp
!              read(iounit1,*) data_rf(1,j,i),data_rf(2,j,i)
!            enddo
!          enddo
!        else
!          write(*,*) 'Open file error: rfdata.in   PROGRAM STOPS'
!          close(iounit1)
!          stop
!        endif

      return
      end
     
! ----------------------------------------------------------
! Sort an array of n real-valued elements, from smallest (most negative)
! to largest (most positive). Actual elements aren't moved -- index
! returns new positions. Uses a basic algorithm (insertion sort)
! since we're sorting few elements. m is the number of elements in a;
! n is the number of elements in index (n <= m).
!
! Index should be pre-initialized with initial locations.
        
      subroutine sort_b(a,indx,m,n)
      
        implicit none
        integer m,n
        integer indx(n),j,pos,i
        real a(m),value
        
        if (n .le. 1) then
          return
        end if
        
        do j=2,n
          i=indx(j)
          value=a(i)          
          pos=j
          do while ((pos .gt. 1) .and. (value .lt. a(indx(pos-1))))
            indx(pos)=indx(pos-1)
            pos=pos-1
          end do
          indx(pos)=i
        end do
      end
!========================================================================

      SUBROUTINE ray3dn(str,di,zh,alp,bet,rh,bazs,ps,
     *   dts,duras,delays,nlr,npts,spik,iqq,ifail)

!========================================================================
!
!------------------------------------------------------------------------
!    calculates travel times, azimuthal anomalies, ray parameter
!      anomalies for primary and multiple converted waves
!      in a dipping structure, implementation of the method of
!      Langston (1977; bssa)
!
!    Written by T.J. Owens, march 1982, revised innumerable times
!                                       since then
!    Version A-1; Revised May 1987
!                 Revised Oct 1989 for round-off problems around
!                   Loop 32 - mod by tjo
!------------------------------------------------------------------------

      USE ray3d_par
      !include 'ray3d_param.inc'
      dimension str(maxlayd),di(maxlayd),zh(maxlayd),alp(maxlayd),
     *  bet(maxlayd),rh(maxlayd),spik(maxsamp,3)
      dimension strike(maxlayd),dip(maxlayd),z(maxlayd),
     *alpha(maxlayd),beta(maxlayd),rho(maxlayd),eta(3,maxlayd),
     *q(3,500000),q0(3),v(2,maxlayd),qv(500000),
     *qloc1(3),qloc2(3),a(3,3),iface(500000),layer(maxlayd),
     *mulyr(maxlayd),amag(3,500000),hmag(3,500000),raymag(3,500000),
     *rayhil(3,500000),exmuls(maxlayd),raytim(500000),spike(500000,3),
     *direct(3)
      logical yes,pors,again,amps,ppps(maxlayd),free,
     *mormul(maxlayd),amps1,pps2,mormu2,tmp
      integer trans,refl,type,exmuls,ippps(maxlayd),ounit

      common /cord/ a
      common /amcal/ qloc1,qloc2,vb,va,sinib,sinia,vp1,vs1,rho1,
     *               vp2,vs2,rho2,free,type
      common /transm/ q,qv,v,alpha,beta,rho,strike,dip,iface,jhilb,
     *                amag,hmag,layer,amps,trans,refl,nlyrs
      common /raywrt/ eta,z,raymag,rayhil,raytim,ntim,p0r,pors,
     *                oldlyr,q0,direct,tdirec,baz1
      common /innout/ inunit,ounit
!
!  ray3d generates specific output file names for synthetics
!    spike series are named "spike name"_sp.[zrt]
!    if synthetic is convolved with a source function,
!    the synthetic is "synthetic name"_sy.[zrt]
!    where "spike name" and "synthetic name" are requested by
!    the program.
!
      rad(deg)=deg/57.2957795

      ifail = 0   ! Added by SED

      inunit=5
      ounit=6
      nlyrs=nlr

      do 1819 i1=1,nlyrs
         strike(i1)=str(i1)
         dip(i1)=di(i1)
         alpha(i1)=alp(i1)
         beta(i1)=bet(i1)
         z(i1)=zh(i1)
         rho(i1)=rh(i1)
1819  continue
! 
!  signal is a RIDGE routine used to request that
!  floating point underflows be ignored, this will
!  probably be deleted on most machines
!
!     call iniocm
!     write(ounit,120)
!
!   all output is to file ray3d.out
!
!     open(unit=9,file='ray3d.out',form='formatted')
!     rewind 9
!     write(9,120)
      again=.false.
! 120 format(' ray tracer for dipping structures',/)
!
!     adjust input values from rdlyrs to necessary form
!
      tmpz1=z(1)
      z(1)=0.
      tmps1=strike(1)
      strike(1)=0.
      tmpd1=dip(1)
      dip(1)=0.
      do 48 i48=2,nlyrs
      tmps2=strike(i48)
      tmpd2=dip(i48)
      strike(i48)=tmps1
      dip(i48)=tmpd1
      tmps1=tmps2
      tmpd1=tmpd2
      tmpz2=z(i48)
      z(i48)=z(i48-1)+tmpz1
      tmpz1=tmpz2
   48 continue

!
!     ask all initial questions
!
    6  p0=ps
       p0r=p0
!      baz=ask('Back azimuth of incident ray: ')
       baz=bazs
       baz1=baz
       pors=.true.
!      pors=yesno('P-wave (y or n) ? ')
       if(pors) go to 16
         amps=.false.
         go to 15
!  16 amps=yesno('Calculate any amplitudes (y or n) ? ')
   16 amps=.true.
      if(.not.amps) go to 15
!         pps2=yesno('Pp and Ps only (y or n) ? ')
!         pamp=ask('Incident p amplitude = ')
      pps2=.false.
!     pps2=.true.
      pamp=10.0
   15 sini=p0*alpha(nlyrs)
      if(.not.pors) sini=p0*beta(nlyrs)
      numint = nlyrs -1
      do 22 i22=1,numint
         layer(i22)=i22+1
         mulyr(i22)=0
         ppps(layer(i22))=.false.
         if(pps2) ppps(layer(i22))=.true.
         mormul(layer(i22))=.false.
   22 continue
      ifl=0
   64 ifl=ifl+1
!     write(ounit,107)
! 107 format(' your layer ray tracing parameters are: ',//,
!    *       'interface  ppps  mormul ')
!     do 21 i21=1,numint
!       write(ounit,105) layer(i21),ppps(layer(i21)),mormul(layer(i21))
!  21 continue
! 105 format(5x,i3,5x,l1,5x,l1)
!      if(yesno('OK ? (y or n) ')) go to 18
      if(ifl.eq.1)tmp=.false.
      if(ifl.eq.2)tmp=.true.
      if(tmp) go to 18
!     write(ounit,102)
! 102 format(' enter the # of interfaces to trace from (i2)')
!     read(inunit,103) numin2
      numin2=nlyrs-1
      if(numin2.le.0) go to 60
      numint=numin2
!     write(ounit,101)
! 101 format(' enter the interface numbers (40i2)')
!     read(inunit,103) (layer(i),i=1,numint)
! 103 format(40i2)
      do 132 ji=1,numint
  132 layer(ji)=ji+1
!  60 if(.not.yesno('Change PpPs options (y or n) ? ')) go to 70
      tmp=.false.
   60 if(.not.tmp) go to 70
!        write(ounit,108)
! 108 format('enter interface #s which need ppps changed from current',
!    *       ' value (40i2) ')
!        read(inunit,103) (ippps(i),i=1,40)
         do 71 i71=1,numint
            if(ippps(i71).eq.0) go to 70
            if(ppps(ippps(i71))) then
              ppps(ippps(i71))=.false.
            else
              ppps(ippps(i71))=.true.
            endif
   71 continue

c   70 mormu2=yesno('Calculate extra multiples ? (y or n) ')
   70 mormu2=.true.
      if(.not.mormu2) go to 69
c     write(ounit,104)
c 104 format(' enter interface numbers for extra multiple',
c    *       '  calculations (40i2)')
c     read(inunit,103) (mulyr(i),i=1,30)
c
      do 20 i20=1,100
         mulyr(i20)=layer(i20)
         if(mulyr(i20).ne.0) go to 20
         nmults=i20-1
         go to 61
   20 continue
   61 if(nmults.eq.0) go to 60
c     if(yesno('Calculate extra mults for all rays ? ')) go to 62
c        write(ounit,106)
c 106 format(' enter only interfaces which have rays that need',
c    *       ' extra mults tacked on')
c     read(inunit,103) (exmuls(i),i=1,40)
c
      tmp=.false.
      if(tmp)go to 62
      do 72 i72=1,nlyrs
      exmuls(i72)=layer(i72)
   72 mormul(i72)=.false.
      do 63 i63=1,40
         if(exmuls(i63).eq.0) go to 64
         mormul(exmuls(i63))=.true.
   63 continue
   69 go to 64
   62 do 65 i65=1,numint
   65    mormul(layer(i65))=.true.
      go to 64
   18 nrays=1
      do 181 i181=1,numint
         nr2=9
         if(ppps(layer(i181))) nr2 = 1
         if(mormul(layer(i181))) nr2 = nr2 + nr2*4*nmults
         nrays=nrays + nr2
  181 continue
c       write(ounit,*) "Number of rays =",nrays
c
      if(nrays.le.500000) go to 182
        write(ounit,183) nrays
  183   format(' nrays = ',i5,' is too big - try again ')
        go to 64
  182 if(again) go to 14
c
c   calculate layer interface unit normal vectors in global coordinates
c

      do 1 i1=1,nlyrs
         strike(i1)=rad(strike(i1))
         dip(i1)=rad(dip(i1))
         call norvec(strike(i1),dip(i1),eta(1,i1))
    1 continue

c
c   define incident ray unit vector in global coordinates
c
   14 q0(1)=-sini*cos(rad(baz))
      q0(2)=-sini*sin(rad(baz))
      q0(3)=-sqrt(1. - sini*sini)

c
c   set up velocity arrays and other initital conditions
c
      do 2 i2=1,nlyrs
         if(.not.pors) go to 3
            v(1,i2)=alpha(i2)
            v(2,i2)=beta(i2)
            go to 2
    3    v(1,i2)=beta(i2)
         v(2,i2)=alpha(i2)
    2 continue

      trans=1
      refl=-1
      qv(1)=v(1,nlyrs)
      iface(1)=0
      do 17 i17=1,3
!      call zeerow(amag(i17,1),1,500000)
!c     call zeerow(amag(i17,1),1,500000)
!      call zeerow(hmag(i17,1),1,500000)
      amag(i17,:) = 0.
      hmag(i17,:) = 0.
      q(i17,1)=q0(i17)
      if(.not.amps) go to 17
      amag(i17,1)=pamp*q0(i17)
      hmag(i17,1)=0.
      if(i17.lt.3) go to 17
         vp1=alpha(nlyrs)
         vs1=beta(nlyrs)
         rho1=rho(nlyrs)
         free=.false.
         ntim=1
   17 continue

c
c   s t a r t   r a y   t r a c i n g   s e c t i o n
c
c   find ray unit vectors for the direct ray
c
      ihilb=0
      jhilb=0
      iq=1
      iqq=iq
      call trnsmt(1,nlyrs,iq,1,.true.)
      iqq=iq
c
c   if iq= -999 then a head wave has been generated and the run will bomb
c
      if(iqq.eq.-999) then
csed         write(ounit,133)
csed  133    format(' Immediate problems with head waves ',
csed     *          'on direct wave pass - Check velocity model !!')
c
c ...  Modified to avoid the stop statement here Oct, 29, Serdar
c
c     stop
      ifail = 1  !added by SED
      return 
      else
      endif
      nlr=nlyrs
      call rayfin(nlr,1,0,0,0,.true.,.false.)
      amps1=amps
c
c   calculate the other rays, first all the unconverted rays & their
c     multiples, then the converted waves & their multiples
c     loops 50,52, & 53 do extra multiples, if necessary
c
      do 4 i4=1,2
         do 8 i8=1,numint
            amps=amps1
c
c           if doing the converted waves
c               recalculate the necessary q-vectors
c
            if(i4.eq.1) go to 13
               iq=nlyrs - layer(i8) + 1
               if(.not.amps) go to 28
               vp1=alpha(nlyrs-iq+1)
               vs1=beta(nlyrs-iq+1)
               rho1=rho(nlyrs-iq+1)
   28       loopst=iq
            call trnsmt(loopst,nlyrs,iq,i4,.true.)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
            if(iq.eq.-999) go to 8
c
c           print results for direct converted waves
c
            call rayfin(nlr,i4,0,0,layer(i8),.false.,.false.)
                  if(.not.mormul(layer(i8))) go to 13
                  iqmul=iq
                  do 66 i66=1,nmults
                     if(mulyr(i66).eq.layer(i8).and.
     *                 (.not.ppps(layer(i8)))) go to 66
                     iqi=iqmul
                     do 67 i67=1,2
                        vs1=beta(1)
                        vp1=alpha(1)
                        rho1=rho(1)
                        rho2=0.0
                        vp2=0.
                        vs2=0.
                        call raydwn(iqi,i67,mulyr(i66),iq)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
                        if(iq.eq.-999) go to 66
                        miqdwn=iq
                        do 68 i68=1,2
                           call rayup(miqdwn,i68,mulyr(i66),iq)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
                           if(iq.eq.-999) go to 68
                           call rayfin(iq,0,i67,i68,mulyr(i66),
     *                                 .false.,.true.)
   68                   continue
   67                continue
   66             continue
   13       if(ppps(layer(i8))) amps=.false.
            do 10 i10=1,2
               vs1=beta(1)
               vp1=alpha(1)
               rho1=rho(1)
               rho2=0.0
               vp2=0.
               vs2=0.
               call raydwn(nlyrs,i10,layer(i8),iq)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
               if(iq.eq.-999) go to 10
               iqdown=iq
               do 11 i11=1,2
                  call rayup(iqdown,i11,layer(i8),iq)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
                  if(iq.eq.-999) go to 11
                  call rayfin(iq,i4,i10,i11,layer(i8),.false.,.true.)
                  if(.not.mormul(layer(i8))) go to 11
                  iqmul=iq
                  do 50 i50=1,nmults
                     iqi=iqmul
                     do 52 i52=1,2
                        vs1=beta(1)
                        vp1=alpha(1)
                        rho1=rho(1)
                        rho2=0.0
                        vp2=0.
                        vs2=0.
                        call raydwn(iqi,i52,mulyr(i50),iq)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
                        if(iq.eq.-999) go to 52
                        miqdwn=iq
                        do 53 i53=1,2
                           call rayup(miqdwn,i53,mulyr(i50),iq)
c
c   if iq = -999, then problem phases exist -- this and all subsequent rays
c                      are skipped
c
                           if(iq.eq.-999) go to 53
                           call rayfin(iq,0,i52,i53,mulyr(i50),
     *                                 .false.,.true.)
   53                   continue
   52                continue
   50             continue
   11          continue
   10       continue
    8    continue
    4 continue
      amps=amps1
      if(.not.amps) go to 29
      ntim=ntim-1

c
c  ray3d.amps can be a big file if many rays are traced
c     use with caution
c
c      yes=yesno('Create ray3d.amps ? ')
      yes=.false.
      if(.not.yes) go to 180
c     open(unit=8,file='ray3d.amps',form='formatted')
c     write(8,788) struc,title,nlyrs,p0r,baz
c 788 format(' file: ',a10,' model ',a10,1x,i2,' layers ',
c    *       ' ray parameter ',f7.5,' back az. ',f6.2)
c
c
  180 do 27 i27=1,ntim
         if(i27.gt.500000)then
         write(6,*)'ARRAY DIMENSION IS EXCEEDED27=',i27
         iqq=-999
         return
         else
         endif
         call rtoi(raymag(1,i27),cos(rad(baz)),sin(rad(baz)),-1.,
     *             .false.)
c        call rtoi(rayhil(1,i27),cos(rad(baz)),sin(rad(baz)),-1.,
c    *             .false.)
c        rayhil(3,i27)=-rayhil(3,i27)
         raymag(3,i27)=-raymag(3,i27)
c        if(yes) write(8,122) i27,(raymag(j,i27),j=1,3),
c    *                 (rayhil(j,i27),j=1,3),raytim(i27)
   27 continue
c
c     write(6,*)' N T I M =',ntim
c 122 format(1x,i3,1x,7e15.7)
c     if(jhilb.eq.1) write(ounit,781)
c 781 format(' phase shifted arrivals exist ')
c     yes=yesno('Save this spike ? ')
      yes=.true.
      if(.not.yes) go to 29
       dt=dts
c      dt=ask('Sampling rate (sec): ')
c      dura=ask('Signal duration (secs): ')
c      delay=ask('First arrival delay: ')
      dura=duras
      delay=delays
      npts=ifix(dura/dt + .5) + 1.
      begin = 0.
      spik = 0.  !! JD change
      spike = 0. !! JD change
!      do 30 i30=1,3
!         call zeerow(spik(1,i30),1,1200)
!         call zeerow(spike(1,i30),1,500000)
!c        if(jhilb.eq.1) call zeerow(hilbt(1,i30),1,500000)
!   30 continue


c
c Section below modified on 10/12/89 on Sun-4 to avoid roundoff
c problems identified by John Cassidy at UBC and known to occur
c on Sun-3 versions of ray3d
c Modifications taken from ray3d subroutine in timinv.f
c
      dtby2 = dt/2.
      do 32 j32=1,ntim
         itinc=0
         rtpdel=raytim(j32) + delay
         isampi=rtpdel/dt
         raytm=dt*float(isampi) +dtby2
         if(raytm.lt.rtpdel) itinc=1
         irayl=isampi + itinc + 1
c      write(6,*) 'Here18.5',irayl,isampi,rtpdel,raytim(j32),delay !*******************
         if (raytim(j32)+10. .eq. raytim(j32)) then   ! SED fix for NaN
c            write(6,*) 'NaN Trapped',raytim(j32)
            ifail = 1
            return
         endif
         if(irayl.gt.npts) goto 32
c        write(6,*)'ARRAY DIMENSION IS EXCEEDED=',irayl
c        go to 3332
c        iqq=-999
c         return
c         else
c         endif
c        write(6,*)'baz,ntim,j32,irayl=',baz,ntim,j32,
c    *irayl,raymag(1,j32)
         do 33 i33=1,3
            spike(irayl,i33)=spike(irayl,i33) + raymag(i33,j32)
c           if(jhilb.eq.0) go to 33
c              hilbt(irayl,i33)=hilbt(irayl,i33) + rayhil(i33,j32)
   33       continue
c
c  END of 10/12/89 modifications
c
   32 continue
c
3332     continue
!         do 8191 ii =1,npts
!         spik(ii,1)=spike(ii,1)
!         spik(ii,2)=spike(ii,2)
!         spik(ii,3)=spike(ii,3)
!8191     continue
          spik(1:npts,:) = spike(1:npts,:)
c 
   29 again=.false.
      if(again) go to 6
      if(amps) close(unit=8)
c     close(unit=9)

      return
      end
c
      subroutine anom(q,v,az,p,sini)
c
c   calculates the azimuth and ray parameter of a ray defined by q
c     in a medium of velocity v, assuming the surface is horizontal
c
c fix by H.Kao for Linux f77.
c nfix=1
c q(1) -> q(nfix)
c q(2) -> q(nfix+1)
c q(3) -> q(nfix+2)

csed      dimension q(1)
      dimension q(3)

      deg(rad)=rad*57.2957795
      nfix=1
      cosi=-q(nfix+2)
      sini=sqrt(1.-cosi*cosi)
      p=sini/v
c
c as always vertical incidence case is special
c
      if(sini.gt..0001) then
         sinb=-q(nfix+1)/sini
         cosb=-q(nfix)/sini
         az=atan2(sinb,cosb)
         az=deg(az)
      else
         az=0.0
      endif
  101 format(1x,5e15.7)
      return
      end
c
      function timcor(x1,x2,q0,v)
c
c  finds the time diference between a ray which enters the
c   layering at point x2 to one which enters the layering at
c   x1 if the half space unit ray vector is q0 and the half
c   space velocity is v
c
csed      dimension x1(1),x2(1),q0(1),r(3)
      dimension x1(3),x2(3),q0(3),r(3)


      do 1 i=1,3
   1  r(i)=x2(i)-x1(i)
      corr=dot(r,q0)
      timcor=corr/v
      return
      end
c
      subroutine norvec(strike,dip,eta)
c
c  calculates the interface unit normal vector, given the layer
c    strike and dip in radians
c
      dimension eta(3)
      sins=sin(strike)
      coss=cos(strike)
      sind=sin(dip)
      cosd=cos(dip)
      eta(1)=sind*sins
      eta(2)=-sind*coss
      eta(3)=cosd
      return
      end
c
      subroutine timdis(dist,q,ii,jj,vel,n,time,
     *                  iface,eta,kk,ll,z,lnumbr,dislyr,deplyr)
c
c   calculates the point a ray, specified by the n ray unit normals
c     given in q, enters the layered medium and its travel-time in
c     the layered system
c
      USE ray3d_par
      !include 'ray3d_param.inc'
      dimension q(3,500000),vel(500000),iface(500000),eta(3,maxlayd),
     *z(maxlayd),dist(3)
      time=0.
c
c   calculates time & dist for the nth to 2nd q-vectors since vector
c     #1 is the incident ray
c
      do 1 i1=1,n-1
         j1=n - i1 + 1
         unum=eta(3,iface(j1))*(z(iface(j1))-dist(3))
     *       -eta(2,iface(j1))*dist(2)
     *       -eta(1,iface(j1))*dist(1)
         u=unum/dot(eta(1,iface(j1)),q(1,j1))
         do 2 i2=1,3
    2       dist(i2)=dist(i2) + u*q(i2,j1)
         time=abs(u)/vel(j1)  + time
         if(iface(j1).eq.lnumbr) then
            dislyr=sqrt(dist(1)**2 + dist(2)**2)
            deplyr=dist(3)
         endif
    1 continue
      if(lnumbr.eq.0) then
          deplyr=dist(3)
          dislyr=sqrt(dist(1)**2 + dist(2)**2)
      endif
      return
      end
c
      subroutine snell(qb,vb,qa,va,itype,sinib,sinia)
c
c   calculates the ray unit normal vector, qa resulting from an
c     incident unit normal vector, qb interacting with a velocity
c     interface.  the medium velocity of qb is vb, the medium
c     velocity of qa is va
c

csed      dimension qb(1),qa(1)
      dimension qb(3),qa(3)

      integer ounit
      common /innout/ inunit,ounit
      torr=float(itype)
      nfix=1
      sinib=sqrt(1.-qb(nfix+2)*qb(nfix+2))
c
c   check for near-vertical incidence.  If sinib < 0.002, then ray is
c     set to true vertical incidence, to avoid instabilities in the
c     calculation of the factor "a" below.  This corresponds to angles
c     of incidence of less than 0.11 degrees, so this manipulation should
c     not cause any significant errors
c***********************************
      if(sinib.le..002) then
         sinib=0.
         qb(nfix+2)=abs(qb(nfix+2))/qb(nfix+2)
         qb(nfix+1)=0.
         qb(nfix)=0.
         sinia=0.
         qa(nfix)=qb(nfix)
         qa(nfix+1)=qb(nfix+1)
         qa(nfix+2)=torr*qb(nfix+2)
         return
      endif
c************************************
c  process all other rays 
c
      sinia=va*sinib/vb
c
c   check for problems with head waves +/or post critically reflected converted
c             phases
c    if any exist, flag the ray and return
c
      if(sinia.ge.1.00) then
c        write(ounit,100) vb,va,qb(nfix+2)
c        if(torr.lt.0..and.qb(nfix+2).lt.0.) write(ounit,101)
c        if(torr.lt.0..and.qb(nfix+2).ge.0.) write(ounit,102)
c        if(torr.gt.0.) write(ounit,103)
         itype = -999
         return
      endif
      if(sinia.lt..0001) then
        a=0.
      else
         a=sinia/sqrt(qb(nfix)*qb(nfix) + qb(nfix+1)*qb(nfix+1))
         qa(nfix)=a*qb(nfix)
         qa(nfix+1)=a*qb(nfix+1)
         qa(nfix+2)=
     +    torr*(qb(nfix+2)/abs(qb(nfix+2)))*sqrt(1. - sinia*sinia)
      endif
      return
!  100 format(' For vb => va of',f6.3,' => ',f6.3,' and qb = ',f6.4)
!  101 format('    ===>  A free surface s-to-p reflection is critical ')
!  102 format('    ===>  An internal s-to-p reflection is critical ')
!  103 format('    ===>  A head wave has been generated')
      end
c
      subroutine wrtray(lyr,az,p,time,baz,p0,pors,init,i4,i10,i11
     *                  ,oldlyr,sini,tdirec)
c
c  writes the results of a ray tracing loop into unit 10
c
c     dimension type(2),wave(4,2),prim(2,2)
      logical pors,emult
c     character type*1,prim*2,wave*3
c     data type/'p','s'/,prim/'pp','ss','ps','sp'/,
c    *     wave/'pmp','pms','smp','sms','sms','smp','pms','pmp'/
      angle=asin(sini)
      angle=angle*57.2957795
      emult=.false.
      if(i4.ne.0) go to 8
         emult=.true.
         go to 1
    8 if(init.ne.0) go to 1
         iprim=1
         if(.not.pors) iprim=2
         if(i4.ne.1) go to 6
            itype=1
            if(.not.pors) itype=2
c           write(9,100) type(itype),baz,p0,tdirec
            t1=0.
c           write(9,102)
c           write(9,101) prim(iprim,i4),t1,az,p,angle
            oldlyr=0
            return
    6    continue
c   6    write(9,103) lyr
c        write(9,102)
c        write(9,104) prim(iprim,i4),lyr,time,az,p,angle
         oldlyr=lyr
         return
    1 ip=1
      if(.not.pors) ip=2
      if(i10.ne.1) goto 2
         if(i11.eq.1) go to 3
            iwave=2
            go to 5
    3       iwave=1
            go to 5
    2 if(i11.eq.1) go to 4
         iwave=4
         go to 5
    4    iwave=3
    5 if(emult) go to 9
      if(lyr.eq.oldlyr) go to 7
c       write(9,103) lyr
c        write(9,102)
         oldlyr=lyr
    7    continue
c   7 write(9,105) prim(ip,i4),wave(iwave,ip),lyr,time,az,p,angle
      return
c   9 if(iwave.eq.1.and.ip.eq.1) write(9,106) lyr
    9 continue
c     write(9,107) wave(iwave,ip),time,az,p,angle
      return
c  100 format(///' incident ',a1,'-wave, back azimuth: ',f6.2,
c    *       ' ray parameter: ',f7.4,/,' direct arrival spends ',f7.3,
c    *       ' secs in layering',/,' all times relative to direct ray'
c    *       ,/)
c 101 format(5x,a2,5x,'direct',2x,f7.3,3x,f7.2,7x,f7.4,6x,f5.2)
c 102 format(' wave type   layer    time     azimuth     ray param.',
c    *       '   angle')
c 103 format(1x,/,' layer ',i2)
c 104 format(5x,a2,7x,i2,4x,f7.3,3x,f7.2,7x,f7.4,6x,f5.2)
c 105 format(3x,a2,a3,6x,i2,4x,f7.3,3x,f7.2,7x,f7.4,6x,f5.2)
c 106 format(63x,'extra multiples from layer ',i2,/,
c    *       63x,' type    time        az.          p         angle')
c 107 format(64x,a3,3x,f7.3,4x,f7.2,6x,f7.4,6x,f5.2)
      end
c
      subroutine ampcal(amagb,hmagb,amaga,hmaga,strike,dip,ihilb)
c
c subroutine to calculate amplitudes for rays from ray3d
c
c   i n p u t
c
c
csed      dimension qb(3),qa(3),amagb(1),amaga(1),r3(3),rt(3),at(3),ai(3),
csed     *          a(3,3),hmaga(1),hmagb(1),rth(3),ht(3)
      dimension qb(3),qa(3),amagb(3),amaga(3),r3(3),rt(3),at(3),ai(3),
     *          a(3,3),hmaga(3),hmagb(3),rth(3),ht(3)
      logical free
      integer type,ounit
      common /cord/ a
      common /amcal/ qb,qa,vb,va,sinib,sinia,vp1,vs1,rho1,
     *               vp2,vs2,rho2,free,type
      common /innout/ inunit,ounit
      nfix=1
!      call zeerow(r3,1,3)
!      call zeerow(rt,1,3)
!      call zeerow(at,1,3)
!      call zeerow(ai,1,3)
!      call zeerow(ht,1,3)
!      call zeerow(rth,1,3)
      r3 = 0.
      rt = 0.
      at = 0.
      ai = 0.
      ht = 0.
      rth = 0.
      rshph=0.
      rph=0.
      rphx=0.
      rphy=0.
      rphz=0.
      rmag=0.
      ncode=0
      eps=.0001
      ihilb=0
      rshmag=0.
      pi=3.14159
c
c vertical incidence case requires special treatment
c
      if(sinib.gt.eps) then
         cosphi=-qb(1)/sinib
         sinphi=-qb(2)/sinib
      else
         cosphi=-1.
         sinphi=0.
      endif
      nd=0
      if(abs(qb(3))/qb(3).gt.0) nd=1
      p=sinib/vb
      if(free) go to 10
      ro2=rho2
c
c   find ncode for non-free surface case
c
      if(abs(vb-vp1).gt.eps) go to 1
        call rcomp(ai,1,nd,sinib,.true.)
        if(type.lt.0) go to 2
        if(abs(va-vp2).lt.eps) ncode=3
        if(abs(va-vs2).lt.eps) ncode=4
        go to 3
    2   if(abs(va-vp1).lt.eps) ncode=1
        if(abs(va-vs1).lt.eps) ncode=2
        go to 3
    1 if(abs(vb-vs1).gt.eps) go to 4
        call rcomp(ai,2,nd,sinib,.true.)
        if(type.lt.0) go to 5
        if(abs(va-vs2).lt.eps) ncode=8
        if(abs(va-vp2).lt.eps) ncode=7
        go to 3
    5   if(abs(va-vp1).lt.eps) ncode=5
        if(abs(va-vs1).lt.eps) ncode=6
    3 ncase=0
      if(ncode.eq.0) go to 4
      if(ncode.le.4) go to 7
         ncase=4
         go to 7
c
c  find ncode for free surface case
c
   10 ro2=0.0
      vp2=0.
      vs2=0.
      if(type.eq.0) go to 15
      if(abs(vb-vp1).gt.eps) go to 12
         call rcomp(ai,1,nd,sinib,.true.)
         if(abs(va-vs1).lt.eps) ncode=2
         if(abs(va-vp1).lt.eps) ncode=1
         go to 13
   12    if(abs(vb-vs1).gt.eps) go to 4
         call rcomp(ai,2,nd,sinib,.true.)
         if(abs(va-vs1).lt.eps) ncode=4
         if(abs(va-vp1).lt.eps) ncode=3
   13 ncase=0
      if(ncode.eq.0) go to 4
      if(ncode.le.2) go to 7
         ncase=2
         go to 7
c
c
c  f i n d  f r e e  s u r f a c e  e f f e c t
c
c
   15 if(abs(vb-vp1).lt.eps) go to 16
      if(abs(vb-vs1).lt.eps) go to 17
      go to 4
   16 call rcomp(ai,1,nd,sinib,.true.)
      call coef8(p,vp1,vs1,rho1,vp2,vs2,0.0,5,nd,rx,rphx)
      call coef8(p,vp1,vs1,rho1,vp2,vs2,0.0,6,nd,rz,rphz)
      ry=0.
      rphy=0.
      go to 18
   17 call rcomp(ai,2,nd,sinib,.true.)
      call coef8(p,vp1,vs1,rho1,vp2,vs2,0.0,7,nd,rx,rphx)
      call coef8(p,vp1,vs1,rho1,vp2,vs2,0.0,8,nd,rz,rphz)
      call coefsh(p,vs1,rho1,vs2,0.0,2,ry,rphy)
   18 if(abs(rphx+pi).gt.eps) go to 22
        rphx=0.
        rx=-rx
   22 if(abs(rphy+pi).gt.eps) go to 23
        rphy=0.
        ry=-ry
   23 if(abs(rphz+pi).gt.eps) go to 24
        rphz=0.
        rz=-rz
   24 do 19 i19=1,3
         rth(i19)=hmagb(i19)
   19    rt(i19)=amagb(i19)
c
c  rt is in global coordinates, but this is equivalent to interface
c    coordinates for the free surface. so transform rt directly to
c    the ray coordinate system
c
      call rtoi(rt,cosphi,sinphi,qb(3),.false.)
      call rtoi(rth,cosphi,sinphi,qb(3),.false.)
      phck=0.
      phck=abs(rphz)+abs(rphx)+abs(rphy)
      if(phck.gt.eps) ihilb=1
      dotar=dot(ai,rt)
c
c vertical incidence can sometimes blow up this step
c    check first
c
      if(abs(dotar).lt.eps) go to 56
      dotar=abs(dotar)/dotar
   56 doth=dot(ai,rth)
      if(abs(doth).lt.eps) go to 26
      doth=abs(doth)/doth
   26 amh=sqrt(rth(1)*rth(1) + rth(3)*rth(3))*doth
      amb=sqrt(rt(1)*rt(1) + rt(3)*rt(3))*dotar
      amaga(nfix)=rx*(amb*cos(rphx) - amh*sin(rphx))
      amaga(nfix+1)=ry*(rt(2)*cos(rphy) - rth(2)*sin(rphy))
      amaga(nfix+2)=rz*(amb*cos(rphz) - amh*sin(rphz))
      hmaga(nfix)=rx*(amh*cos(rphx) + amb*sin(rphx))
      hmaga(nfix+1)=ry*(rth(2)*cos(rphy) + rt(2)*sin(rphy))
      hmaga(nfix+2)=rz*(amh*cos(rphz) + amb*sin(rphz))
      call rtoi(amaga,cosphi,sinphi,qb(3),.true.)
      call rtoi(hmaga,cosphi,sinphi,qb(3),.true.)
      return
c
c
c  g e n e r a l  c o e f i c i e n t  c a l c u l a t i o n
c
c  first find rt, the incident displacement vector in ray coordinates
c        &    rth, the distorted displacement vector in ray coordinates
c
    7 call coord(amagb,strike,dip,rt,'local',.true.)
      call coord(hmagb,strike,dip,rth,'local',.true.)

      call rtoi(rt,cosphi,sinphi,qb(3),.false.)
      call rtoi(rth,cosphi,sinphi,qb(3),.false.)
      call coef8(p,vp1,vs1,rho1,vp2,vs2,ro2,ncode,nd,rmag,rph)
      call rcomp(r3,ncode-ncase,nd,sinia,.false.)
      if(abs(rph + pi).gt.eps) go to 20
         rph=0.
         rmag=-rmag
   20 at(2)=0.0
c
c  if incident & resulting waves are both s-waves, find sh coeficient
c
      if(ncode.le.4) go to 9
      if(ncode.eq.6) then
          ncodsh=1
      elseif(ncode.eq.8) then
          ncodsh=2
      else
          go to 9
      endif
      call coefsh(p,vs1,rho1,vs2,ro2,ncodsh,rshmag,rshph)
      at(2)=rshmag*(rt(2)*cos(rshph)-rth(2)*sin(rshph))
      ht(2)=rshmag*(rth(2)*cos(rshph)-rt(2)*sin(rshph))
      if(abs(rshph+pi).lt.eps) go to 9
      if(rshph.gt.eps) ihilb=1
    9 dotar=dot(ai,rt)
c
c vertical incidence can sometimes blow up this step
c    check first
c
      if(abs(dotar).lt.eps) go to 55
      dotar=abs(dotar)/dotar
   55 amb=sqrt(rt(1)*rt(1) + rt(3)*rt(3))*dotar
      doth=dot(ai,rth)
      if(abs(doth).lt.eps) go to 25
      doth=abs(doth)/doth
   25 amh=sqrt(rth(1)*rth(1) + rth(3)*rth(3))*doth
      atmag=rmag*(amb*cos(rph)-amh*sin(rph))
      htmag=rmag*(amh*cos(rph)+amb*sin(rph))
      if(rph.gt.eps) ihilb=1
      at(1)=atmag*r3(1)
      at(3)=atmag*r3(3)
      ht(1)=htmag*r3(1)
      ht(3)=htmag*r3(3)
      call rtoi(at,cosphi,sinphi,qb(3),.true.)
      call rtoi(ht,cosphi,sinphi,qb(3),.true.)
      call coord(at,strike,dip,amaga,'globe',.true.)
      call coord(ht,strike,dip,hmaga,'globe',.true.)
      return
c    4 write(ounit,102) va,vb,vp1,vs1,vp2,vs2
c  102 format(' ncode = 0 for ',6f6.2)
   4  continue
      return
      end
c
      subroutine rtoi(r,cosp,sinp,qb,dirtcn)
c
c   transforms a vector r from the ray coordinate system
c     to the interface coordinate system and vice versa
c
c   if dirtcn = .true.  ray => interface
c      dirtcn = .false. interface => ray
c
c   qb is the z component of the ray in the interface system
c

csed      dimension r(1)
      dimension r(3)

      logical dirtcn
      nfix=1
      q=abs(qb)/qb
      r(nfix+2)=r(nfix+2)*(-q)
      if(dirtcn) go to 1
      xr=r(nfix)*cosp + r(nfix+1)*sinp
      yr=r(nfix)*sinp - r(nfix+1)*cosp
      r(nfix)=xr*q
      r(nfix+1)=yr
      return
    1 xr=r(nfix)*q
      xl=+xr*cosp + r(nfix+1)*sinp
      yl= xr*sinp - r(nfix+1)*cosp
      r(nfix)=xl
      r(nfix+1)=yl
      return
      end
c
      subroutine rcomp(r3,ncode,nd,sini,incdnt)
c
c   resolves a reflection coeficient r from s/r coef8 into
c     x and z components (in the ray coordinate system)
c     given the resulting ray type:
c       reflected p => ncode = 1
c       reflected s => ncode = 2
c       transmitted p => ncode = 3
c       transmitted s => ncode = 4
c
      dimension r3(3)
      logical incdnt
      cosi=sqrt(1. - sini*sini)
      r3(2)=0.
      if(incdnt) go to 10
      if(nd.ne.0) go to 5
      go to (1,2,3,4) ncode
    1 r3(3)=cosi
      r3(1)=sini
      return
    2 r3(3)=sini
      r3(1)=-cosi
      return
    3 r3(3)=-cosi
      r3(1)=sini
      return
    4 r3(3)=sini
      r3(1)=cosi
      return
    5 go to (6,7,8,9) ncode
    6 r3(3)=cosi
      r3(1)=-sini
      return
    7 r3(3)=-sini
      r3(1)=-cosi
      return
    8 r3(3)=-cosi
      r3(1)=-sini
      return
    9 r3(3)=-sini
      r3(1)=cosi
      return
   10 if(nd.ne.0) go to 11
      go to (3,4) ncode
   11 go to (8,9) ncode
      return
      end
c
      subroutine rayfin(iq,i4,i10,i11,lnumbr,dflag,mflag)
      USE ray3d_par
      !include 'ray3d_param.inc'
      dimension strike(maxlayd),dip(maxlayd),z(maxlayd),alpha(maxlayd),
     *beta(maxlayd),
     *eta(3,maxlayd),q(3,500000),q0(3),v(2,maxlayd),qv(500000),
     *raydis(3),qloc1(3),qloc2(3),a(3,3),iface(500000),layer(maxlayd),
     *amag(3,500000),hmag(3,500000),raymag(3,500000),rayhil(3,500000),
     *raytim(500000),direct(3),rho(maxlayd)
      logical pors,amps,free,dflag,mflag
      integer trans,refl,type
      common /cord/ a
      common /amcal/ qloc1,qloc2,vb,va,sinib,sinia,vp1,vs1,rho1,
     *               vp2,vs2,rho2,free,type
      common /transm/ q,qv,v,alpha,beta,rho,strike,dip,iface,jhilb,
     *                amag,hmag,layer,amps,trans,refl,nlyrs
      common /raywrt/ eta,z,raymag,rayhil,raytim,ntim,p0,pors,
     *                oldlyr,q0,direct,tdirec,baz
      !call zeerow(raydis,1,3)
      raydis = 0.
      call timdis(raydis,q,3,500000,qv,iq,time,
     *            iface,eta,3,maxlayd,z,lnumbr,dislyr,deplyr)
      tmpdis=sqrt(raydis(1)**2 + raydis(2)**2)
      init=1
      if(.not.dflag) go to 1
         tdirec=time
         time=0.
         init=0.
         do 3 i3=1,3
    3    direct(i3)=raydis(i3)
         go to 2
    1 time=time + timcor(direct,raydis,q0,v(1,nlyrs))-tdirec
      if(.not.dflag.and..not.mflag) init=0
    2 call anom(q(1,iq),qv(iq),azanom,panom,angle)
      if(.not.amps) go to 26
         free=.true.
         type=0
         do 52 i52=1,3
   52    qloc1(i52)=q(i52,iq)
         vb=qv(iq)
         sinib=angle
         va=0.
         sinia=0.
         vp1=alpha(1)
         vs1=beta(1)
         rho1=rho(1)
         rho2=0.
         vp2=0.
         vs2=0.
         call ampcal(amag(1,iq),hmag(1,iq),
     *               raymag(1,ntim),rayhil(1,ntim),
     *               0.,0.,ihilb)
         if(ihilb.eq.1) jhilb=1
         free=.false.
         raytim(ntim)=time
         ntim=ntim+1
   26 continue
c     write(67,*) ntim,lnumbr,dislyr,deplyr
c     call wrtray(lnumbr,azanom,panom,time,baz,p0,pors,
c    *            init,i4,i10,i11,oldlyr,angle,tdirec)
      return
      end



c
      subroutine trnsmt(loopst,looped,iq,iv,up)
c
c *******************
c
c     calculates the amplitude of a wave transmitted through
c     a stack of layers
c
c *******************
      USE ray3d_par
      !include 'ray3d_param.inc'
      dimension strike(maxlayd),dip(maxlayd),alpha(maxlayd),
     *beta(maxlayd),rho(maxlayd),q(3,500000),v(2,maxlayd),qv(500000),
     *qloc1(3),qloc2(3),a(3,3),iface(500000),layer(maxlayd),
     *amag(3,500000),hmag(3,500000)
      logical amps,free,up
      integer trans,refl,type,itype
      common /cord/ a
      common /amcal/ qloc1,qloc2,vb,va,sinib,sinia,vp1,vs1,rho1,
     *               vp2,vs2,rho2,free,type
      common /transm/ q,qv,v,alpha,beta,rho,strike,dip,iface,jhilb,
     *                amag,hmag,layer,amps,trans,refl,nlyrs
      do 7 i7=loopst,looped-1
         j7=looped - i7 + 1
         if(.not.up) j7=i7
         k7=j7-1
         if(.not.up) k7=k7+1
        call coord(q(1,iq),strike(j7),dip(j7),qloc1,'local',
     *              .false.)


         vb=qv(iq)
         va=v(iv,k7)
         itype=trans
         call snell(qloc1,vb,qloc2,va,itype,sinib,sinia)
c
c   if itype returns as -999, then a problem phase exists
c      iq is flagged for return to main program -- ray will be skipped
c
         if(itype.eq.-999) then
            iq=-999
            return
         endif
         if(.not.amps) go to 19
            vp2=alpha(k7)
            vs2=beta(k7)
            rho2=rho(k7)
            type=trans
            call ampcal(amag(1,iq),hmag(1,iq),
     *                  amag(1,iq+1),hmag(1,iq+1)
     *                 ,strike(j7),dip(j7),ihilb)
            if(ihilb.eq.1) jhilb=1
            vp1=vp2
            vs1=vs2
            rho1=rho2
   19    qv(iq+1)=va
         call coord(qloc2,strike(j7),dip(j7),q(1,iq+1),'globe',
     *              .true.)
         iq=iq+1
         iface(iq)=j7
    7 continue
      return
      end
c
      subroutine raydwn(iqref,i10,lyref,iq)
c
c **************
c
c     subroutine to reflect a ray from the free surface then
c                propagate it down to a designated interface
c
c **************
      USE ray3d_par
      !include 'ray3d_param.inc'
      dimension strike(maxlayd),dip(maxlayd),alpha(maxlayd),
     *beta(maxlayd),rho(maxlayd),q(3,500000),v(2,maxlayd),qv(500000),
     *qloc1(3),qloc2(3),a(3,3),iface(500000),layer(maxlayd),
     *amag(3,500000),hmag(3,500000)
      logical amps,free
      integer trans,refl,type,itype
      common /cord/ a
      common /amcal/ qloc1,qloc2,vb,va,sinib,sinia,vp1,vs1,rho1,
     *               vp2,vs2,rho2,free,type
      common /transm/ q,qv,v,alpha,beta,rho,strike,dip,iface,jhilb,
     *                amag,hmag,layer,amps,trans,refl,nlyrs
      iq=iqref
c
c  take ray down to the reflecting interface --
c
c   do reflection from free surface first
c
      call coord(q(1,iq),strike(1),dip(1),qloc1,'local',
     *           .false.)
      vb=qv(iq)
      va=v(i10,1)
      type=refl
      itype=type
      call snell(qloc1,vb,qloc2,va,itype,sinib,sinia)
c
c   if itype returns as -999, then a problem phase exists
c      iq is flagged for return to main program -- ray will be skipped
c
      if(itype.eq.-999) then
         iq=-999
         return
      endif
      if(.not.amps) go to 20
         free=.true.
         call ampcal(amag(1,iq),hmag(1,iq),
     *               amag(1,iq+1),hmag(1,iq+1),
     *               strike(1),dip(1),ihilb)
         if(ihilb.eq.1) jhilb=1
         free=.false.
   20    qv(iq+1)=va
         call coord(qloc2,strike(1),dip(1),q(1,iq+1),'globe'
     *              ,.true.)
         iq=iq+1
         iface(iq)=1
c
c   now transmit wave down to reflecting interface
c
      if(lyref.eq.2) return
c
c  iq could be returned as -999 from s/r trnsmt -- ray would be skipped
c
      call trnsmt(2,lyref,iq,i10,.false.)
      return
      end
c
      subroutine rayup(iqref,i11,lyref,iq)
c
c **************
c
c     subroutine to reflect a ray off an interface at depth then
c                transmit it back up to the free surface
c
c **************
      USE ray3d_par
      !include 'ray3d_param.inc'
      dimension strike(maxlayd),dip(maxlayd),alpha(maxlayd),
     *beta(maxlayd),rho(maxlayd),q(3,500000),v(2,maxlayd),qv(500000),
     *qloc1(3),qloc2(3),a(3,3),iface(500000),layer(maxlayd),
     *amag(3,500000),hmag(3,500000)
      logical amps,free
      integer trans,refl,type,itype
      common /cord/ a
      common /amcal/ qloc1,qloc2,vb,va,sinib,sinia,vp1,vs1,rho1,
     *               vp2,vs2,rho2,free,type
      common /transm/ q,qv,v,alpha,beta,rho,strike,dip,iface,jhilb,
     *                amag,hmag,layer,amps,trans,refl,nlyrs
      iq=iqref
      vp1=alpha(lyref-1)
      vs1=beta(lyref-1)
      rho1=rho(lyref-1)
c
c  do the reflection off the interface first
c
      j12=lyref
      call coord(q(1,iq),strike(j12),dip(j12),qloc1,
     *           'local',.false.)
      vb=qv(iq)
      va=v(i11,j12-1)
      type=refl
      itype=type
      call snell(qloc1,vb,qloc2,va,itype,sinib,sinia)
c
c   if itype returns as -999, then a problem phase exists
c      iq is flagged for return to main program -- ray will be skipped
c
      if(itype.eq.-999) then
         iq=-999
         return
      endif
      if(.not.amps) go to 22
         vp2=alpha(j12)
         vs2=beta(j12)
         rho2=rho(j12)
         call ampcal(amag(1,iq),hmag(1,iq),
     *               amag(1,iq+1),hmag(1,iq+1),
     *               strike(j12),dip(j12),ihilb)
         if(ihilb.eq.1) jhilb=1
   22 qv(iq+1)=va
      call coord(qloc2,strike(j12),dip(j12),q(1,iq+1),
     *           'globe',.true.)
      iq=iq+1
      iface(iq)=j12
c
c now transmit wave back to surface
c
      if(lyref.eq.2) return
c
c   iq could be returned as -999 from s/r trnsmt -- ray would be skipped
c
      call trnsmt(2,lyref,iq,i11,.true.)
      return
      end
      subroutine coef8(p,vp1,vs1,ro1,vp2,vs2,ro2,ncode,nd,rmod,rph)
c
c     the routine coef8 is designed for the computation of reflection
c     and transmission coefficients at a plane interface between two
c     homogeneous solid halfspaces or at a free surface of a homogeneous
c     solid halfspace.
c
c     the codes of individual coefficients are specified by the
c     following numbers
c     a/ interface between two solid halfspaces
c     p1p1...1       p1s1...2       p1p2...3       p1s2...4
c     s1p1...5       s1s1...6       s1p2...7       s1s2...8
c     b/ free surface (for ro2.lt.0.00001)
c     pp.....1       px.....5       px,pz...x- and z- components of the
c     ps.....2       pz.....6       coef.of conversion,incident p wave
c     sp.....3       sx.....7       sx,sz...x- and z- components of the
c     ss.....4       sz.....8       coef.of conversion,incident s wave
c
c     i n p u t   p a r a m e t e r s
c           p...ray parameter
c           vp1,vs1,ro1...parameters of the first halfspace
c           vp2,vs2,ro2...parameters of second halfspace. for the free
c                    surface take ro2.lt.0.00001,eg.ro2=0., and
c                    arbitrary values of vp2 and vs2
c           ncode...code of the computed coefficient
c           nd...=0  when the positive direction of the ray
c                    and the x-axis make an acute angle
c                =1  when the wave impinges on the interface
c                    against the positive direction of the x-axis
c
c     o u t p u t   p a r a m e t e r s
c           rmod,rph...modul and argument of the coefficient
c
c     n o t e s
c     1/ positive p...in the direction of propagation
c     2/ positive s...to the left from p
c     3/ time factor of incident wave ... exp(-i*omega*t)
c     4/ formulae are taken from cerveny ,molotkov, psencik, ray method
c        in seismology, pages 30-35. due to the note 2, the signs at
c        certain coefficients are opposite
c
c       written by v.cerveny,1976
c       modified by t.j. owens, 3/22/82 see comments in code
c
      complex b(4),rr,c1,c2,c3,c4,h1,h2,h3,h4,h5,h6,h,hb,hc
      dimension prmt(4),d(4),dd(4)
c
      if(ro2.lt.0.000001)go to 150
      prmt(1)=vp1
      prmt(2)=vs1
      prmt(3)=vp2
      prmt(4)=vs2
      a1=vp1*vs1
      a2=vp2*vs2
      a3=vp1*ro1
      a4=vp2*ro2
      a5=vs1*ro1
      a6=vs2*ro2
      q=2.*(a6*vs2-a5*vs1)
      pp=p*p
      qp=q*pp
      x=ro2-qp
      y=ro1+qp
      z=ro2-ro1-qp

      g1=a1*a2*pp*z*z
      g2=a2*x*x
      g3=a1*y*y
      g4=a4*a5
      g5=a3*a6
      g6=q*q*pp
      do 21 i=1,4
      dd(i)=p*prmt(i)
   21 d(i)=sqrt(abs(1.-dd(i)*dd(i)))

      if(dd(1).le.1..and.dd(2).le.1..and.dd(3).le.1..and.dd(4).le.1.)
     1go to 100
c
c     complex coefficients
      do 22 i=1,4
      if(dd(i).gt.1.)go to 23
      b(i)=cmplx(d(i),0.)
      go to 22
   23 b(i)= cmplx(0.,d(i))
   22 continue
      c1=b(1)*b(2)
      c2=b(3)*b(4)
      c3=b(1)*b(4)
      c4=b(2)*b(3)
      h1=g1
      h2=g2*c1
      h3=g3*c2
      h4=g4*c3
      h5=g5*c4
      h6=g6*c1*c2
      h=1./(h1+h2+h3+h4+h5+h6)
      hb=2.*h
      hc=hb*p
      go to (1,2,3,4,5,6,7,8),ncode
    1 rr=h*(h2+h4+h6-h1-h3-h5)
      go to 26
    2 rr=vp1*b(1)*hc*(q*y*c2+a2*x*z)
      if(nd.ne.0)rr=-rr
      go to 26
    3 rr=a3*b(1)*hb*(vs2*b(2)*x+vs1*b(4)*y)
      go to 26
    4 rr=-a3*b(1)*hc*(q*c4-vs1*vp2*z)
      if(nd.ne.0)rr=-rr
      go to 26
    5 rr=-vs1*b(2)*hc*(q*y*c2+a2*x*z)
      if(nd.ne.0)rr=-rr
      go to 26
    6 rr=h*(h2+h5+h6-h1-h3-h4)
      go to 26
    7 rr=a5*b(2)*hc*(q*c3-vp1*vs2*z)
      if(nd.ne.0)rr=-rr
      go to 26
    8 rr=a5*b(2)*hb*(vp1*b(3)*y+vp2*b(1)*x)
      go to 26
c     real coefficients
  100 e1=d(1)*d(2)
      e2=d(3)*d(4)
      e3=d(1)*d(4)
      e4=d(2)*d(3)
      s1=g1
      s2=g2*e1
      s3=g3*e2
      s4=g4*e3
      s5=g5*e4
      s6=g6*e1*e2
      s=1./(s1+s2+s3+s4+s5+s6)
      sb=2.*s
      sc=sb*p
      go to (101,102,103,104,105,106,107,108),ncode
  101 r=s*(s2+s4+s6-s1-s3-s5)
      go to 250
  102 r=vp1*d(1)*sc*(q*y*e2+a2*x*z)
      if(nd.ne.0)r=-r
      go to 250
  103 r=a3*d(1)*sb*(vs2*d(2)*x+vs1*d(4)*y)
      go to 250
  104 r=-a3*d(1)*sc*(q*e4-vs1*vp2*z)
      if(nd.ne.0)r=-r
      go to 250
  105 r=-vs1*d(2)*sc*(q*y*e2+a2*x*z)
      if(nd.ne.0)r=-r
      go to 250
  106 r=s*(s2+s5+s6-s1-s3-s4)
      go to 250
  107 r=a5*d(2)*sc*(q*e3-vp1*vs2*z)
      if(nd.ne.0)r=-r
      go to 250
  108 r=a5*d(2)*sb*(vp1*d(3)*y+vp2*d(1)*x)
      go to 250
c
c     earths surface,complex coefficients and coefficients of conversion
c
c   n o t e :
c
c    signs of coefficients at loops 162, 166, & 168 have been changed
c    from the originnal version of coef8 due to inconsistencies in
c    notation from the cerveny, et al book
c    3/22/82
c
  150 a1=vs1*p
      a2=a1*a1
      a3=2.*a2
      a4=2.*a1
      a5=a4+a4
      a6=1.-a3
      a7=2.*a6
      a8=2.*a3*vs1/vp1
      a9=a6*a6
      dd(1)=p*vp1
      dd(2)=p*vs1
      do 151 i=1,2
  151 d(i)=sqrt(abs(1.-dd(i)*dd(i)))
      if(dd(1).le.1..and.dd(2).le.1.)go to 200
      do 154 i=1,2
      if(dd(i).gt.1.)go to 155
      b(i)=cmplx(d(i),0.)
      go to 154
  155 b(i)= cmplx(0.,d(i))
  154 continue
      h1=b(1)*b(2)
      h2=h1*a8
      h=1./(a9+h2)
      go to (161,162,163,164,165,166,167,168),ncode
  161 rr=(-a9+h2)*h
      go to 26
  162 rr=-a5*b(1)*h*a6
      if(nd.ne.0)rr=-rr
      go to 26
  163 rr=a5*b(2)*h*a6*vs1/vp1
      if(nd.ne.0)rr=-rr
      go to 26
  164 rr=-(a9-h2)*h
      go to 26
  165 rr=a5*h1*h
      if(nd.ne.0)rr=-rr
      go to 26
  166 rr=-a7*b(1)*h
      go to 26
  167 rr=a7*b(2)*h
      go to 26
  168 rr=a5*vs1*h1*h/vp1
      if(nd.ne.0)rr=-rr
   26 z2=real(rr)
      z3=aimag(rr)
      if(z2.eq.0..and.z3.eq.0.)go to 157
      rmod=sqrt(z2*z2+z3*z3)
      rph=atan2(z3,z2)
      return
  157 rmod=0.
      rph=0.
      return
c
c     earths surface,real coefficients and coefficients of conversion
c   n o t e :
c
c    signs of coeficients at loops 202, 206, & 208 have been reversed
c    by t.j. owens because of inconsistencies w/sign conventions
c    3/22/82
c
  200 s1=d(1)*d(2)
      s2=a8*s1
      s=1./(a9+s2)
      go to (201,202,203,204,205,206,207,208),ncode
  201 r=(-a9+s2)*s
      go to 250
  202 r=-a5*d(1)*s*a6
      if(nd.ne.0)r=-r
      go to 250
  203 r=a5*d(2)*s*a6*vs1/vp1
      if(nd.ne.0)r=-r
      go to 250
  204 r=(s2-a9)*s
      go to 250
  205 r=a5*s1*s
      if(nd.ne.0)r=-r
      go to 250
  206 r=-a7*d(1)*s
      go to 250
  207 r=a7*d(2)*s
      go to 250
  208 r=a5*vs1*s1*s/vp1
      if(nd.ne.0)r=-r
  250 if(r.lt.0.)go to 251
      rmod=r
      rph=0.
      return
  251 rmod=-r
      rph=-3.14159
      return
      end

      subroutine coefsh(p,vs1,rho1,vs2,rho2,ncode,rmod,rph)
c
c  calculates sh-wave reflection and transmission coeficients
c    at a solid-solid interface of a free surface
c
c   solid-solid:  ncode =>
c
c    s1s1 = 1   s1s2 = 2
c
c   free surface: rho2 = 0.0 ncode =>
c
c    s1s1 = 1   free surface correction = 2
c
      complex p2,p4,d,h1,h2,rr
      if(rho2.lt..00001) go to 5
      a1=vs1*p
      a2=vs2*p
      b1=rho1*vs1
      b2=rho2*vs2
      g1=sqrt(abs(1. - a1*a1))
      g2=sqrt(abs(1. - a2*a2))
      p2=cmplx(g1,0.)
      p4=cmplx(g2,0.)
      if(a1.gt.1.) p2=cmplx(0.,g1)
      if(a2.gt.1.) p4=cmplx(0.,g2)
      h1=cmplx(b1,0.)*p2
      h2=cmplx(b2,0.)*p4
      d= h1 + h2
      go to (1,2) ncode
    1 rr=(h1-h2)/d
      go to 3
    2 rr=2.*h1/d
    3 z1=real(rr)
      z2=aimag(rr)
      if(z2.eq.0.) go to 4
      rmod=sqrt(z1*z1 + z2*z2)
      rph=atan2(z2,z1)
      return
    4 rmod=z1
      rph=0.
      return
c
c    free surface problem
c
    5 go to (6,7) ncode
    6 rmod=1.
      rph=0.
      return
    7 rmod=2.
      rph=0.
      return
      end

      function dot(x,y)
c
c  calculates the dot product of two vectors
c
csed      dimension x(1),y(1)
      dimension x(3),y(3)
      z=0.
      do 1 i=1,3
    1    z=z + x(i)*y(i)
      dot=z
      return
      end

        subroutine pwaveqn(npts,dt,r1,z1,recv,dely,c,agauss)
c
c  ********************************************************************
c
c    fortran 77 program to perform a source equalization deconvolution
c     on a three component seismogram - using the method of langston (1979)
c
c   f77 -g -o pwaveqn pwaveqn.f ../Subs/subs.a /usr/local/sac/lib/sac.a
c
c  *************************************************************************
c
c     parameter(MAXPOINTS=10000, MAXPOINTS2=MAXPOINTS*2)
c     dimension d2(100000),caz(3),recv(1),r1(1),z1(1)


csed2      dimension d2(100000),recv(1),r1(1),z1(1)
      dimension d2(100000),recv(2*npts),r1(npts),z1(npts)


c     character eqfile*32,outfil*32,comp(3,2)*6,outc(2)*5,knm(2)*8
c     character comp(3,2)*6,outc(2)*5
      character knm(2)*8
      complex dat(100000,2)
c     double precision gnorm
      integer ounit
c **********************************************************************
c
c common block info for link with subroutine sacio
c
c   scaio is in Tom Owens' SAC input/output routines in Subs library
c	
c     real instr
c     integer year,jday,hour,min,isec,msec
      character*8 cmpnm
c
c **************************************************************************
c
c   parameter definitions may be found in sacio comments
c
      common /innout/ inunit,ounit
c     data
c    *     comp/'.z    ','.n    ','.e    ','_sp.z','_sp.r','_sp.t'/
c    *    ,outc/'.eqr ','.eqt '/
      data knm/'radial  ','tangentl'/
      inunit=5
      ounit=6
      zero1=cmplx(0.,0.)
      pi=3.141592654
c     call iniocm
      do 1 i=1,2
         call zer(dat(1,i),1,10000)
    1 continue
      do 88 i=1,npts/2
        i2=2*i
        i1=i2-1
        dat(i,1)=cmplx(z1(i1),z1(i2))
 88     dat(i,2)=cmplx(r1(i1),r1(i2))
c
c    *******************************************************************
c
      nft=npowr2(npts)
      nfpts=nft/2 + 1
      fny=1./(2.*dt)
      delf=fny/float(nft/2)
c     write(ounit,102) npts,nft,fny,delf,dt
c 102 format(1x,'npts=',i5,1x,'nft=',i5,1x,'fny=',f7.4,1x,
c    *          'delf=',f8.4,1x,'dt=',f6.3)
c
c  change delf to normalize deconvolution so dfftr works
c  dt from forward step is cancelled in decon, so no delf
c  on inverse, just 1/nft
c
      cdelf = 1. / float( nft )
      do 5 i=1,2
         call dfftr(dat(1,i),nft,'forward',dt)
         if(i.ne.1) go to 5
         d2max=0.
         do 7 j=1,nfpts
            d2(j)=real(dat(j,i)*conjg(dat(j,i)))
            if(d2(j).gt.d2max) d2max=d2(j)
    7    continue
    5 continue
      decon=1.
c
c ... Decon parameters
c
      tdelay=dely
      t0=dely
c
      phi1=c*d2max
      gnorm = 0.0d0
c      gnorm = 0.0
c
      do 8 i=1,1
         do 9 j=1,nfpts
            freq=float(j-1)*delf
            w=2.*pi*freq
            phi=phi1
            if(d2(j).gt.phi) phi=d2(j)
            gauss=-w*w/(4.*agauss*agauss)
            gnorm = gnorm + exp(gauss)
            dat(j,i+1)=dat(j,i+1)*conjg(dat(j,1))*
     *                   cmplx(exp(gauss)/phi,0.)
            dat(j,i+1)=dat(j,i+1)*exp(cmplx(0.,-w*tdelay))
    9    continue
         call dfftr(dat(1,i+1),nft,'inverse',cdelf)
    8 continue
c
c     deconvolve the vertical from itself using the
c         specified water-level parameter and gaussian
c
      do 19 j=1,nfpts
         freq=float(j-1)*delf
         w=2.*pi*freq
         phi=phi1
         gauss=-w*w/(4.*agauss*agauss)
         if(d2(j).gt.phi) phi=d2(j)
         dat(j,1)=dat(j,1)*conjg(dat(j,1))*
     *                   cmplx(exp(gauss)/phi,0.)
         dat(j,1)=dat(j,1)*exp(cmplx(0.,-w*tdelay))
19    continue
c
c*************************************************************
c
c     inverse transform the equalized vertical component
c
      call dfftr(dat(1,1),nft,'inverse',cdelf)

c
c     compute the maximum value of the vertical component
c        to include it in the normalization factor gnorm
c
      call minmax(dat(1,1),npts,dmin,dmax,dmean)
c
c     normalize the dmax for the transforms and gaussian
c     Not really necessary, horizontals have the same factors.
c
c     dmax = dmax *  float(nfpts) / gnorm
c
c************************************************************* 
c
c     gnorm = dmax * gnorm / nfpts
c
c     note that (dmax * gnorm / nfpts) = the unormalized dmax
c
      gnorm = dmax
      !gnorm =  gnorm / nfpts !JSMALE: I changed here
      do 111 j = 1,npts/2
       dat(j,2) = dat(j,2) / gnorm
       !dat(j,2) = dat(j,2)
 111  continue
           az=0.0
           cinc=90.
           cmpnm=knm(1)
c          call minmax(dat(1,2),npts,dmin,dmax,dmean)
c          
c     outfil='o.eqr'
c	   call sacio(outfil,dat(1,2),npts,dt,-1)
c
      do 99 j=1,nft/2
      i2=2*j
      i1=i2-1
      recv(i1)=real(dat(j,2))
      recv(i2)=aimag(dat(j,2))
99    continue
      return 
      end
c
      subroutine zer(x,start,end)
csed2      dimension x(1)
      integer start,end
      dimension x(end-start+1)
      do 1 i=start,end
   1   x(i)=0.
       return
       end

      subroutine zeerow(x,start,end)
c
c formerly called zero, name changed 10/89 to avoid
c conflicts with SAC routine names
c
csed2      dimension x(1)
      integer start,end
      dimension x(end-start+1)


      do 1 i=start,end
    1 x(i)=0.

      return
      end

      subroutine coord(x,theta,delta,y,trans,same)
c
c  transforms a vector x in one coordinate system to a vector y
c    in another coordinate system defined by strike of theta
c    and a dip of delta where y' is the strike direction and
c    z' is the dip direction of the new system wrt to x and z of
c    the old system respectively
c
c    trans defines the direction of the transform --
c     if trans = 'local' then y will be in the primed system
c     if trans = 'globe' then y will be in the original system
c
c    if same = .true. then the transformation matrix is not recalculted
c                          from the previous call

c
csed      dimension x(1),y(1),a(3,3)
      dimension x(3),y(3),a(3,3)

      character trans*5
      logical same
      integer ounit
      common /innout/ inunit,ounit
      common /cord/ a

      nfix=1
      if(same) go to 4
      cost=cos(theta)
      sint=sin(theta)
      cosd=cos(delta)
      sind=sin(delta)
      a(1,1)=cost
      a(2,1)=-cosd*sint
      a(3,1)=sind*sint
      a(1,2)=sint
      a(2,2)=cosd*cost
      a(3,2)=-sind*cost
      a(1,3)=0.
      a(2,3)=sind
      a(3,3)=cosd
   4  if(trans.eq.'globe') go to 1
      if(trans.ne.'local') go to 5
      do 2 i=1,3
    2    y(i)=a(i,1)*x(nfix)+a(i,2)*x(nfix+1)+a(i,3)*x(nfix+2)
      return
    1 do 3 i=1,3
    3    y(i)=a(1,i)*x(nfix)+a(2,i)*x(nfix+1)+a(3,i)*x(nfix+2)
      return
    5 write(ounit,101) trans
  101 format(' trans = ',a5,' in coord, no transformation done')
      return
      end

      function npowr2(n)
c
c finds the next power of 2 .ge.n
c
      ipowr=alog10(2.*float(n)-1.)/.301029996
      if(n.eq.1) ipowr=1
      npowr2=2**ipowr
      return
      end

      subroutine dfftr (x,nft,dirctn,delta)
c                                              a.shakal, 1/78, 15 jul 80
c           this subroutine does a fast fourier transform on a real
c        time series.  it requires 1/2 the storage and e1/2 the time
c        required by a complex fft.
c
c     forward transform, "call dfftr(x,nft,'forward',dt)":
c           input = x(1),x(2),..,x(nft) = real time series of nft points
c          output = x(1),x(2),..,x(nft+2) = nft/2+1 complex spectral poi
c        these spectral points are identical to the first nft/2+1 return
c        by subroutine myfft (i.e., pos freq terms).  thus, the coefficien
c        at fj, the j-th frequency point (where fj = (j-1)*delf, j=1,nft
c        and delf = 1/(nft*dt)), is in x(i-1),x(i), where i=2j.  x(1) is
c        dc term, x(2) = 0 (because real time series), x(nft+1) is real
c        of nyquist coef, and x(nft+2) is imaginary part (0 because real
c        series).
c
c     inverse transform, "call dfftr(x,nft,'inverse',delf)":
c        input and output are interchanged.
c
c           if this subroutine is called with 'forward', and then with '
c        and delf of 1/(nft*dt), the original time series is recovered.
c        identical results (but for scaling) can be obtained by calling
c        myfft(x,nft,isign), but in fft a real time series must be stored
c        complex array with zero imaginary parts, which requires 2*nft p
c        of array x.  also, the coefs returned by the fft will differ by
c        n-scaling, since fft's leave out the dt,delf of the approximate
c        integrations.  this subroutine calls fft.
c           this subroutine is a modification of the subroutine 'fftr',
c        written by c.frasier.  the principal modifications are:
c             1) the delt,delf of the integrations are included to make
c                a discrete approximation to the fourier transform.
c             2) the storage of the spectrum (on output if forward, or i
c                if inverse) has x(2) = zero, with the nyquist component
c                x(nft+1), with x(nft+2) = 0.
c
      logical forwrd, invrse
      character dirctn*7
      complex  csign, c1, c2, c3, speci, specj
csed2      real x(1)
      real x(2*nft)
      integer ounit
      common /innout/ inunit,ounit
      pi = 3.1415927
      nfix=1
c
      call locast(dirctn,invrse,forwrd)
c
      nftby2 = nft/2
      if (.not.(forwrd)) go to 20001
c            forward transform..
      call myfft (x,nftby2,-1)
      x1 = x(nfix)
      x(nfix) = x1 + x(nfix+1)
      x(nfix+1) = x1 - x(nfix+1)
      sign = -1.
      go to 20002
20001 if (.not.(invrse)) go to 10001
c            adjust nyquist element storage for inverse transform
      x(nfix+1) = x(nft+1)
      x(nft+1) = 0.
      sign = +1.
      go to 20002
10001 stop 'dirctn bad to dfftr'
c
c           manipulate elements as appropropriate for a 1/2 length
c        complex fft, after the forward fft, or before the inverse.
20002 piovrn = pi*sign/float(nftby2)
      csign = cmplx(0.,sign)
      do 10 i = 3,nftby2,2
      j = nft-i+2
      c1 = cmplx(x(i)+x(j), x(i+1)-x(j+1))
      c2 = cmplx(x(i)-x(j), x(i+1)+x(j+1))
      w = piovrn*float(i/2)
      c3 = cmplx(cos(w),sin(w))*c2
      speci = c1 + csign*c3
      x(i) = real(speci)/2.
      x(i+1) = aimag(speci)/2.
      specj = conjg(c1) + csign*conjg(c3)
      x(j) = real(specj)/2.
      x(j+1) = aimag(specj)/2.
   10 continue
      x(nftby2+2) = -x(nftby2+2)
      if (.not.(forwrd)) go to 20004
c            include dt of integration, for forward transform...
      dt = delta
      do 9000  i = 1,nft
 9000 x(i) = x(i)*dt
c            adjust storage of the nyquist component...
      x(nft+1) = x(nfix+1)
      x(nft+2) = 0.
      x(nfix+1) = 0.
      go to 20005
20004 if (.not.(invrse)) go to 10002
      x1 = x(nfix)
      x(nfix) = (x1+x(nfix+1))/2.
      x(nfix+1) = (x1-x(nfix+1))/2.
c            do the inverse transform...
      call myfft (x,nftby2,+1)
c            in the inverse transform, include the df of the integration
c            and a factor of 2 because only doing half the integration
c            (i.e., just over the positive freqs).
      twodf = 2.*delta
      do 9002  i = 1,nft
 9002 x(i) = x(i)*twodf
10002 continue
20005 return
      end

      subroutine minmax(x,npts,min,max,mean)
csed      dimension x(1)
      dimension x(npts)
      real min,max,mean
      min=9.0e+19
      max=-9.0e+19
      mean=0.
      do 1 i=1,npts
           if(x(i).gt.max) max=x(i)
           if(x(i).lt.min) min=x(i)
           mean=mean + x(i)
    1 continue
      mean=mean/float(npts)
      return
      end

      subroutine locast(dirctn,invrse,forwrd)
      character dirctn*7
      logical forwrd,invrse
      integer ounit
      common /innout/ inunit,ounit
      if(dirctn.eq.'forward') go to 1
      if(dirctn.eq.'inverse') go to 2
      write(ounit,100)dirctn
  100 format(1x,a7,2x,'is meaningless to dfftr, use forward or inverse
     *only')
      invrse=.false.
      forwrd=.false.
      return
    1 invrse=.false.
      forwrd=.true.
      return
    2 invrse=.true.
      forwrd=.false.
      return
      end

      subroutine myfft(data,nn,isign)
c  Formerly s/r fft, name changed in 10/89 to avoid conflicts
c    with SAC routine names
c
c                                              a.shakal, 1/78, 10 jul 80
c        cooley-tukey 'fast fourier trnasform' in ansi fortran 77.
c
c           transform(j) = sum {data(i)*w**u(i-1)*(j-1)e}, where i and
c        j run from 1 to nn, and w = exp(sign*twopi*sqrtu-1e/nn).
c        data is a one-dimensional complex array (i.e., the real and
c        imaginary parts of the data are located immediately adjacent
c        in storage, such as fortran places them) whose length nn is
c        a power of two.  isign is +1 or -1, giving the sign of the
c        transform.  transform values are returned in array data,
c        replacing the input data.  the time is proportional to
c        n*log2(n), rather than the non-fft n**2.  modified from the
c        fortran ii coding from n.brenner's mit-ll tech rept.
c
csed      real data(1)
      real data(2*nn)
      pi = 3.1415926
c
      n = 2*nn
      j = 1
      do 5 i = 1,n,2
      if (.not.(i .lt. j)) go to 20001
      tempr = data(j)
      tempi = data(j+1)
      data(j) = data(i)
      data(j+1) = data(i+1)
      data(i) = tempr
      data(i+1) = tempi
20001 m = n/2
    3 if (.not.(j .gt. m)) go to 20004
      j = j-m
      m = m/2
      if (m .ge. 2) go to 3
20004 j = j+m
   5  continue
c
c
      mmax = 2
    6 if (.not.(mmax .ge. n)) go to 20007
      return
20007 if (.not.(mmax .lt. n)) go to 10001
      istep = 2*mmax
      pibymx = pi*float(isign)/float(mmax)
c
      do 8 m = 1,mmax,2
      theta = pibymx*float(m-1)
      wr = cos(theta)
      wi = sin(theta)
      do 8 i = m,n,istep
      j = i + mmax
      tempr = wr*data(j) - wi*data(j+1)
      tempi = wr*data(j+1) + wi*data(j)
      data(j) = data(i) - tempr
      data(j+1) = data(i+1) - tempi
      data(i) = data(i) + tempr
      data(i+1) = data(i+1) + tempi
   8  continue
      mmax = istep
      go to 6
10001 continue
20008 return
      end
