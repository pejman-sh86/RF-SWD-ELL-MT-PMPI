******************************************************************************
c
      subroutine itd(npts, dt, radial, vertical, maxbumps, theshift,
     +               tol, gwidth, idum, ioutput, MAXPTS, MAXG,
     +               decon, nshifts, fit, ierror) 
c
******************************************************************************
c
c     Chuck Ammon - Saint Louis University
c     November 1996 / September 1998 
c     VERSION 1.04
c
c     Based on the Kikuchi and Kanamori (1981) iterative
c      deconvolution algorithm. The output of this code
c      is the deconvolution of the "denominator" from
c      the "numerator"
c
c     The final deconvolution is called "decon.out"
c       the observed time series (original convolved with Gaussian)
c        is called "observed"
c       the predicted is call "predicted"
c
c     Header values from the "numerator" are copied into the
c       decon.out file. The gwidth is stored in user0.
c
c******************************************************************************
c
c     if you choose verbose output
c       the deconvolutions are called d???
c       the predictions are called p???
c       the residuals are called r???
c
c      where ??? is corresponds to the number of "bumps"
c
******************************************************************************
c
c     Because the correlation is computed in the 
c     Frequency Domain, the maximum lag computed is npts / 2
c
c     That is, negative lags are not used - so make sure you wavelet starts
c       before the signal.
c
c     A low-pass Gaussian filter is applied to the signals before the
c      deconvolution is performed. The predicted will not match the
c      original signal, but the filtered signal, which is stored in the
c      file called "numerator".
c
******************************************************************************
c
c     numerator: numerator file
c     denominator: denominator file
c     maxbumps:  max number of iterations
c     theshift: the phase shift (secs) for the output
c     tol: minimum percent error increase to accept
c     gwidth: Gaussian filter width factor
c     idum: Allow negative pulses? (1->y, 0->no)
c     ioutput: Minimal (0) or verbose output(1)
c     MAXPTS: maximum number of data points
c     MAXG
c
******************************************************************************
c      character*256 numerator, denominator
      integer npts, idum, ioutput, maxbumps
      real radial(npts), vertical(npts)
      real dt, theshift, tol, gwidth
      integer MAXPTS, MAXG


      real decon(MAXPTS)
      real fit
      integer nshifts, ierror

      real p(MAXPTS), r(MAXPTS)
      real amps(MAXG)
      integer shifts(MAXG)
      real f(MAXPTS), g(MAXPTS), f1(MAXPTS), g1(MAXPTS)


      character*12 resfile, filename
      integer stdin,stdout,ounit,inunit,forward,inverse
      logical lpositive,verbose
      
      stdin = 5
      stdout = 6
      ounit = 9
      inunit = 10
      forward =  1
      inverse = -1
c      gwidth = 2.5
c
c      write(stdout,*) ' '
c      write(stdout,*) 'Program iterdeconfd - Version 1.04, 1997-98'
c      write(stdout,*) 'Chuck Ammon, Saint Louis University'
c      write(stdout,*) ' '
c
      if(maxbumps .gt. MAXG)then
c       write(stdout,*) 'Maximum Number of bumps is ',MAXG
        maxbumps = MAXG
      end if
      if(idum .eq. 1) then
        lpositive = .false.
      else
        lpositive = .true.
      end if 
c
      if(ioutput .eq. 0) then
        verbose = .false.
      else
        verbose = .true.
      end if 
c
      ierror = 0
c
******************************************************************************
c     
c      call rsac1(numerator,f,npts,beg,delta,MAXPTS,nerr)
c      if(nerr .ne. 0)then
c        write(stdout,*)'Problem reading the numerator file'
c        stop
c         return
c      end if
c
c      call rsac1(denominator,g,nptsd,b,dt,MAXPTS,nerr)
c      if(nerr .ne. 0)then
c        write(stdout,*)'Problem reading the denominator file'
c        stop
c         return
c      end if      
      call zero(f,MAXPTS)
      call zero(g,MAXPTS)
      call copy_array(npts, radial, f)
      call copy_array(npts, vertical, g)
c
******************************************************************************
c
c     Find the next power of two greater than the data
c       dimensions - use the numerator, zero pad
c
      n = 1
119   continue
      if(n .ge. npts) go to 120
      n = n * 2
      go to 119
c      
120   continue
      if(n .gt. MAXPTS)then
c        write(stdout,*) 'Too many points needed.'
c        write(stdout,*) 'n = ', n
         ierror = 1
c        stop
         return
      end if
c
******************************************************************************
c     zero-pad the data
c
      n2 = npts
      npts = n
c
******************************************************************************
c     FINISHED READING FILES
c      
c     Now begin the cross-correlation procedure
c
c      Put the filter in the signals
c
c     Be careful about wsac1 and rsac1 subroutines
c     The wsac1 subroutine in the following stores the original unfiltered arrays f and g with n2 (original npts) non-zero elements 
      call copy_array(n2, f, f1)
      call copy_array(n2, g, g1)
      call gfilter(f,gwidth,npts,dt)
      call gfilter(g,gwidth,npts,dt)
c      call wsac1('numerator',f,npts,beg,dt,nerr)
c      call wsac1('observed',f,npts,beg,dt,nerr)
c      call wsac1('denominator',g,npts,beg,dt,nerr)
c      call copy_array(n2, f, f1)
c      call copy_array(n2, g, g1)
c
c     compute the power in the "numerator" for error scaling
c
      power = 0
      do 100 i = 1, npts
        power = power + f(i)*f(i)
100   continue
c
c     correlate the signals
c 
      call fcorrelate(f,g,npts,MAXPTS,dt)
c     call wsac1('ccor0',g,npts,beg,dt,nerr)
c
c     find the peak in the correlation
c
      maxlag = npts/2
c      write(stdout,'(/,a27,f10.5)') 'The maximum spike delay is ', 
c     &   real(maxlag) * dt
c
      if(lpositive) then
        call getmax(g,maxlag,amps(1),shifts(1))
      else
        call getabsmax(g,maxlag,amps(1),shifts(1))
      end if
      amps(1) = amps(1) / dt
c
      nshifts = 1
c
c     read in the signals again
c
      call zero(f,MAXPTS)
      call zero(g,MAXPTS)
c      call rsac1(numerator,f,ndummy,beg,delta,MAXPTS,nerr)
c      call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
c      ndummy equals to the n2 (original npts)
      call copy_array(n2, f1, f)
      call copy_array(n2, g1, g)
c
c     compute the predicted deconvolution result
c
      call zero(p,MAXPTS)
      call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
c      if(verbose) then
c          call phs_shift(p,theshift,npts,dt)      
c        call wsac1('d001',p,npts,-theshift,dt,nerr)
c          call phs_shift(p,-theshift,npts,dt)      
c      end if
c
c     convolve the prediction with the denominator signal
c      
      call convolve(p,g,npts,dt)
c
c      if(verbose) then
c        call wsac1('p001',p,npts,beg,dt,nerr)
c      end if
c
c     filter the signals
c     
      call gfilter(f,gwidth,npts,dt)
      call gfilter(g,gwidth,npts,dt)
c
c      if(verbose)then
c       write(resfile,'(a1,i3.3)') 'r',0
c        call wsac1(resfile,f,npts,beg,dt,nerr)
c      end if
c      
c     compute the residual (initial error is 1.0)
c
      call getres(f,p,npts,r,sumsq_ip1)
c
      sumsq_i = 1.0
      sumsq_ip1 = sumsq_ip1 / power
      d_error = 100*(sumsq_i - sumsq_ip1) 
c
c      write(resfile,'(a1,i3.3)') 'r',1
c      if(verbose)then
c        call wsac1(resfile,r,npts,beg,dt,nerr)
c      end if
c     
c      write(stdout,1000)
c      write(stdout,1001)
c     &  resfile, dt*amps(1),(shifts(1)-1)*dt,100*sumsq_ip1,
c     &  d_error
c1000  format(/,1x,'File',9x,
c     & 'Spike amplitude   Spike delay   Misfit   Improvement')
c1001  format(1x,a10,2x,e16.9,2x,f10.3,3x,f7.2,'%',3x,f9.4,'%')
c
******************************************************************************
c    
      do while(d_error .gt. tol .and. nshifts .lt. (maxbumps))
c
        nshifts = nshifts + 1
        sumsq_i = sumsq_ip1
c
        call zero(g,MAXPTS)
c        call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
        call copy_array(n2, g1, g)
        call gfilter(g,gwidth,npts,dt)
        call fcorrelate(r,g,npts,MAXPTS,dt)
        if(lpositive)then
          call getmax(g,maxlag,amps(nshifts),shifts(nshifts))
        else
         call getabsmax(g,maxlag,amps(nshifts),shifts(nshifts))
        end if
        amps(nshifts) = amps(nshifts) / dt
c
        call zero(p,MAXPTS)
        call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
c        if(verbose)then
c          write(filename,'(a1,i3.3)') 'd',nshifts
c          call phs_shift(p,theshift,npts,dt)      
c          call wsac1(filename,p,npts,-theshift,dt,nerr)
c          call phs_shift(p,-theshift,npts,dt)      
c        end if
c        
        call zero(g,MAXPTS)
c        call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
        call copy_array(n2, g1, g)
        call convolve(p,g,npts,dt)
c        if(verbose)then
c          write(filename,'(a1,i3.3)') 'p',nshifts
c          call wsac1(filename,p,npts,beg,dt,nerr)
c        end if
c                
        call zero(f,MAXPTS)
c        call rsac1(numerator,f,ndummy,beg,delta,MAXPTS,nerr)
        call copy_array(n2, f1, f)
        call gfilter(f,gwidth,npts,dt)
        call getres(f,p,npts,r,sumsq_ip1)
        
        sumsq_ip1 = sumsq_ip1/ power
c        write(resfile,'(a1,i3.3)') 'r',nshifts
c        if(verbose)then
c          call wsac1(resfile,r,npts,beg,dt,nerr)
c        end if
        d_error = 100*(sumsq_i - sumsq_ip1)

c        write(stdout,1001)
c     &   resfile,dt*amps(nshifts),(shifts(nshifts)-1)*dt,
c     &   100*sumsq_ip1,d_error
c    
      enddo
c
******************************************************************************
c      
c      write(stdout,1010) d_error
c1010  format(/,1x,'Last Error Change = ',f9.4,'%',/)
c
c     if the last change made no difference, drop it
c      
      fit = 100 - 100*sumsq_ip1
c
      if(d_error .le. tol)then
         nshifts = nshifts - 1
         fit = 100 - 100*sumsq_i
c         write(stdout,*)'Hit the min improvement tolerance - halting.'
      end if
c
c      if(nbumps .ge. maxbumps)then
c         write(stdout,*)'Hit the max number of bumps - halting.'
c      end if
c
c      write(stdout,*)'Number of bumps in final result: ', nshifts
c      write(stdout,1011) fit
c1011  format(1x,'The final deconvolution reproduces ',
c     &    f6.1,'% of the signal.',/)
c
******************************************************************************
c
c     compute the final prediction
c
c      call zero(p,MAXPTS)
c      call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
c      call zero(g,MAXPTS)
c      call rsac1(denominator,g,ndummy,b,dt,MAXPTS,nerr)
c      call convolve(p,g,npts,dt)
c      call wsac1('predicted',p,npts,beg,dt,nerr)
c      call zero(g,MAXPTS)
c
c     write out the answer
c
c      call zero(p,MAXPTS)
c      call build_decon(amps,shifts,nshifts,p,npts,gwidth,dt)
c      call phs_shift(p,theshift,npts,dt)      
      call zero(decon,MAXPTS)
      call build_decon(amps,shifts,nshifts,decon,npts,gwidth,dt)
      call phs_shift(decon,theshift,npts,dt)      
c
c      call newhdr
c      call rsac1(numerator,g,ndummy,b,dt,MAXPTS,nerr)
c      call setnhv('NPTS',npts,nerr)
c      call setfhv('B',-theshift,nerr)
c      theend = -thshift + (npts-1)*dt
c      call setfhv('E',theend,nerr)
c      call setnhv('NZSEC',-12345,nerr)
c      call setfhv('USER0',gwidth,nerr)
c     call setkhv('KUSER0','Rftn',nerr)
c     call setkhv('KUSER1','IT_DECON',nerr)
c      call wsac0('decon.out',xdummy,p,nerr)
c
c     write out the gaussian filter
c
c      if(verbose)then
c        call newhdr
c        call zero(p,MAXPTS)
c        p(1) = 1 / dt
c        call phs_shift(p,theshift,npts,dt)      
c        call gfilter(p,gwidth,npts,dt)
c        call wsac1('thefilter',p,npts,beg,dt,nerr)
c      end if
      
c      stop
      end
c
********************************************************************************

      subroutine copy_array(n, x, x1)
          IMPLICIT NONE
          integer n
          real x(n), x1(n)

          integer I

          DO I = 1, n
              x1(I) = x(I)
          END Do

      end subroutine copy_array

   
