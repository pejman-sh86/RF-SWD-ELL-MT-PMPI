      SUBROUTINE CORREL(DATA1,DATA2,N,ANS,NMAX)
      INTEGER NMAX
      DIMENSION DATA1(N),DATA2(N)
      COMPLEX ANS(N)
      COMPLEX,ALLOCATABLE :: FFT(:)
      ALLOCATE(FFT(NMAX))
      CALL TWOFFT(DATA1,DATA2,FFT,ANS,N)
      NO2=FLOAT(N)/2.0
      DO 11 I=1,N/2+1
        ANS(I)=FFT(I)*CONJG(ANS(I))/NO2
11    CONTINUE
      ANS(1)=CMPLX(REAL(ANS(1)),REAL(ANS(N/2+1)))
      CALL REALFT(ANS,N/2,-1)
      RETURN
      END
      SUBROUTINE REALFT(DATA,N,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      THETA=6.28318530717959D0/2.0D0/DBLE(N)
      C1=0.5
      IF (ISIGN.EQ.1) THEN
        C2=-0.5
        CALL FOUR1(DATA,N,+1)
      ELSE
        C2=0.5
        THETA=-THETA
      ENDIF
      WPR=-2.0D0*DSIN(0.5D0*THETA)**2
      WPI=DSIN(THETA)
      WR=1.0D0+WPR
      WI=WPI
      N2P3=2*N+3
      DO 11 I=2,N/2+1
        I1=2*I-1
        I2=I1+1
        I3=N2P3-I2
        I4=I3+1
        WRS=SNGL(WR)
        WIS=SNGL(WI)
        H1R=C1*(DATA(I1)+DATA(I3))
        H1I=C1*(DATA(I2)-DATA(I4))
        H2R=-C2*(DATA(I2)+DATA(I4))
        H2I=C2*(DATA(I1)-DATA(I3))
        DATA(I1)=H1R+WRS*H2R-WIS*H2I
        DATA(I2)=H1I+WRS*H2I+WIS*H2R
        DATA(I3)=H1R-WRS*H2R+WIS*H2I
        DATA(I4)=-H1I+WRS*H2I+WIS*H2R
        WTEMP=WR
        WR=WR*WPR-WI*WPI+WR
        WI=WI*WPR+WTEMP*WPI+WI
11    CONTINUE
      IF (ISIGN.EQ.1) THEN
        H1R=DATA(1)
        DATA(1)=H1R+DATA(2)
        DATA(2)=H1R-DATA(2)
      ELSE
        H1R=DATA(1)
        DATA(1)=C1*(H1R+DATA(2))
        DATA(2)=C1*(H1R-DATA(2))
        CALL FOUR1(DATA,N,-1)
      ENDIF
      RETURN
      END
      SUBROUTINE FOUR1(DATA,NN,ISIGN)
      REAL*8 WR,WI,WPR,WPI,WTEMP,THETA
      DIMENSION DATA(*)
      N=2*NN
      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=SNGL(WR)*DATA(J)-SNGL(WI)*DATA(J+1)
            TEMPI=SNGL(WR)*DATA(J+1)+SNGL(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF
      RETURN
      END
      SUBROUTINE TWOFFT(DATA1,DATA2,FFT1,FFT2,N)
      DIMENSION DATA1(N),DATA2(N)
      COMPLEX FFT1(N),FFT2(N),H1,H2,C1,C2
      C1=CMPLX(0.5,0.0)
      C2=CMPLX(0.0,-0.5)
      DO 11 J=1,N
        FFT1(J)=CMPLX(DATA1(J),DATA2(J))
11    CONTINUE
      CALL FOUR1(FFT1,N,1)
      FFT2(1)=CMPLX(AIMAG(FFT1(1)),0.0)
      FFT1(1)=CMPLX(REAL(FFT1(1)),0.0)
      N2=N+2
      DO 12 J=2,N/2+1
        H1=C1*(FFT1(J)+CONJG(FFT1(N2-J)))
        H2=C2*(FFT1(J)-CONJG(FFT1(N2-J)))
        FFT1(J)=H1
        FFT1(N2-J)=CONJG(H1)
        FFT2(J)=H2
        FFT2(N2-J)=CONJG(H2)
12    CONTINUE
      RETURN
      END
*
******************************************************************************
******************************************************************************
*
*
*
*
************************************************************************
*
*     correlation routine - correlates f and g and replaces the 
*       g with the cross-correlation the value is normalized
*       by the zero-lag autocorrelation of g
*
************************************************************************
*
      subroutine fcorrelate(f,g,n,MAXPTS,dt)
      real f(MAXPTS), g(MAXPTS)
      real sum0, temp
      integer i,n,n2,n2o2

      real, allocatable :: c(:)
      allocate(c(MAXPTS))

c
c     compute the zero-lag autocorrelation of g
c
      sum0 = 0
      do 1 i = 1, n
        sum0 = sum0 + g(i)*g(i)
1     continue
      sum0 = sum0 * dt
c
c     compute the next power of 2 greater than n
c
      n2 = 1
5     n2 = n2 * 2
      if(n2 .lt. n) go to 5
c     
6     continue
      n2o2 = n2 / 2
c
c     Use the Numerical Recipes routine to compute the cross correlation
c
      call correl(f,g,n2,c,MAXPTS)
c 
      temp = dt / sum0
c
      do 20 i = 1,n2
        g(i) = c(i) * temp
20    continue
      
      return
      end   
*
************************************************************************
*
*     zero a real array
*
************************************************************************
*
      subroutine zero(x,n)
      real x(n)
      integer i,n
      
      do 1 i = 1,n
        x(i) = 0
1     continue
      return
      end
*
************************************************************************
*
*     get max value of array and its index
*
************************************************************************
*
      subroutine getmax(x,n,maxvalue,maxindex)
      real x(n), maxvalue
      integer i,n,maxindex
      
      maxvalue = x(1)
      maxindex = 1
      do 20 i = 2, n
        if(x(i) .gt. maxvalue) then
           maxvalue = x(i)
           maxindex = i
        end if
20    continue

      return
      end
*
************************************************************************
*
*     find max absolute value of array and its index
*
************************************************************************
*
      subroutine getabsmax(x,n,thevalue,maxindex)
      real x(n), maxvalue, thevalue
      integer i,n,maxindex
      
      maxvalue = abs(x(1))
      maxindex = 1
      thevalue = x(1)
      do 20 i = 2, n
        if(abs(x(i)) .gt. maxvalue) then
           maxvalue = abs(x(i))
           thevalue = x(i)
           maxindex = i
        end if
20    continue

      return
      end
*
*
************************************************************************
*
*     getres
*
************************************************************************
*
      subroutine getres(x,y,n,r,sumsq)
      real x(n), y(n), r(n), sumsq
      integer i,n
      
      sumsq = 0 
      do 20 i = 1, n
       r(i) = x(i) - y(i)
       sumsq = sumsq + r(i)*r(i) 
20    continue

      return
      end
*
************************************************************************
*
*     compute the predicted time series from a set of
*       amplitudes and shifts
*
************************************************************************
*
      subroutine build_decon(amps,shifts,nshifts,p,n,gwidth,dt)
      real p(n), amps(nshifts)
      integer shifts(nshifts)
      integer i, n, nshifts
      
      call zero(p,n)
      do 1 i = 1, nshifts
        p(shifts(i)) = p(shifts(i)) + amps(i)
1     continue

      call gfilter(p,gwidth,n,dt)

      return
      end
*
************************************************************************
*
*     convolve a function with a unit-area Gaussian filter.
*
************************************************************************
*
      subroutine gfilter(x,gwidth_factor,n,dt)
      real x(n), pi, two_pi, gauss, d_omega, omega
      real gwidth, gwidth_factor, sum
      integer i, j, n, n2, halfpts
      integer forward, inverse
c      
      forward = 1
      inverse = -1
      pi = acos(-1.0)
      two_pi = 2 * pi
      sum = 0
c      
      n2 = 1
1     n2 = n2 * 2
      if(n2 .ge. n) goto 2
      go to 1
c      
2     continue
      halfpts = n2 / 2
c
      call realft(x,halfpts,forward)
c
      df = 1 / (float(n2) * dt)
      d_omega = two_pi * df
      gwidth = 4.0*gwidth_factor*gwidth_factor
c 
c     Handle the nyquist frequency
c
      omega = two_pi/(2.0*dt)
      gauss = exp(-omega*omega / gwidth)
      x(2) = x(2) * gauss
c 
      do 5 i = 2, halfpts
          j = i*2
          omega = (i-1) * d_omega
          gauss = exp(-omega*omega / gwidth)
          x(j-1) = x(j-1) * gauss
          x(j)   = x(j)   * gauss
5     continue 
c
      call realft(x,halfpts,inverse)
c      
      scalefactor = dt * (2 * df)
      do 10 i = 1, n
        x(i) = x(i) * scalefactor
10    continue    
c
      return
      end


*
************************************************************************
*
*     replace x with x convolved with y
*
************************************************************************
*
      subroutine convolve(x,y,n,dt)
      real x(n), y(n), dt, scale0, scale1
      integer i, j, n, n2, halfpts
      integer forward, inverse, stderr
c      
      forward = 1
      inverse = -1
      stderr = 6
c      
      if(mod(n,2) .ne. 0) then
        write(stderr,*) 'Error in convolve - n .ne. power of two.'
        stop
      endif
      n2 = n
      halfpts = n2 / 2
c
      call realft(x,halfpts,forward)
      call realft(y,halfpts,forward)
c
      df = 1 / (float(n2) * dt)
c 
c     Handle the zero & nyquist frequency
c
      x(1) = x(1) * y(1)
      x(2) = x(2) * y(2)
c 
      do 5 i = 2, halfpts
          j = i*2
          a = x(j-1)
          b = x(j)
          c = y(j-1)
          d = y(j)
          x(j-1) = a * c - b * d
          x(j)   = a * d + b * c
5     continue  
c
      call realft(x,halfpts,inverse)
      call realft(y,halfpts,inverse)
c      
      scale0 = dt * (2 * df)
      scale1 = dt * scale0
      do 10 i = 1, n
        y(i) = y(i) * scale0
        x(i) = x(i) * scale1
10    continue    
c
      return
      end
*
************************************************************************
*
*     phase shifts a signal
*      
*
************************************************************************
*
      subroutine phs_shift(x,theshift,n,dt)
      real x(n), pi, two_pi, theshift, d_omega, omega
      integer i, j, n, n2, halfpts
      integer forward, inverse
c      
      forward = 1
      inverse = -1
      pi = acos(-1.0)
      two_pi = 2 * pi
c      
      n2 = 1
1     n2 = n2 * 2
      if(n2 .ge. n) goto 2
      go to 1
c      
2     continue
      halfpts = n2 / 2
c
      call realft(x,halfpts,forward)
c
      df = 1 / (float(n2) * dt)
      d_omega = two_pi * df
c 
c     Handle the nyquist frequency
c
      omega = two_pi/(2.0*dt)
      x(2) = x(2) * cos(omega * tshift)
c 
      do 5 i = 2, halfpts
          j = i*2
          omega = (i-1) * d_omega
          a = x(j-1)
          b = x(j)
          c = cos(omega*theshift)
          d = sin(omega*theshift)
          x(j-1) = a*c-b*d
          x(j)   = a*d+b*c
5     continue 
c
      call realft(x,halfpts,inverse)
c      
      scalefactor = dt * (2 * df)
      do 10 i = 1, n
        x(i) = x(i) * scalefactor
10    continue    
c
      return
      end

