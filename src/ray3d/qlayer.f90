subroutine qlayer( &
    lc, ang, n, h2, valpha, vbeta, rho2, qa, qb, rs, &
    w, up, wp, usv, wsv, vsh )
!**********************************************                         
!     crustal response for p and s waves      *                         
!     in layered dissipative medium           *                         
!     codedby T. Kurita & T. Mikumo           *                         
!     revised by T. Shibutani                 *                         
!**********************************************                         
!                                                                       
!
!include 'rfi_param.inc'
USE ray3d_com

parameter     ( nlmx = maxlayd )
parameter     ( nb  = maxsamp )
parameter     ( nb2 = nb/2+1 )
!
!
real w(nb2), h2(nlmx), qa(nlmx), qb(nlmx), rs(nlmx), rho2(nlmx)
!
complex alpha2(nlmx), ralpha(nlmx), valpha(nlmx),&
        beta2(nlmx), rbeta(nlmx), vbeta(nlmx),&
        gamma(nlmx), vn, c, imag, &
        en(4,4), coef1, coef2, coef3, coef4,&
        dca, dea, gca, gea,&
        up(nb2), wp(nb2), usv(nb2), wsv(nb2),&
        bm, bn, rvs, vsh(nb2)
!                                                                       
complex*16 a(4,4,nlmx), b(4,4,nlmx), aj(4,4),&
           s, da, dc, de, gc, ge, &
           al(2,2,nlmx), bl(2,2,nlmx), sl, ab, &
           p(nlmx), q(nlmx), cp, sp2, cq, sq
!
parameter ( imag=(0.0,1.0) )                                      
!                                                                       
!     *** complex velocity for dissipative medium ***
!
do 1 m=1,n                                                        
  alpha2(m)=valpha(m)+imag*valpha(m)/(2.0*qa(m))&
          +valpha(m)/(8.0*qa(m)**2)
  beta2(m)=vbeta(m)+imag*vbeta(m)/(2.0*qb(m))& 
         +vbeta(m)/(8.0*qb(m)**2)
1 continue
!
      if (lc.eq.1) then
         vn=alpha2(n)
      else
         vn=beta2(n)
      end if
!
!     *** apparent velocity ***
!
      c=vn/sin(ang)
!
      do 2 m=1,n
!
        cal=real((c/alpha2(m))**2-1. )
        if (cal.ge.0.) then
           ralpha(m)=csqrt((c/alpha2(m))**2-1. )
        else
           ralpha(m)=-imag*csqrt(1. -(c/alpha2(m))**2)
        end if
!
        cbe=real((c/beta2(m))**2-1. )                                    
        if (cbe.ge.0.) then
           rbeta(m)=csqrt((c/beta2(m))**2-1. )
        else
           rbeta(m)=-imag*csqrt(1. -(c/beta2(m))**2)
        end if
!
        gamma(m)=2.*(beta2(m)/c)**2
!
    2 continue
!
      if (lc.le.2) then
!                                                                       
!     *** inverse matrix at the lowermost interface for P-SV ***
!
         en(1,1)=-2.*(beta2(n)/alpha2(n))**2
         en(1,2)=0.
         en(1,3)=1./(rho2(n)*alpha2(n)**2)
         en(1,4)=0.
         en(2,1)=0.
         en(2,2)=c**2*(gamma(n)-1.)/(alpha2(n)**2*ralpha(n))
         en(2,3)=0.
         en(2,4)=1./(rho2(n)*alpha2(n)**2*ralpha(n))
         en(3,1)=(gamma(n)-1.)/(gamma(n)*rbeta(n))
         en(3,2)=0.
         en(3,3)=-1./(rho2(n)*c**2*gamma(n)*rbeta(n))
         en(3,4)=0.
         en(4,1)=0.
         en(4,2)=1.
         en(4,3)=0.
         en(4,4)=1./(rho2(n)*c**2*gamma(n))
         coef1=2.*c**2/alpha2(n)**2
         coef2=coef1/ralpha(n)
         coef3=2.0/gamma(n)
         coef4=coef3/rbeta(n)
!
      end if
!
!     *** loop on frequency ***
!
      do 5 kk=1,nb2                                                     
!
        nn=n-1                                                          
        do 3 m=1,nn                                                     
          p(m)=w(kk)*ralpha(m)*h2(m)*rs(m)/c                             
          q(m)=w(kk)*rbeta(m)*h2(m)*rs(m)/c                              
          cp=cdcos(p(m))                                                
          sp2=cdsin(p(m))                                                
          cq=cdcos(q(m))                                                
          sq=cdsin(q(m))                                                
          bm=beta2(m)**2*rho2(m)*rbeta(m)/rs(m)
!                           
!         **  matrix at layer interfaces  **                             
!
    if (lc.le.2) then
!
!            ** for P-SV problem
!
       a(1,1,m)=gamma(m)*cp-(gamma(m)-1.0)*cq
       a(1,2,m)=imag*((gamma(m)-1.0)*sp2/ralpha(m)&
            +gamma(m)*rbeta(m)*sq)
       a(1,3,m)=-(cp-cq)/(rho2(m)*c**2)
       a(1,4,m)=imag*(sp2/ralpha(m)+rbeta(m)*sq)/(rho2(m)*c**2)
       a(2,1,m)=-imag*(gamma(m)*ralpha(m)*sp2&
            +(gamma(m)-1.0)*sq/rbeta(m))
       a(2,2,m)=-(gamma(m)-1.0)*cp+gamma(m)*cq
       a(2,3,m)=imag*(ralpha(m)*sp2+sq/rbeta(m))/(rho2(m)*c**2)
       a(2,4,m)=a(1,3,m)
       a(3,1,m)=rho2(m)*c**2*gamma(m)*(gamma(m)-1.0)*(cp-cq)
       a(3,2,m)=imag*(rho2(m)*c**2*((gamma(m)-1.0)**2*sp2/ralpha(m)&
            +gamma(m)**2*rbeta(m)*sq))
       a(3,3,m)=a(2,2,m)
       a(3,4,m)=a(1,2,m)
       a(4,1,m)=imag*(rho2(m)*c**2*(gamma(m)**2*ralpha(m)*sp2 &
            +(gamma(m)-1.0)**2*sq/rbeta(m)))
       a(4,2,m)=a(3,1,m)
       a(4,3,m)=a(2,1,m)
       a(4,4,m)=a(1,1,m)
    else
!
!      ** for SH problem **
!
       al(1,1,m)=cq
       al(2,2,m)=al(1,1,m)
       al(1,2,m)=imag*sq/bm
       al(2,1,m)=imag*sq*bm
!
    end if
!                                          
    3   continue
!
        if (lc.le.2) then
!
!       **  matrix product for P-SV case **
!
           m=n-1
    4      continue
           if(m<=1)go to 40
           do 20 i=1,4
           do 20 k=1,4
             s=0.
             do 10 j=1,4
               s=s+a(i,j,m)*a(j,k,m-1)
   10        continue
             b(i,k,m-1)=s
   20      continue
           do 30 i=1,4
           do 30 k=1,4
             a(i,k,m-1)=b(i,k,m-1)
   30      continue
           m=m-1
           go to 4
   40      continue
           do 60 i=1,4
           do 60 k=1,4
             s=0.
             do 50 j=1,4
               s=s+en(i,j)*a(j,k,1)
50        continue
       aj(i,k)=s
60      continue
     da=(aj(1,1)-aj(2,1))*(aj(3,2)-aj(4,2))-(aj(1,2)-aj(2,2)) &
          *(aj(3,1)-aj(4,1))
  if (lc.eq.1) then
!
!       **  up horizontal/total, wp vertical/total of P  **             
!
           dc=aj(4,2)-aj(3,2)
           de=aj(4,1)-aj(3,1)
           dca=dc/da
           dea=de/da
           up(kk)=dca*coef1*alpha2(n)/c
           wp(kk)=dea*coef2*alpha2(n)*ralpha(n)/c
!                           
  else if (lc.eq.2) then
!
!     **  usv horizontal/total, wsv vertical/total of SV  **          
!
           gc=aj(1,2)-aj(2,2)
           ge=aj(2,1)-aj(1,1)
           gca=gc/da
           gea=ge/da
           usv(kk)=gca*coef4*beta2(n)*rbeta(n)/c
           wsv(kk)=gea*coef3*beta2(n)/c
!
  end if
!
else
!
!       **  matrix product for SH waves  **                             
!
         m=n-1
  114    continue
         if(m.eq.1) go to 140
         do 120 i=1,2
         do 120 k=1,2
           sl=0.0
           do 110 j=1,2
             sl=sl+al(i,j,m)*al(j,k,m-1)
  110      continue
           bl(i,k,m-1)=sl
  120    continue
         do 130 i=1,2
         do 130 k=1,2
           al(i,k,m-1)=bl(i,k,m-1)
  130    continue
         m=m-1
         go to 114
  140    continue
         bn=beta2(n)**2*rho2(n)*rbeta(n)/rs(n)
         ab=al(2,1,1)+bn*al(1,1,1)
!                                       
!       **  vsh horizontal/total of SH  **                              
!
         rvs=2.*bn/ab
         vsh(kk)=rvs
!
end if
!                                                     
5 continue
! 
return
end
