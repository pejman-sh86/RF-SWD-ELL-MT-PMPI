           
;================================================================================             
          pro ray_jan,nh,ns,xs,xh,z,c,csq,nz,t,iA,A            
;================================================================================          

;  izs ... index of z2 for source depth
;  izh ... index of z2 for h/ph depth

t=dblarr(nh*ns)
it=0
for ih=0,nh-1 do begin 
for is=0,ns-1 do begin            


;--------------------------------------------------------------------------------
;  Adjust depths of the ssp to correspond to source and hydrophone depths. 
;  Form vectors dz,cc,ccsq and ccdz once initially for efficiency.
;-------------------------------------------------------------------------------- 

r = sqrt((xs(is*3)-xh(ih*3))^2+(xs(is*3+1)-xh(ih*3+1))^2)
       
dz = z(izs+1:izh)-z(izs:izh-1)
cc = c(izs:izh-1)
ccsq = csq(izs:izh-1)
ccsq = cc^2
ccdz = cc*dz  
ch = (z2(izh)-z2(izs))/total(dz/cc) 
   
p0 = r/sqrt(r^2+(xs(is*3+2)-xh(ih*3+2))^2)/ch
rp0 = p0*total(ccdz/sqrt(1.-p0^2*ccsq))

;------------------------------------------------------------------------------- 
;  Newton's method to refine the ray-parameter value.
;------------------------------------------------------------------------------- 

for i=0,20 do begin       
   drdp=total(ccdz/sqrt(1.-p0^2*ccsq)+p0^2*ccsq*ccdz/(1.-p0^2*ccsq)^1.5)
   p1=p0+(r-rp0)/drdp 
   rp1=p1*total(ccdz/sqrt(1.-p1^2*ccsq))
   if (abs(rp1-r) lt 0.01) then goto, OUT       
   p0=p1
   rp0=rp1            
endfor    
print,'Ray_time has not converged'  
stop 
     

;------------------------------------------------------------------------------- 
;  Once a good ray-parameter has been found, compute the travel time t
;  and partial derivatives of t.
;------------------------------------------------------------------------------- 
        
OUT:
p=p1 
t(it) = total(dz/(cc*sqrt(1.-p^2*ccsq)))
it=it+1

     
endfor
endfor
             
return             
end            
