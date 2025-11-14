%==============================================================================
 
function [t,ang,theta,status,p,izs,izr,c,z] ...
= rayfast_shru(paths,z_in,c_in,zs_in,zr_in,r,r_tol); 

%==============================================================================

%------------------------------------------------------------------------------
%  RAYFAST  --- Stan Dosso, SEOS/UVic March, 2005
%------------------------------------------------------------------------------
%
%  Subroutine to compute travel-times and take-off angles for 
%  direct and reflected rays. 
%
%  Newton's method is applied to determine the ray-parameter 
%  p = cos(theta(z))/c(z) based on an initial estimate of
%  a straight ray through a uniform SSP set to the harmonic 
%  mean, and an analytic (integral) expression for the partial 
%  derivative drdp.  
%
%  Surface and/or bottom reflections are handled using the Method
%  of Images, ie, by reflecting source and/or receiver and SSP about 
%  the boundaries. 
%
%  This provides a very fast method to determine eigenrays to high 
%  precision, but is not applicable when ray curvature is large, 
%  eg, for turning rays.
%
%------------------------------------------------------------------------------
%  
%  INPUT PARAMETERS:
%
%     paths ... 6-element array specifying ray paths. 
%               Set elements equal to 1 to trace the ray, 0 to omit. 
%               paths(1) ... direct arrival
%               paths(2) ... surface reflection
%               paths(3) ... bottom reflection
%               paths(4) ... surface-bottom reflection
%               paths(5) ... bottom-surface reflection
%               paths(6) ... surface-bottom-surface reflection
%
%     z_in  ... Input depth values of SSP
%
%     c_in  ... Input sound-speed values of SSP
%
%     zs_in ... Input source depth
%
%     zr_in ... Input receiver depth
%
%     r     ... Source-reciever range
%
%     r_tol ... Required precision for matching r
%
%------------------------------------------------------------------------------
%  
%  OUTPUT PARAMETERS:
%
%     t     ... 6-element array of computed travel times.
%
%     ang   ... 6-element array of computed take-off angles 
%               (in degrees, relative to the horizontal).
%
%     status ...6-element flag array indicating success:
%               status(ipath) = 0 ... ray-tracing successful 
%               status(ipath) = 1 ... ray-tracing failed 
%
%------------------------------------------------------------------------------
%format('long','g');

%------------------------------------------------------------------------------
%  Initialize parameters.
%------------------------------------------------------------------------------  

c  = c_in;
z  = z_in;
zs = zs_in;
zr = zr_in;
len = 0;
len2 = 0;
len3 = 0;
rr = 0;
theta = 0;

nz = length(z);

t      = zeros(1,6,'double');
p      = zeros(1,6,'double');
ang    = zeros(1,6,'double');
status = zeros(1,6,'int8');

%------------------------------------------------------------------------------
%  Insert source depth into sound-speed profile.
%------------------------------------------------------------------------------   
[dzs,izs] = min(abs(z-zs)); 
if((zs < z(izs)) & (zs > 0)) izs = izs-1;end
if((dzs ~= 0.) & (izs ~= nz))
   czs = c(izs)+(zs-z(izs))*(c(izs+1)-c(izs))/(z(izs+1)-z(izs));
   c   = [c(1:izs) czs c(izs+1:nz)];
   z   = [z(1:izs)  zs z(izs+1:nz)];
   nz  = nz+1;
   izs = izs+1;
end

%------------------------------------------------------------------------------
%  Insert receiver depth into sound-speed profile.
%------------------------------------------------------------------------------               
                     
[dzr,izr] = min(abs(z-zr));
if((zr < zs) & (zr ~= z(izr))) izs = izs+1;end
if((zr < z(izr))) & ((zr > 0)) izr = izr-1;end
if(dzr ~= 0.)
   czr = c(izr)+(zr-z(izr))*(c(izr+1)-c(izr))/(z(izr+1)-z(izr));
   c   = [c(1:izr),czr,c(izr+1:nz)];
   z   = [z(1:izr), zr,z(izr+1:nz)];
   nz  = nz+1;
   izr = izr+1;
end

%------------------------------------------------------------------------------
%  Begin loop over requested raypaths. 
%------------------------------------------------------------------------------

for ipath=1:6

  if(paths(ipath) == 1)
  
  
%------------------------------------------------------------------------------  
%  Define the "image" SSP (nzi, zi, ci) and source & receiver indices
%  (izsi, izri) reflected about sea-surface and/or sea-floor for each
%  path requested.
%------------------------------------------------------------------------------
  
%------------------------------------------------------------------------------
%  Direct path 
%------------------------------------------------------------------------------ 

      if((ipath == 1) & (paths(ipath) == 1))
         if(zs == zr)
            fprintf(1,'\n Horizontal direct ray not permitted\n\n');
            status(ipath) = 1;
            break;
         end
         ci   = c;
         zi   = z;
         nzi  = nz;
         izsi = izs;
         izri = izr;
      end
      
      
%------------------------------------------------------------------------------
%  Surface reflection 
%------------------------------------------------------------------------------ 

      if((ipath == 2) & (paths(ipath) == 1))
         ci   = [ fliplr(c) c(2:nz)];
         zi   = [-fliplr(z) z(2:nz)];
         nzi  = length(zi);
         izsi = (nz+1)-izs;
         izri = (nzi)-(nz-izr);
      end


%------------------------------------------------------------------------------
%  Bottom reflection 
%------------------------------------------------------------------------------

      if((ipath == 3) & (paths(ipath) == 1))
         ci   = [c fliplr(c(1:nz-1))];
         zi   = [z 2*z(nz)-fliplr(z(1:nz-1))];
         nzi  = length(zi);
         izsi = izs;
         izri = (nzi+1)-izr;
      end


%------------------------------------------------------------------------------
%  Surface-Bottom reflection 
%------------------------------------------------------------------------------

      if((ipath == 4) & (paths(ipath) == 1))
         ci   = [ fliplr(c) c(2:nz) fliplr(c(1:nz-1))];
         zi   = [-fliplr(z) z(2:nz) 2*z(nz)-fliplr(z(1:nz-1))];
         nzi  = length(zi);
         izsi = (nz+1)-izs;
         izri = (nzi+1)-izr;
      end


%------------------------------------------------------------------------------
%  Bottom-Surface reflection 
%------------------------------------------------------------------------------

      if((ipath == 5) & (paths(ipath) == 1))
         ci   = [c fliplr(c(1:nz-1)) c(2:nz)];
         zi   = [z 2*z(nz)-fliplr(z(1:nz-1)) 2*z(nz)+z(2:nz)];
         nzi  = length(zi);
         izsi = izs;
         izri = (nzi)-(nz-izr);
      end


%------------------------------------------------------------------------------
%  Surface-Bottom-Surface reflection 
%------------------------------------------------------------------------------

      if((ipath == 6) & (paths(ipath) == 1))
         ci   = [fliplr(c) c(2:nz-1) fliplr(c),c(2:nz)];
         zi   = [-fliplr(z) z(2:nz-1) 2*z(nz)- ...
                  fliplr(z) 2*z(nz)+z(2:nz)];
         nzi  = length(zi);
         izsi = (nz+1)-izs;
         izri = (nzi)-(nz-izr);
      end
          

%------------------------------------------------------------------------------
%  Define "image" source & receiver depths (zsi, zri). Form SSP gradient (gi).
%  Store SSP quantities for vector ray integrations. Compute harmonic mean 
%  sound speed (ch).
%------------------------------------------------------------------------------
     
      zsi = zi(izsi);
      zri = zi(izri);
         
      gi = (ci(2:nzi)-ci(1:nzi-1))./(zi(2:nzi)-zi(1:nzi-1));
      gi = [gi gi(end)];
     
      zstart  = min([zsi zri]);
      zend    = max([zsi zri]);
      izstart = min([izsi izri]);
      izend   = max([izsi izri]);
      if(ipath == 3)
        izr;
        izs;
        iizri = izr - izs;
      end
      if(ipath == 4)
        iizri = izr - izs + 2 * (izs - 1);
      end
      c1      = ci(izstart:izend-1);
      c2      = ci(izstart+1:izend);
      g1      = gi(izstart:izend-1);
      z1      = zi(izstart:izend-1);
      z2      = zi(izstart+1:izend);
      dz      = z2-z1;

%      This is exatly c2, no ???
%      abs(g1.*z2+c1-g1.*z1)


      ch = abs(zi(izend)-zi(izstart))/ ...
           sum((log(abs(g1.*z2+c1-g1.*z1)./c1))./g1);
           
%------------------------------------------------------------------------------
%  Check for special case of vertical ray.
%------------------------------------------------------------------------------

      if(r == 0.)
         t(ipath)   = abs(zsi-zri)/ch;
         p(ipath)   = 0.;
         ang(ipath) = 90.;
         break;
      end
      

%------------------------------------------------------------------------------
%  Compute initial ray-parameter estimate (p0) from straight-line path
%  with constant sound speed equal to harmonic mean (ch). Call subroutine
%  rayint_direct to compute the ray range (r0), time (t0) and partial 
%  derivative (drdp) for p0.
%------------------------------------------------------------------------------

      p0 = r/(sqrt(r^2+(zri-zsi)^2)*ch);
      
      [rr0,r0,t0,drdp,s0,status(ipath)] = rayint_direct(p0,c1,c2,g1,dz);
      
      if(status(ipath) == 1)
         fprintf(1,'\n Rayfast failed to find eigenray. ');
         fprintf(1,'zs, zr, r =%10.4f\t%10.4f\t%10.4f\n',zs,zr,r);
         return;
      end
      
      
%------------------------------------------------------------------------------
%  Apply Newton's method to refine rayparamter (p1) until ray 
%  range (r1) agrees with desired range (r) to required accuracy
%  (r_tol).
%------------------------------------------------------------------------------
   
      conv = 0;
      for iter=1:20
           
         p1 = p0+(r-r0)/drdp;
               
         [rr1,r1,t1,drdp,s,status(ipath)] = rayint_direct(p1,c1,c2,g1,dz);
         
         if(status(ipath) == 1)
            fprintf(1,'\n Rayfast failed to find eigenray. ');
            fprintf(1,'zs, zr, r =%10.4f\t%10.4f\t%10.4f\n',zs,zr,r);
            return
         end
         
         if(ipath == 3)
             iizri;
             rr = sum(rr1(1:iizri));
             len2 = sum(s(1:iizri));
             theta = (acos(p1*c(nz))*180/pi);
         end
         if(ipath == 4)
             iizri;
             len3 = sum(s(1:iizri));
         end

         if(abs(r-r1) < r_tol) 
           conv = 1;
           break;
         end
         
         p0 = p1;
         r0 = r1;
 
      end
      
      if(conv == 0)
        fprintf(1,'\n Rayfast has not converged. zs, zr =');
        fprintf(1,'zs, zr, r =%10.4f\t%10.4f\t%10.4f\n',zs,zr,r);
        status(ipath) = 1
%      else
%        fprintf(1,'\n\n Converged!\n');
      end
  

%------------------------------------------------------------------------------
%  Store computed ray time (t) and ray take-off angle (ang). 
%------------------------------------------------------------------------------

      p(ipath)  = p1;
      t(ipath)   = t1;
      ang(ipath) = acos(p1*ci(izsi))*180.0/pi;
      len(ipath) = sum(s);

   end
end

return;


%===============================================================================        
 
function [rr,r,t,drdp,s,status] = rayint_direct(p,c1,c2,g1,dz);

%===============================================================================

%------------------------------------------------------------------------------
%  Integration of range r, travel-time t, and partial drdp for a 
%  non-turning ray given ray-parameter p. For a turning ray, set 
%  flag status=1.
%------------------------------------------------------------------------------

pc1   = p*c1;
pc2   = p*c2;
pc1sq = pc1.^2;;
pc2sq = pc2.^2;
sq1   = sqrt(1.0-pc1sq);
sq2   = sqrt(1.0-pc2sq);

r    = sum((sq1-sq2)./(p*g1));
rr   = (sq1-sq2)./(p*g1);
t    = sum((log(pc2.*(1.0+sq1)./(pc1.*(1.0+sq2))))./g1);
drdp = -r/p-sum((c1.^2./sq1-c2.^2./sq2)./g1);


theta1 = acos(p.*c1);
theta2 = acos(p.*c2);
rho = 1./(p*g1);
%dis = sqrt(dz.^2 + ((sq1-sq2)./(p*g1)).^2);

%beta = acos((rho-(cos((pi/2)-theta1).*dis))./rho)
beta = abs(theta2-theta1);
s = abs(rho.*beta);


status = 0;
if(isfinite(r) == 0) 
  status = 1;
end

return;

