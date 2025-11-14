%=======================================================================
 
function [t,rp1,status] = rayfast3(paths,z_in,c_in,zs_in,zr_in,r,p1,r_tol);

%=======================================================================

%-----------------------------------------------------------------------
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
%------------------------------------------------------------------------
%  
%  OUTPUT PARAMETERS:
%
%     t     ... 6-element array of computed travel times.
%
%     status ...6-element flag array indicating success:
%               status(ipath) = 0 ... ray-tracing successful 
%               status(ipath) = 1 ... ray-tracing failed 
%
%------------------------------------------------------------------------
format('long','g');

%------------------------------------------------------------------------
%  Initialize parameters.
%------------------------------------------------------------------------  

c  = c_in;
z  = z_in;
zs = zs_in;
zr = zr_in;

nz = length(z);

t      = zeros(1,6,'double');
status = zeros(1,6,'int8');
  
%-------------------------------------------------------------------------
%  Insert source depth into sound-speed profile.
%-------------------------------------------------------------------------   

[dzs,izs] = min(abs(z-zs)); 
if((zs < z(izs)) & (zs > 0)) izs = izs-1;end
if((dzs ~= 0.) & (izs ~= nz))
   czs = c(izs);
   c   = [c(1:izs) czs c(izs+1:nz)];
   z   = [z(1:izs)  zs z(izs+1:nz)];
   nz  = nz+1;
   izs = izs+1;
end

%-------------------------------------------------------------------------
%  Insert receiver depth into sound-speed profile.
%-------------------------------------------------------------------------               
                     
[dzr,izr] = min(abs(z-zr));
if((zr < zs) & (zr ~= z(izr))) izs = izs+1;end
if((zr < z(izr))) & ((zr > 0)) izr = izr-1;end
if(dzr ~= 0.)
   czr = c(izr);
   c   = [c(1:izr),czr,c(izr+1:nz)];
   z   = [z(1:izr), zr,z(izr+1:nz)];
   nz  = nz+1;
   izr = izr+1;
end

%-------------------------------------------------------------------------
%  Begin loop over requested raypaths. 
%-------------------------------------------------------------------------

for ipath=1:6

  if(paths(ipath) == 1)
  
  
%-------------------------------------------------------------------------  
%  Define the "image" SSP (nzi, zi, ci) and source & receiver indices
%  (izsi, izri) reflected about sea-surface and/or sea-floor for each
%  path requested.
%-------------------------------------------------------------------------
  
%-------------------------------------------------------------------------
%  Direct path 
%------------------------------------------------------------------------- 

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
      
      
%-------------------------------------------------------------------------
%  Surface reflection 
%------------------------------------------------------------------------- 

      if((ipath == 2) & (paths(ipath) == 1))
         ci   = [ fliplr(c) c(2:nz)];
         zi   = [-fliplr(z) z(2:nz)];
         nzi  = length(zi);
         izsi = (nz+1)-izs;
         izri = (nzi)-(nz-izr);
      end


%-------------------------------------------------------------------------
%  Bottom reflection 
%-------------------------------------------------------------------------

      if((ipath == 3) & (paths(ipath) == 1))
         ci   = [c fliplr(c(1:nz-1))];
         zi   = [z 2*z(nz)-fliplr(z(1:nz-1))];
         nzi  = length(zi);
         izsi = izs;
         izri = (nzi+1)-izr;
      end


%-------------------------------------------------------------------------
%  Surface-Bottom reflection 
%-------------------------------------------------------------------------

      if((ipath == 4) & (paths(ipath) == 1))
         ci   = [ fliplr(c) c(2:nz) fliplr(c(1:nz-1))];
         zi   = [-fliplr(z) z(2:nz) 2*z(nz)-fliplr(z(1:nz-1))];
         nzi  = length(zi);
         izsi = (nz+1)-izs;
         izri = (nzi+1)-izr;
      end


%-------------------------------------------------------------------------
%  Bottom-Surface reflection 
%-------------------------------------------------------------------------

      if((ipath == 5) & (paths(ipath) == 1))
         ci   = [c fliplr(c(1:nz-1)) c(2:nz)];
         zi   = [z 2*z(nz)-fliplr(z(1:nz-1)) 2*z(nz)+z(2:nz)];
         nzi  = length(zi);
         izsi = izs;
         izri = (nzi)-(nz-izr);
      end


%-------------------------------------------------------------------------
%  Surface-Bottom-Surface reflection 
%-------------------------------------------------------------------------

      if((ipath == 6) & (paths(ipath) == 1))
         ci   = [fliplr(c) c(2:nz-1) fliplr(c),c(2:nz)];
         zi   = [-fliplr(z) z(2:nz-1) 2*z(nz)- ...
                  fliplr(z) 2*z(nz)+z(2:nz)];
         nzi  = length(zi);
         izsi = (nz+1)-izs;
         izri = (nzi)-(nz-izr);
      end
          

%------------------------------------------------------------------------
%  Define "image" source & receiver depths (zsi, zri). Form SSP gradient 
%  (gi). Store SSP quantities for vector ray integrations. Compute 
%  harmonic mean sound speed (ch).
%------------------------------------------------------------------------

      zsi = zi(izsi);
      zri = zi(izri);

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
      z1      = zi(izstart:izend-1);
      z2      = zi(izstart+1:izend);
      dz      = z2-z1;

      c1sq = c1.^2;
      c1dz = c1.*dz;
      ch = abs(zi(izend)-zi(izstart))/sum(dz./c1);

%------------------------------------------------------------------------
%  Compute initial ray-parameter estimate (p0) from straight-line path
%  with constant sound speed equal to harmonic mean (ch). Call subroutine
%  rayint_direct to compute the ray range (r0), time (t0) and partial 
%  derivative (drdp) for p0.
%------------------------------------------------------------------------

    rp1(ipath) = p1*sum(c1dz./sqrt(1.-p1^2*c1sq));

    t(ipath) = sum(dz./(c1.*sqrt(1.-p1^2*c1sq)));

   end
   
end

return;
