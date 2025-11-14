function [ref] = forward_ref_nlay3(m,angobs,bands);

NAVEF = 1;
IGA = 1;
FBW = 50;
frbw = 1/15;
ispher = 0;
cw = 1512;
rw = 1.029;
      model = m;
      NFP=length(m);
      k = (NFP-3)/4;
      clear geo vp alf1dB refl r_ave;
      geo = zeros(k+2,4);
      geo(1,:) = [NaN, cw, 0., rw];
      for i=1:k;
         geo(i+1,:) = [model((i-1)*4+1),model((i-1)*4+2),...
                       model((i-1)*4+4),model((i-1)*4+3)];
      end;
      geo(k+2,:) = [NaN, model(NFP-2),model(NFP),...
                       model(NFP-1)];
%      k(l) = 2;
%      geo = [ NaN, 1511., 0., 1.029;...
%              0.633004308, 1593.9545,  0.61620242,  1.97112622;...
%              0.144615054, 1601.31862, 0.618422141, 2.01869688;...
%              NaN, 1589.17373, 0.334276375, 2.00589373];
%      afdep(l) = 0.584093125;

      for iband=1:length(bands);
         if(NAVEF > 1)
            if(IGA == 1)
               flo = bands(iband) - bands(iband)*frbw;
               fhi = bands(iband) + (iband)*frbw;
            elseif(IGA == 2)
               flo = bands(iband) / (2.^(1./6.)); % Use 1/3 octave
               fhi = bands(iband) * (2.^(1./6.)); %
            else
               flo = bands(iband) - FBW;
               fhi = bands(iband) + FBW;
            end
            fstep = (fhi-flo)/(NAVEF-1);
            fr = flo + ([1:NAVEF]-1) .* fstep;
         else
            fr = bands(iband);
         end
         if(ispher == 0)
            [refl]=ref_nlay3(angobs,geo,fr);
            refl = abs(refl);
         else
            [refl, Rp] = spherical_refl(geo,z_t,fr,angobs);
            refl = abs(refl);
         end;
         %% Do freq average
         if(NAVEF > 1)
%            if(IGA == 1)
%               refl = abs(refl);
%               nang = length(angobs);
%               df = mean(diff(fr));
%               f0 = fr;
%               for iang = 1:nang
%                  r_ave(iang) = sum(refl(:,iang) .* ...
%                                exp(-(fr'-f0).^2/(frbw*f0)^2)*df)/...
%                                sum(exp(-(fr'-f0).^2/(frbw*f0)^2)*df);
%               end;
%            else
               r_ave = sqrt(sum(refl.*refl,1)/length(fr));
%            end;
         else
            r_ave = abs(refl);
         end;
         ref(:,iband) = r_ave;
      end;

return;
