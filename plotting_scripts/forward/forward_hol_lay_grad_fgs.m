%================================================================
%   FORWARD.M
%================================================================
function [E,Fm] = forward_hol_lay_grad_fgs(m,F,frbw);
global nlayfix nlay count;
%
% Setting up the environment in the (nlay+2)x4 Array geo_sin:
%
%
% Gradient:
%
nglay = 2;
nglay2 = 10;
ilay = 1;
for iglay=1:nglay

    dcdz = (m((iglay-1)*6+3)-m((iglay-1)*6+2))/(nglay2-1);
    drdz = (m((iglay-1)*6+5)-m((iglay-1)*6+4))/(nglay2-1);
    for i=1:nglay2

           lay(ilay,1) = (m((iglay-1)*6+1)/nglay2);
           lay(ilay,2) = m((iglay-1)*6+2)+dcdz*(i-1);
           lay(ilay,4) = m((iglay-1)*6+4)+drdz*(i-1);

           lay(ilay,3) = m((iglay-1)*6+6)/((m((iglay-1)*6+3)+m((iglay-1)*6+2))/2000);
           ilay = ilay + 1;
    end;
end;
nlay2 = nlay+(nglay*nglay2)-nglay;
iipar = [1 2 4 3];
idx = 6*nglay+1;
for ilay = (nglay*nglay2)+1:nlay2

    for ipar = 1:4

        lay(ilay,iipar(ipar)) = m(idx);

        % Conversion to dB/m/kHz
        if iipar(ipar) == 3

            lay(ilay,3) = m(idx)/(m(idx-2)/1000);

        end

        idx = idx + 1;
end; end;
% Halfspaces:
uhs = [NaN F(1).cw 0 F(1).rw];
lhs = [NaN m(end-2) m(end)/(m(end-2)/1000) m(end-1)];
geo_sin=[uhs; lay; lhs];
%save bla geo_sin;

%
% Compute reflectivity:
%
  E = 0;
  for ifreq=1:F(1).nfreq
    [ref] = ref_nlay_gfav3(F(ifreq).ang,geo_sin,F(1).freq(ifreq),frbw);% compute Reflection
%    ref         = -20*log10(abs(ref));
    ref         = abs(ref);
    Fm(ifreq).dat = reshape(ref,length(ref),1);
    F(ifreq).dat = reshape(F(ifreq).dat,length(F(ifreq).dat),1);
%    E = E + sum((Fm(ifreq).dat-F(ifreq).dat).^2)/(2*F(ifreq).sd^2);

    E = E + (Fm(ifreq).dat-F(ifreq).dat)'*F(ifreq).c_inv*...
       (Fm(ifreq).dat-F(ifreq).dat)/2;
    Rtmp(ifreq,:) = ref;

  end;

%Rtmp(ifreq+1,:) = F(1).ang;
%save('rep.txt','Rtmp','-ascii');

return;
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF
