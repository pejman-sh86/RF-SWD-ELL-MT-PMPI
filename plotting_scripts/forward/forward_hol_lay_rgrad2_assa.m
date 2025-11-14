%================================================================
%   FORWARD.M
%================================================================
function [E] = forward_hol_lay_rgrad_assa(m,F,frbw);
global nlay nlayfix fix_hs;
%
% Setting up the environment in the (nlay+2)x4 Array geo_sin:
%
nfreq = length(F(1).freq);

%
% Gradient:
%
nfpgl = 4;
nglay = 1;
nglay2 = 10;
ilay = 1;
for iglay=1:nglay

    drdz = (m((iglay-1)*nfpgl+7)-m((iglay-1)*nfpgl+3))/(nglay2-1);
    for i=1:nglay2

           lay(ilay,1) = (m((iglay-1)*nfpgl+1)/nglay2);
           lay(ilay,2) = m((iglay-1)*nfpgl+2);
           lay(ilay,4) = m((iglay-1)*nfpgl+3)+drdz*(i-1);
    
           lay(ilay,3) = m((iglay-1)*nfpgl+nfpgl)/((m((iglay-1)*nfpgl+2))/1000);
           ilay = ilay + 1;
    end;
end;
nlay2 = nlay+(nglay*nglay2)-nglay;
iipar = [1 2 4 3];
idx = nfpgl*nglay+1;
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
%save bla geo_sin m;

%return;
%
% Compute reflectivity:
%
if((m > F(1).minlim) & (m < F(1).maxlim))
  E = 0;
  for i=1:nfreq
    [ref] = ref_nlay_gfav3(F(i).ang,geo_sin,F(1).freq(i),frbw);% compute Reflection
%    ref         = -20*log10(abs(ref));
    ref         = abs(ref);
    Fm(i).dat = reshape(ref,length(ref),1);
    F(i).dat = reshape(F(i).dat,length(F(i).dat),1);
    E = E + length(find(F(1).Rex(:,1))~=0)/2*log(sum(((Fm(i).dat-F(i).dat).*F(1).Rex(:,i)).^2));
%    E = E + sum((Fm(i).dat-F(i).dat).^2)/(2*F(i).sd^2);
    Rtmp(i,:) = ref;
  end
%  E = E/10^3;
else
  disp('Out of bounds')
  for i=1:nfreq
    [ref] = ref_nlay3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
%    ref         = -20*log10(abs(ref));
    ref         = abs(ref);
    Fm(i).dat = ref';
  end
  E = 10^50.;
end

%Rtmp(i+1,:) = F(1).ang;
%save('replica.txt','Rtmp','-ascii');
return;
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF
