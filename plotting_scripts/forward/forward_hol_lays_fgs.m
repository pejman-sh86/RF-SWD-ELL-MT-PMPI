%================================================================
%   FORWARD.M
%================================================================
function [E,Fm] = forward_hol_lays_fgs(m,F,frbw);
global nlayfix nlay count;
%
% Setting up the environment in the (nlay+2)x4 Array geo_sin:
%
uhs = [NaN F(1).cw 0 F(1).rw];
lhs = [NaN m(end-2) m(end)/(m(end-2)/1000) m(end-1)];

idx = 1;
iipar = [1 2 4 3];
if nlayfix > 0
    maxlay = nlayfix+nlay;
else
    maxlay = nlay;
end
for ilay = 1:maxlay

    for ipar = 1:4

        lay(ilay,iipar(ipar)) = m(idx);

        if iipar(ipar) == 3

            lay(ilay,3) = m(idx)/(m(idx-2)/1000);

        end

        idx = idx + 1;

    end; end;

if (nlayfix > 0)
    for ipar = 1:(nlayfix*4)
        if(ipar <= (nlayfix*4))
            if(m(ipar)>F(1).lim(:,ipar))
                idx = find(m(ipar)>F(1).lim(:,ipar));
            else
                idx = 1;
                count = count + 1;
            end;
            prior(ipar) = F(1).n1(idx(end),ipar);
        else
            if(m(ipar+1)>F(1).lim(:,ipar))
                idx = find(m(ipar+1)>F(1).lim(:,ipar));
            else
                idx = 1;
                count = count + 1;
            end;
            prior(ipar) = F(1).n1(idx(end),ipar);
        end
    end

    P = prod(prior);

end;

geo_sin=[uhs; lay; lhs];
%save bla geo_sin;
%return;

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

  end;
  if (nlayfix > 0)
      E = E - log(P);
  end;

return;
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF
