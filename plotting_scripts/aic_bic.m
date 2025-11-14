%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [minL,BIC,AIC,C] = aic_bic(filebase,len,len2);

ifirst = 0;
isd = 1;
resfile = strcat(filebase,'_res.txt');
data    = strcat(filebase(1:len2),'.mat');
sdfile  = strcat(filebase(1:len),'exstd.txt');
if(isd == 0)
   disp('WARNING: RUNNING EXPERIMENTAL ERRORS');
   sd = load(sdfile);
else
   disp('WARNING: ASSUMING UNKNOWN STD-DEV');
end;


res = load(resfile);
ii=1;

load(data);
C(ii).file = filebase;

M  = length(F(1).minlim)	% Number of model parameters
N  = size(F(1).dat,1);		% Number of data per frequency
Nf = length(F(1).freq);		% Number of data per frequency

%
% Compute misfit:
% --------------------------------------------------------------------------------
res = F(1).Rex(:,1:Nf)'.*res(1:Nf,:);

%
% Normalization for f(m_i)
% --------------------------------------------------------------------------------
fnorm(ii) = 0;
fnorm2(ii) = 0;
minL(ii) = 0.;
minL2(ii) = 0.;
for i = 1:length(F(1).freq)
  
  idx = find(F(1).Rex(1:N,i) == 1);
  N2(i) = length(idx);

  if(isd == 0)
    minL(ii) = minL(ii) + sum(res(i,:).^2/(2*sd(i)^2.));
    fnorm(ii) = fnorm(ii) - N2(i)/2*log(2*pi) - N2(i)*log(sd(i));
  else
    sd2(i) = sqrt(1/(N2(i)-(M/Nf))*sum(res(i,:).^2));
    minL(ii) = minL(ii) + N2(i)/2*log(sum(res(i,:).*res(i,:),2));
    fnorm(ii) = fnorm(ii) - N2(i)/2*log(2*pi) - N2(i)*log(sd2(i));
  end;

end
%minL(ii)

%
% Bayesian Information Criterion (for comparison)
% --------------------------------------------------------------------------------
if(isd == 0)
   BIC(ii) = -2 * (fnorm(ii) - minL(ii)) + M * log(sum(N2));
else
   BIC(ii) = -2 * (fnorm(ii) - minL(ii)) + M * log(sum(N2));
end;
%
%  Akaike Information Criterion (for comparison)
% --------------------------------------------------------------------------------
if(isd == 0)
   AIC(ii) = -2 * (fnorm(ii) - minL(ii)) + 2 * M;
else
   AIC(ii) = -2 * (fnorm(ii) - minL(ii)) + 2 * M;
end;

fprintf(1,'%s,\t %4i, %12.4f, %12.4f, %4i, %12.4f \n',filebase, M,minL(ii),M*log(sum(N2)),Nf,fnorm2(ii));

return;
