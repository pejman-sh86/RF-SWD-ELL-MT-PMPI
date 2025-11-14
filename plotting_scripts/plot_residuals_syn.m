function plot_residual_syn
%
% Estimates sd from ML model (output of ASSA) and calculates
% residuals of ML model and real data.
%
% 1) Saves residuals.dat
%

paper = 0;
set(0, 'DefaultFigurePaperPosition', [0 0 6 6]);
global c1 rho1 znorm lay_thick sz

sample = 'syn_sample.dat';
resfile = strrep(sample,'sample.dat','residuals0.mat');
plotfile1 = strrep(sample,'sample.dat','residuals0.eps');
plotfile2 = strrep(sample,'sample.dat','autocorr.eps');

npar = 7;
xml = [1.6052915      1.3276892      1.4473407 0.77455655      1472.8899      1474.7546     0.30502054];
%xml = [1.881867      1.3199531      1.5188286     0.81879086      1472.5972 1466.5915     0.33494602];
%xml = [ 1.8656428      1.3327651      1.5171949     0.95574165      1472.9799 1460.1604     0.28694253];
%xml = [1.8380438       1.319941      1.5167713     0.71579598      1470.7303 1460.8863     0.30391302];
%xml = [1.3770171      1.3396415      1.4722703     0.95119215       1470.546 1468.9141     0.31469968];
%xml = [1.9124603      1.3106415      1.5160627     0.70002737      1471.8001 1467.6332     0.32856802];
%xml = [1.9 1.32 1.50 0.80 1472.83  1465.87 0.30];

load ../blde4_syn.mat;
nfreq    =  8;  % Nr of frequencies
%
% Forward model parameters:
%
c1 = 1511; rho1 = 1.029;

%---------------------------------------------------------------------------
%  CODE STARTS
%---------------------------------------------------------------------------
lay_thick = 0.1;    % discretization for layer thickness - controls 
                      % upper frequency limit
znorm=lay_thick/2:lay_thick:1-lay_thick/2;
sz=length(znorm);

%---------------------------------------------------------------------------
%  Load synthetic or real data.
%---------------------------------------------------------------------------

freq = F(1).freq;
for i = 1:nfreq
  F(i).c_inv = inv(F(i).csave);
  nang(i) = length(F(i).ang);
end

nang

[Frep] = forward(xml,F,freq,nfreq,nang);

ml_sd = zeros(nfreq,1);
for i = 1:nfreq
  ml_sd(i) = sqrt(1/(nang(i)-npar) * sum((Frep(i).dat-...
          F(i).dat).^2));
end

stdv_dB = [freq' ml_sd]

save ml_std.mat stdv_dB;

for i = 1:nfreq
  F(i).diff = (Frep(i).dat-F(i).dat);
  B(i).diff = F(i).diff;
end

save(resfile,'B');

llw = [0.10 0.40 0.7  0.10 0.40 0.7  0.10 0.40 0.7 ];
llh = [0.69 0.42 0.15];
figure(1);
for k=1:nfreq
   if k <= 3
     subplot('Position',[llw(k) llh(1) 0.26 0.26]);
     i = i + 1;
   elseif k <= 6
     subplot('Position',[llw(k) llh(2) 0.26 0.26]);
     i = i + 1;
   else
     subplot('Position',[llw(k) llh(3) 0.26 0.26]);
     i = i + 1;
   end
   hold on;box on
   plot([0 90],[0 0],'-k');
   plot(F(k).ang, F(k).diff/ml_sd(k),'kx');

   text(52,2.3,[num2str(freq(k)) ' Hz'])
   axis([0 90 -3 3]);
   if k==6;xlabel('Angle [deg]');end
   if k>6;xlabel('Angle [deg]');end
   if k == 1;ylabel('Residual [\sigma]');end
   if k == 4;ylabel('Residual [\sigma]');end
   if k == 7;ylabel('Residual [\sigma]');end
   set(gca,'XTick',[0 30 60 90]);
   set(gca,'XTickLabel',[]);
   set(gca,'YTick',[-2 0 2]);
   set(gca,'YTickLabel',[]);
   if k>=6 ; set(gca,'XTickLabel',[0 30 60 90]);end;
   if k==1 ; set(gca,'YTickLabel',[-2 0 2]);end;
   if k==4 ; set(gca,'YTickLabel',[-2 0 2]);end;
   if k==7 ; set(gca,'YTickLabel',[-2 0 2]);end;

   set(gca,'YGrid','on');
end

saveas(gca,plotfile1,'epsc2');

%
% Autocorrelation of residuals
%
gc2 = figure(2);
if(paper == 1)
  width = 0.26;
  height = 0.26;
  llw = [0.10 0.40 0.7  0.10 0.40 0.7  0.10 0.40 0.7 ];
  llh = [0.75 0.44 0.13];
else
  width = 0.26;
  height = 0.24;
  llw = [0.10 0.40 0.70 0.10 0.40 0.70 0.10 0.40 0.70];
  llh = [0.75 0.44 0.13];
end
for k=1:nfreq

  hold on;
  if k <= 3
    subplot('Position',[llw(k) llh(1) width height]);
  elseif k <= 6
    subplot('Position',[llw(k) llh(2) width height]);
  else
    subplot('Position',[llw(k) llh(3) width height]);
  end
  Aresres = xcov(F(k).diff)/nang(k);
  lag = -floor(length(Aresres)/2):1:floor(length(Aresres)/2);
  plot(lag,Aresres/max(Aresres),'-k.');
  text(length(Aresres)/6,0.9,[num2str(freq(k)) ' Hz'],'FontSize',8);
  set(gca,'FontSize',10);
  set(gca,'XLim',[lag(1) lag(end)],'YLim',[-0.5 1.1]);
  if(k>=6);xlabel('Lag');end;
  if ((k==1) | (k==4) | (k==7)); ylabel('A_{nn}');end;
  set(gca,'XTick',[-100 -50 0 50 100],'XTickLabel',{'-100';'-50';'0';'50';'100'});
  set(gca,'YTick',[-0.5 0 0.5 1],'YTickLabel',{'';'';'';''});
  if((k==1) | (k==4) | (k==7) ); set(gca,'YTickLabel',{'';'0';'';'1'});end;
%
% Perform run test:
%
%   z = runtest(tmp2')
end
saveas(gc2,plotfile2,'epsc2');
%figure(3);
%for k=1:nfreq
%   mxlag = floor((nang(k)-1)/2);
%   [acov,lag]= xcov(F(k).diff,mxlag,'coeff');
%   [acov,lag]= xcov(F(k).diff);
%   subplot(3,3,k);hold on;box on
%   plot(lag, acov,'-xk');
%   hold on;
%   
%   title([num2str(freq(k)) ' Hz'] )
%   axis([0 90 -5 5]); set(gca,'Xtick',[0:30:90]);
%   if k>6;xlabel('Lag');end
%   set(gca,'YGrid','on');
%
% Perform run test:
%
%   z = runtest(tmp2')
%end

%===========================================================================
%   FORWARD.M
%===========================================================================

function [Fm] = forward(m,F,freq,nfreq,nang);

%--------------- GLOBAL VARIABLES ------------------------------------------
global c1 rho1 znorm lay_thick sz
%----------------------------------------------------------------------------
%
% Setting up the environment in the (sz+2)x4 Array geo_sin:
%

% rhos = rhot + sin(znorm*pi/2).^no*(rhob-rhot)
rhos = m(2) + sin(znorm*pi/2).^m(4)*(m(3)-m(2));
% cs=ct+(cb-ct)*znorm;
cs=m(5)+(m(6)-m(5))*znorm;

geo_sin=[NaN c1 0 rho1; ...
lay_thick*m(1)*ones(sz,1) cs' m(7)*ones(sz,1) rhos'; ...
        NaN m(6) m(7) m(3)];
%
% Compute reflectivity:
%
for i=1:nfreq
  [ref] = ref_nlay3(F(i).ang,geo_sin,freq(i));% compute Reflection
  Fm(i).dat = ref';
end 

return

% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
