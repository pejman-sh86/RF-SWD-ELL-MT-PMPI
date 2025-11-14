function plot_residual
%
% Estimates sd from ML model (output of ASSA) and calculates
% residuals of ML model and real data.
% CODE ALWAYS SAVES ABSOLUTE RESIDUALS (IN dB)!!
% BUT PLOTS RESIDUALS IN STDEVs
%
% 1) Saves residuals.dat
%

dex = 0;	% Are we looking at dexp statistics here?
		%
i_f_s = 0;	% Save first residuals (first ASSA run)
		% or are we saving the residuals for second
		% ML model (second, corrected ASSA run?
isd = 1;	% Scale plot to sds? 1 => yes

set(0, 'DefaultFigurePaperPosition', [0 0 6 6]);
global c1 rho1 znorm lay_thick sz

sample = 'real2_sample.dat';
if dex == 0
  if i_f_s == 0
    mlfile = strrep(sample,'sample.dat','ml_std0.mat');
    plotfile1 = strrep(sample,'sample.dat','residuals0.eps');
    resfile = strrep(sample,'sample.dat','residuals0.mat');
  elseif i_f_s == 1
    mlfile = strrep(sample,'sample.dat','ml_std1.mat');
    plotfile1 = strrep(sample,'sample.dat','residuals1.eps');
    resfile = strrep(sample,'sample.dat','residuals1.mat');
  else
    mlfile = strrep(sample,'sample.dat','ml_std2.mat');
    plotfile1 = strrep(sample,'sample.dat','residuals2.eps');
    resfile = strrep(sample,'sample.dat','residuals2.mat');
  end
else
  if i_f_s == 0
    mlfile = strrep(sample,'sample.dat','ml_std_dexp0.mat');
    resfile = strrep(sample,'sample.dat','residuals_dexp0.mat');
    plotfile1 = strrep(sample,'sample.dat','residual_dexp0.eps');
  elseif i_f_s == 1
    mlfile = strrep(sample,'sample.dat','ml_std_dexp1.mat');
    resfile = strrep(sample,'sample.dat','residuals_dexp1.mat');
    plotfile1 = strrep(sample,'sample.dat','residual_dexp1.eps');
  else
    mlfile = strrep(sample,'sample.dat','ml_std_dexp2.mat');
    resfile = strrep(sample,'sample.dat','residuals_dexp2.mat');
    plotfile1 = strrep(sample,'sample.dat','residual_dexp2.eps');
  end
end

npar = 7;
if i_f_s == 0
  % FIRST ASSA RUN:
  xml = [2.87679      1.3428398      1.5435234     0.83889602      1475.3855 1492.5991     0.29226766];
%  xml = [1.8698263 1.3488856 1.5070618 0.8012752 1472.836 1465.8798 0.30029256];
% DEXP
%  xml = [1.8986634 1.3569792 1.4606326 0.90246767 1472.8487 1466.3762 0.24461243];
elseif i_f_s == 1
  % FIRST ASSA RUN:
  xml = [3.0958466      1.4342324      1.2828879 0.18719779      1475.5362      1476.4563     0.13225876];
% DEXP
%  xml = [1.544594 1.3673591 1.4810022 1.2538181 1472.8164 1468.3417 0.30431891];
else
  % CORRECTED (SECOND) ASSA RUN:
  xml = [];
% DEXP
%  xml = [1.6850218 1.3625658 1.5272122 1.0693537 1472.5125 1467.5807 0.35709187];
end

load hol_ida;

freq = F(1).freq;
nfreq = length(freq);
F(1).lay_thick = 0.1; % discretization for layer thickness - controls
                      % upper frequency limit
F(1).znorm=F(1).lay_thick/2:F(1).lay_thick:1-F(1).lay_thick/2;
F(1).sz=length(F(1).znorm);

%---------------------------------------------------------------------------
%  Load synthetic or real data.
%---------------------------------------------------------------------------

[E, Frep] = forward_hol_grad(xml,F);

ml_sd = zeros(nfreq,1);
for i = 1:nfreq
  ml_sd(i) = sqrt(1/(length(F(i).ang)-npar) * sum((Frep(i).dat-...
          F(i).dat).^2));
end

stdv_dB = [freq ml_sd]

save(mlfile,'stdv_dB');
%save bla F Frep;

for i = 1:nfreq
  F(i).diff = (Frep(i).dat-F(i).dat);
  B(i).diff = F(i).diff;
  B(i).ang = F(i).ang;
  if(isd == 1)
    F(i).diff = F(i).diff/ml_sd(i);
  end
end

save(resfile,'B');

llw = [0.10 0.40 0.7  0.10 0.40 0.7  0.10 0.40 0.7 ];
llh = [0.69 0.42 0.15];

figure(1);
for k=1:nfreq
%   subplot(3,3,k);
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
   plot(F(k).ang, F(k).diff,'kx');
   plot(0,0,0,90,'k');
   hold on;
   
%   title([num2str(freq(k)) ' Hz'] )
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

return;

% ------------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
