function plot_residual_hol
%
% Estimates sd from ML model (output of ASSA) and calculates
% residuals of ML model and real data.
% CODE ALWAYS SAVES ABSOLUTE RESIDUALS (IN dB)!!
% BUT PLOTS RESIDUALS IN STDEVs
%
% 1) Saves residuals.dat
%

isd = 1;	% Scale plot to sds? 1 => yes

set(0, 'DefaultFigurePaperPosition', [0 0 6 6]);
global c1 rho1 znorm lay_thick sz

plotfile1 = 'site20b_res0.eps';
resfile = 'site20b_res0.mat';

%plotfile1 = 'residuals_1.eps';
%resfile = 'residuals_1.mat';

%load('site19_mapdata.mat');
load('site20b_rep0.mat');
Frep = F;
load('hol_site20b.mat');
nfreq = length(F(1).freq);
nang = length(F(1).ang);
nsubpfig = 12;
nfig = ceil(nfreq/nsubpfig);
if nfreq < nsubpfig; nsubpfig = nfreq;end;

%clear B;

for i = 1:nfreq
  F(i).diff = (Frep(i).dat-F(i).dat);
  B(i).diff = F(i).diff;
  ml_sd(i) = sqrt(1/(nang-F(1).nmod) * sum((Frep(i).dat-...
             F(i).dat).^2));
  if(isd == 1)
    F(i).diff = F(i).diff/ml_sd(i);
  end
end

save(resfile,'B');

llw = [0.05 0.285 0.52 0.755 0.05 0.285 0.52 0.755 ... 
       0.05 0.285 0.52 0.755];
llh = [0.69 0.42 0.15];

spw = 0.195
sph = 0.26

kdat = 1;
for ifig = 1:nfig

figure(ifig);hold on;
i = 1;
for k=1:nsubpfig
%   subplot(3,3,k);
   if k <= 4
     subplot('Position',[llw(k) llh(1) spw sph]);
     i = i + 1;
   elseif k <= 8
     subplot('Position',[llw(k) llh(2) spw sph]);
     i = i + 1;
   else
     subplot('Position',[llw(k) llh(3) spw sph]);
     i = i + 1;
   end
   hold on;box on
   plot(F(kdat).ang, F(kdat).diff,'kx');
   plot(0,0,0,90,'k');
   hold on;
   
%   title([num2str(freq(k)) ' Hz'] )
%   text(52,2.3,[num2str(freq(k)) ' Hz'])
   axis([16 86 -2 2]);
%   if k==6;xlabel('Angle [deg]');end
%   if k>6;xlabel('Angle [deg]');end
%   if k == 1;ylabel('Residual [\sigma]');end
%   if k == 4;ylabel('Residual [\sigma]');end
%   if k == 7;ylabel('Residual [\sigma]');end
%   set(gca,'XTick',[0 30 60 90]);
%   set(gca,'XTickLabel',[]);
%   set(gca,'YTick',[-2 0 2]);
%   set(gca,'YTickLabel',[]);
%   if k>=6 ; set(gca,'XTickLabel',[0 30 60 90]);end;
%   if k==1 ; set(gca,'YTickLabel',[-2 0 2]);end;
%   if k==4 ; set(gca,'YTickLabel',[-2 0 2]);end;
%   if k==7 ; set(gca,'YTickLabel',[-2 0 2]);end;

   set(gca,'YGrid','on');
  kdat = kdat+1;
  if kdat>nfreq;break;end;
end;end;

saveas(gca,plotfile1,'epsc2');

return;

% ---------------------------------------------------------------------
% ...this is the end my fiend.
% EOF   
