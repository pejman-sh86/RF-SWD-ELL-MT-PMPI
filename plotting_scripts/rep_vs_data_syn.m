function rep_vs_data
%
%PLOT reflection model and data 
%
paper = 0;
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);
global c1 rho1 znorm lay_thick sz

sample = 'sample.dat';
plotfile1 = strrep(sample,'sample.dat','rep_vs_data.eps');
ml_file = strrep(sample,'sample.dat','ml_std.mat');

start = 2;	% start and step define the number of 
step = 3;	% data that are plotted (e.g. every 3rd datum)
load ./blde4_syn.mat;
load(ml_file);
aincr    =   3;		%
a_start  =   1;		%
nang_tmp =  45;		% Nr of angles
fincr    =   1;		%
f_start  =   6;		%
nfreq    =   8;		% Nr of frequencies
freq = [315 400 500 630 800 1000 1250 1600];
for i = 1:nfreq
  nang(i) = length(F(i).ang);
end
%
% Forward model parameters:
%
c1 = 1511; rho1 = 1.029;
lay_thick = 0.1;    % discretization for layer thickness - controls 
                      % upper frequency limit
znorm=lay_thick/2:lay_thick:1-lay_thick/2;
sz=length(znorm);
%
% Model:
% h  rhot  rhob  nu  ct  cb  alpha
mlm = [1.993618      1.3172073      1.5574454 0.79358805        1472.13      1462.2936     0.35067362];

%
% REAL DATA
%
ml_sd = stdv_dB(:,2)

[Frep] = forward(mlm,F,freq,nfreq,nang);

figure(1);
if(paper == 1)
  llw = [0.10 0.40 0.7  0.10 0.40 0.7  0.10 0.40 0.7 ];
  llh = [0.69 0.42 0.15];
else
  llw = [0.15 0.43 0.71 0.15 0.43 0.71 0.15 0.43 0.71];
  llh = [0.69 0.42 0.15];
end
i = 1;
j = 1;
for k=1:nfreq
   if k <= 3
     subplot('Position',[llw(i) llh(1) 0.26 0.26]);
     i = i + 1;
   elseif k <= 6
     subplot('Position',[llw(i) llh(2) 0.26 0.26]);
     i = i + 1;
   else
     subplot('Position',[llw(i) llh(3) 0.26 0.26]);
     i = i + 1;
   end
   hold on;box on
   Err = ml_sd(k)*ones(nang(k),1);
   plot(F(k).ang(start:step:end), F(k).dat(start:step:end),'k.');
   errorbar(F(k).ang(start:step:end),F(k).dat(start:step:end),Err(start:step:end),'k.');
   hold on;
   
   if(paper == 1)
     plot(F(k).ang(start:step:end),Frep(k).dat(start:step:end),'k');
   else
     plot(F(k).ang(start:step:end),Frep(k).dat(start:step:end),'b');
   end
%   title([num2str(freq(k)) ' Hz'] )
   text(55,40,[num2str(freq(k)) ' Hz'])
   axis([0 90 0 50]); 
   if(paper == 1)
     set(gca,'XTick',[0 30 60 90],'YTick',[0:20:40],...
     'YLim',[10 45],'FontSize',12);
   else
     set(gca,'XTick',[0 40 80],'YTick',[0:20:40],...
     'YLim',[10 45],'FontSize',16);
   end
   if k >=  6;xlabel('Angle [deg]');end
   if k == 1;ylabel('BL [dB]');end
   if k == 4;ylabel('BL [dB]');end
   if k == 7;ylabel('BL [dB]');end
   if(k < 6) 
     set(gca,'XTickLabel',[]);
   else 
     if(paper == 1)
       set(gca,'XTickLabel',[0 30 60 90]);
     else
       set(gca,'XTickLabel',[0 40 80]);
     end
   end
   if k == 1; set(gca,'YTickLabel',[0 20 40]);end
   if k == 2; set(gca,'YTickLabel',[]);end
   if k == 4; set(gca,'YTickLabel',[0 20 40]);end
   if k == 3; set(gca,'YTickLabel',[]);end
   if k == 5; set(gca,'YTickLabel',[]);end
   if k == 6; set(gca,'YTickLabel',[]);end
   if k == 7; set(gca,'YTickLabel',[0 20 40]);end
   if k == 8; set(gca,'YTickLabel',[]);end
   if k == 9; set(gca,'YTickLabel',[]);end
end
saveas(gca,plotfile1,'epsc2');

%===========================================================================
%   FORWARD.M
%===========================================================================
function [Fm] = forward(m,F,freq,nfreq,nang);
  
%--------------- GLOBAL VARIABLES ------------------------------------------
global c1 rho1 znorm lay_thick sz
global minlim maxlim
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
