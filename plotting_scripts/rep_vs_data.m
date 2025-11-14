function rep_vs_data
%
%PLOT reflection model and data 
%
dex = 0;
paper = 0;
set(0, 'DefaultFigurePaperPosition', [0 0 7 7]);
global c1 rho1 znorm lay_thick sz

sample = 'real2_sample.dat';
if dex == 0
  plotfile1 = strrep(sample,'sample.dat','rep_vs_data3.eps');
  ml_file = strrep(sample,'sample.dat','ml_std1.mat');
else
  plotfile1 = strrep(sample,'sample.dat','rep_vs_data_dexp.eps');
  ml_file = strrep(sample,'sample.dat','ml_std_dexp.mat');
end

load ../blde4_2_ping_ida_inv.mat;
%load ./blde4_syn.mat;
load(ml_file);

fid = fopen('./parameters.txt');
A = fscanf(fid,'%11c%3i\n');
A([12,24,36,48,60,72]);

aincr   = A(12);          %
a_start = A(24);          %
nang_tmp= A(36);          % Nr of angles
fincr   = A(48);          %
f_start = A(60);          %
nfreq   = A(72);  % Nr of frequencies
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
mlm = [1.6665  1.3563  1.5461  1.0035  1472.7412  1465.8553  0.3799];
%mlm = [3.0670062      1.3937878      1.2539207 0.3421364      1477.0239      1470.4547     0.10347888];
%mlm = [2.0650185      1.3553093 1.4392162     0.85836735      1473.0071 1465.6735     0.21551576];
%mlm = [1.9486875      1.3543391      1.4835551 0.85430818       1472.925 1465.4089      0.2609109];
%mlm = [2.0930214      1.3508944      1.4841878     0.73839097 ...
%       1472.8123 1467.5688     0.25152184];

%
% REAL DATA
%
ml_sd = stdv_dB(:,2)
ang = zeros(1,nang_tmp);
freq = zeros(1,nfreq);

j = a_start;
for i = 1:nang_tmp
  ang(i) = xde4.ang(j);
  j = j + aincr;
end
j = f_start;
for i = 1:nfreq
  freq(i) = xde4.pref(j,1);
  j = j + fincr;
end
tmp = xde4.bl(f_start:fincr:f_start+(fincr*nfreq)-1,a_start:aincr:...
              a_start+(aincr*nang_tmp)-1);
for i = 1:nfreq
  jj = 1;
  for j = 1:nang_tmp
    if(isnan(tmp(i,j)) ~= 1)
      dat(jj) = tmp(i,j);
      tmp_ang(jj) = ang(j);
      jj = jj+1;
    end
  end
  F(i).dat = dat';
  nang(i) = jj-1;
  F(i).ang = tmp_ang;
end
disp('I work on real data now!');

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
   stp=3;
   Err = ml_sd(k)*ones(nang(k),1);
   plot(F(k).ang(1:stp:end), F(k).dat(1:stp:end),'k.');
   errorbar(F(k).ang(1:stp:end),F(k).dat(1:stp:end),Err(1:stp:end),'k.');
   hold on;
   
   if(paper == 1)
     plot(F(k).ang,Frep(k).dat,'k');
   else
     plot(F(k).ang,Frep(k).dat,'b');
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
