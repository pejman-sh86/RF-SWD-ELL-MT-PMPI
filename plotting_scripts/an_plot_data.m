function [] = an_plot_data()

set(gcf, 'renderer', 'painters')
%set(gcf, 'renderer', 'OpenGL')
%set(gcf, 'renderer', 'zbuffer')
set(0, 'DefaultFigurePaperPosition', [0 0 8 4.8]);
sample = 'sample.dat';
plotfile1 = strrep(sample,'sample.dat','real_data.eps');
plotfile2 = strrep(sample,'sample.dat','rep_data.eps');
nfreq = 43;
nang  = 91;

load V0Th0.mat;    % Loads Theta0(91)    -- angles
                   %       freq(43)      -- frequencies
                   %       V0(43x91)     -- reflectivities
                   %       Weight(43x91) -- arbitrary weights vs freq & angle
id = 1;
Weight = Weight / max(max(Weight));
for ifreq = 1:nfreq
  for iang = 1:nang
    if(Weight(ifreq,iang) > 0.); % Selects points.
      Weight(ifreq,iang) = 1;
      dat(id) = V0(ifreq,iang);
      id = id+1;
    else
      Weight(ifreq,iang) = 0;
    end
  end
end
ndat = id-1;
dat  = dat';

[Bft] = Beams([.5:.5:16],freq,1500,1024); % Calculate beam information

% model parameters: h, c1,rho1,alpha1, c2,rho2,alpha2
m_ml = [0.92976574      1529.1586      1.2151128 0.12781081 1649.2394 2.5     0.99962241]';
[VB,t1] = forward(m_ml,freq,nfreq,nang,Bft);

gc1 = figure(1);
%subplot(2,1,1);
%colormap(winter);
pcolor(Theta0,freq,-20*log10(V0));
%title('R, real data','FontSize',14);
set(gca,'FontSize',14);
%set(gca,'XTick',20,'XTickLabel',[]);
xlabel('Angle [deg]','FontSize',16);
ylabel('Frequency [Hz]','FontSize',16);
shading interp;
caxis([0 15]);
%caxis([0 1]);
gc2 = colorbar;
set(gc2,'FontSize',14);
text(110,820,'BL [dB]','Rotation',90,'FontSize',16)

gc2 = figure(2);
%subplot(2,1,2);
%colormap(winter);
pcolor(Theta0,freq,-20*log10(VB));
%pcolor(Theta0,freq,VB);
%%title('R, real data','FontSize',14);
set(gca,'FontSize',14);
%set(gca,'XTick',20,'XTickLabel',[]);
xlabel('Angle [deg]','FontSize',16);
ylabel('Frequency [Hz]','FontSize',16);
shading interp;
caxis([0 15]);
%caxis([0 1]);
gc2 = colorbar;
set(gc2,'FontSize',14);
text(110,820,'BL [dB]','Rotation',90,'FontSize',16)

saveas(gc1,plotfile1,'epsc2');
saveas(gc2,plotfile2,'epsc2');

return;

%===========================================================================
%   FORWARD.M
%===========================================================================

function [VB,t1] = forward(m,freq,nfreq,nang,Bft);

global fstart fend astart aend sd

% 2-Layer model
h1  = [   0 m(1)     ];
c1  = [1500 m(2) m(5)];
r1  = [   1 m(3) m(6)];
a1  = [   0 m(4) m(7)];

%  3-Layer model
%h1  = [   0 m(1) m(5)     ];
%c1  = [1500 m(2) m(6) m(9)];
%r1  = [   1 m(3) m(7) m(10)];
%a1  = [   0 m(4) m(8) m(11)];

[t1,V] = arb_layer(c1,h1,r1,a1,freq); 

V2 = fliplr(-20*log10(abs(V)));
t2 = fliplr(round(t1));
ship = 0.;

%
% Simulate smearing of beam forming
%
[RLcalc,angle,VB,newang] = Smudge(freq,t2,V2,1500,1500,1500,ship,Bft,'yes',t2);
VB = 10.^(-VB/20);

return;

