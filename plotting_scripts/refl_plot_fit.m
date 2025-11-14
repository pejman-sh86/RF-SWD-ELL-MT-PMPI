%function [] = refl_plot_fit(filebase,NBAND);
%filebase = 'blde4_b_bl_ping_id2_nb2_15deg_400_1600';
filebase = 'blge5_3_bl_ping_id1_7oct_verbose';
%filebase = 'blge5_3_bl_ping_id1_7oct_verbose_05Ang';
%filebase = 'blge5_3_bl_ping_id1_IvgsEven_5Fres1d5';
NBAND = 6;
obsfile = strcat(filebase,'.txt');
prefile = strcat(filebase,'_rep.dat');

tmp1=dlmread(obsfile);
pre=dlmread(prefile);

obs = tmp1(5:4+NBAND,:);
ang = tmp1(5+NBAND,:);
Rex = tmp1(6+NBAND:5+NBAND+NBAND,:);
NANG = length(ang);

obs1 = obs;
obs1(find(Rex == 0)) = NaN;

figure();
for ibnd=1:NBAND;
    subplot(3,3,ibnd);hold on;box on;
    for iang=1:NANG;
        if(Rex(ibnd,iang) == 1);
          plot(ang(iang),obs(ibnd,iang),'xk','LineWidth',2);
%        else;
%          alpha = .8;
%          plot(ang(iang),obs(ibnd,iang),'x','MarkerEdgeColor',[0,0,0]+alpha,...
%          'LineWidth',2);
        end;
    end;
    plot(ang,pre(ibnd,:),'-r','LineWidth',2);
    xlabel('Angle (deg.)');
    ylabel('Refl. Coeff.');
end;

%return;