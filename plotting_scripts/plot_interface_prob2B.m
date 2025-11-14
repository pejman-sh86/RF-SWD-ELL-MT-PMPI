
% DSlabel = {'M', 'S', 'R', 'SM', 'RS', 'RSM'};
DSlabel = {'M', 'SM', 'RSM'};
% deplim = { [25. 50.], [90. 110.], [120. 250.] }; %sim7
% deplim = { [0. 15.], [15. 60,], [60. 220.] }; %real
% deplim = { [0. 5.], [5. 15.], [15. 60,], [60. 110.], [110. 220.]}; %real
% deplim = { [0. 5.], [5. 15.], [15. 60,], [110. 220.]}; %real
% deplim = { [0. 5.], [5. 15.], [15. 60,], [110. 250.]}; %sim8
% deplim = { [0. 70.], [70. 140.], [140. 250.] }; %sim9
deplim = { [0. 5.], [5. 15.], [15. 60,], [60. 250.]}; %sim10
% nbins = {'400', '400', '100'}; %sim7
% nbins = {'2000', '400', '100'}; %real
% nbins = {'2000', '2000', '400', '100', '100'}; %real
% nbins = {'2000', '2000', '400', '100'}; %real
nbins = {'2000', '2000', '400', '100'}; %real, sim8, sim10
% nbins = {'400', '400', '100'}; %sim9
Npanel = length(deplim);

itrue = 1;
% dep_true = [40., 100., 200.]; %sim7
dep_true = [1., 8., 35., 170.]; %sim8, sim10
% dep_true = [40., 100., 200.]; %sim9

fig=figure;
color = [1. 0. 0.;...
         0. 1. 0.;...
         0. 0. 0.];
% color = [1. 0. 0.;...
%          0. 0. 1.;...
%          .5 .5 0.;...
%          0. 1. 0.;...
%          1. 1. 0.;...
%          0. 0. 0.];
% Pmax = [.4 .2 .03]; %sim7
% Pmax = [1.5 .6 .06]; %sim7
% Pmax = [4. 1.6 .9 .05]; %sim8
% Pmax = [.65 .2 .03]; %sim9
Pmax = [6. 2. 1. .05]; %real_3cases, real_6cases
%Pmax = [.65 .105 .03]; %real_old
% Pmax = [5. 4. .6 .08 .08]; %real
% Pmax = [4.5 4. .6 .08]; %real_3cases, real_6cases
font_size = 25;
nx = 3;
% nx=Npanel/2;
ny=2; 
% nx=1; ny=length(deplim); Npanel=ny;
xim = 0.05; yim = 0.055;%sim8
% xim = 0.05; yim = 0.01;%sim7
% xim = 0.01; yim = 0.05;
% xim = 0.05; yim = .055; %real_3cases, real_6cases
xymarg = [0.1 0.04 0.04 0.1];
[loc2,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
% sph = sph - 0.005;
% spw = spw - 0.005;
% spw2 = spw;
% spw = spw*9./10.;
% for ipanel=1:Npanel
%     loc2(1,ipanel) = loc2(1,ipanel) - (ipanel-1)*(spw2-spw);
% end

% loc = loc2; %sim7, sim9
loc = zeros(2, Npanel); %sim8, sim10, and real
jump = nx - Npanel/2;
loc(:,1:Npanel/2) = loc2(:,1:Npanel/2);
loc(:,Npanel/2+1:Npanel) = loc2(:,Npanel/2+1+jump:Npanel+jump);

% XTick = [0:5:10];

%%
NDS = length(DSlabel);
binsdepA = cell(1, Npanel);
NdepA = cell(1, Npanel);

for idepth = 1:Npanel
    for iDS = 1:NDS
        load( strcat(DSlabel{iDS},'_hist_',nbins{idepth}) );
        binsdepA{idepth} = [binsdepA{idepth} binsdep'];
        NdepA{idepth} = [NdepA{idepth} Ndep'];
    end
end

%%
for ipanel = 1:Npanel

    binsdepA2 = binsdepA{ipanel}; NdepA2 = NdepA{ipanel};
    %
    deplimA2 = deplim{ipanel};
    %
    h1 = subplot('Position',[loc(1,ipanel) loc(2,ipanel) spw sph]);
    hold on; box on;
    for ids = 1:NDS
        % normalize distributions in each panel
        AA = binsdepA2(binsdepA2(:,ids)>=deplimA2(1)&binsdepA2(:,ids)<=deplimA2(2), ids);
        BB = NdepA2(binsdepA2(:,ids)>=deplimA2(1)&binsdepA2(:,ids)<=deplimA2(2), ids);
        AA = [deplimA2(1); AA; deplimA2(2)];
        BB = [BB(1); BB; BB(end)];
        BB = BB / trapz(AA,BB);
        if (ids==1); BBMTmin = min(BB);end
%         if (ipanel==1); BB = BB - 0.03; end %sim7
%         BB = BB - BBMTmin;
%         if (ipanel~=1); BB = BB - BBMTmin; end 
%         if(ipanel~=5); BB = BB - BBMTmin; end
%         if(ipanel~=3); BB = BB - BBMTmin; end
        %
%         stairs(AA,BB); 
%         if (ids==1 || ids==4 || ids==6); stairs(AA,BB,'Color',color(ids,:));else;stairs(AA,BB);end 
        stairs(AA,BB,'Color',color(ids,:)); 
%         stairs(binsdepA2(:,ids),NdepA2(:,ids),'Color',color(ids,:));
%     [xx,yy]=stairs(binsdepA2(:,ids),NdepA2(:,ids));
%     patch(xx,yy,color(ipa:));
    end

    set(gca,'Fontsize',font_size)
%     set(gca,'XLim',[deplim(1) deplim(ipanel+1)])
%     set(gca,'XLim',[deplim(ipanel) deplim(ipanel+1)])
    set(gca,'XLim',deplim{ipanel})
    set(gca,'YLim',[0., Pmax(ipanel)]);
    set(gca,'TickDir','in');
%     if (ipanel==1); ylabel({'Interface', 'probability'}); end
%     if (ipanel==1); ylabel('Interface probability'); end
    if (ipanel==1 || ipanel==3); ylabel('Interface probability'); end
%     ylabel({'Interface', 'probability'});
%     xlabel('Depth (km)');
    if (ipanel==3 || ipanel==4); xlabel('Depth (km)'); end
%     if (ipanel==1); legend('MT','SWD-MT','RF-SWD-MT','FontSize',16); end %real_3cases
%     if (ipanel==1); legend('MT', 'SWD', 'RF', 'MT-SWD', 'SWD-RF', 'MT-SWD-RF','FontSize',16); end %real_6cases
%     if (ipanel==4); set(gca,'XTick', [120.:40.:220.]); end %real_3cases, real_6cases
    if (itrue)
        ylim = get(gca,'Ylim');
%         plot([dep_true(ipanel), dep_true(ipanel)],[ylim(1), ylim(2)],'--y','LineWidth',.5)
        plot([dep_true(ipanel), dep_true(ipanel)],[ylim(1), ylim(2)],'--k','LineWidth',2.)
        if (ipanel==1); legend('MT', 'MT-SWD', 'MT-SWD-RF','True depth','FontSize',16); end
%         if (ipanel==1); legend('MT', 'SWD', 'RF', 'MT-SWD', 'SWD-RF', 'MT-SWD-RF','True depth','FontSize',16); end
    end

end