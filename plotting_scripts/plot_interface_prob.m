
DSlabel = {'M', 'SM', 'RSM'};
deplim = [0. 15. 60. 220.];
nbins = '2500';

NDS = length(DSlabel);
binsdepA = [];
NdepA = [];

for iDS = 1:NDS
    load( strcat(strcat(DSlabel{iDS},'_hist_'),nbins) );
    binsdepA = [binsdepA binsdep'];
    NdepA = [NdepA Ndep'];
end
%%
fig=figure;
color = [1. 0. 0.;...
         0. 1. 0.;...
         0. 0. 0.];
Pmax = [.8 .25 .065];
font_size = 25;
figw = 12;
figh = 6;
set(fig,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
ny = length(deplim) - 1;
nx = 1;
xim = 0.01;
yim = 0.05;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
% spw = spw - 0.005;
% spw2 = spw;
spw = spw / 4.;
% XTick = [0:5:10];
%%
for ipanel = 1:ny

    h1 = subplot('Position',[loc(1,ipanel) loc(2,ipanel) spw sph]);
    hold on; box on;
    for ids = 1:NDS
        stairs(binsdepA(:,ids),NdepA(:,ids),'Color',color(ids,:));
%     [xx,yy]=stairs(binsdepA(:,ids),NdepA(:,ids));
%     patch(xx,yy,color(ipa:));
    end

    set(gca,'Fontsize',font_size)
%     set(gca,'XLim',[deplim(1) deplim(ipanel+1)])
    set(gca,'XLim',[deplim(ipanel) deplim(ipanel+1)])
%     set(gca,'YLim',[0., Pmax(ipanel)]);
    set(gca,'TickDir','in');
    ylabel({'Interface', 'probability'});
    xlabel('Depth (km)');
    if (ipanel==1); legend('MT','SWD-MT','RF-SWD-MT','FontSize',16); end

end