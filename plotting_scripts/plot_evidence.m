%----------------------------------------------------------
% Compute Evidence from PPD
%----------------------------------------------------------
function [] = plot_evidence(efile);

load(efile);
i_save = 1;
plotfile1 = strrep(efile,'.mat','1.eps');
plotfile2 = strrep(efile,'.mat','2.eps');
plotfile3 = strrep(efile,'.mat','1.fig');
plotfile4 = strrep(efile,'.mat','2.fig');

efile
plotfile1
plotfile2

log10EN = log10E - max(log10E);
occ = occ/log(10);
occn    = occ - max(occ);
BICN    = BIC - min(BIC);
AICN    = AIC - min(AIC);

%log10EN2 = log10E - (max(log10E) - (max(log10E)-min(log10E))/2);
%EN2 = 10.^log10EN2;
%EN2N = EN2/sum(EN2)
%occn2    = occ - (max(occ) - (max(occ)-min(occ))/2);
%BICN2    = BIC - (min(BIC) - (max(BIC)-min(BIC))/2);
%AICN2    = AIC - (min(AIC) - (max(AIC)-min(AIC))/2);

%%
%% Sim B:
%%
%x1 = [3];
%x2 = [7];
%x3 = [3];
%x4 = [7];
%col = {'-*b' '-*k' '*r'};
%col = char(col);
%xticks = [3,4,5,6,7];


%%
%% Sim C:
%%
%x1 = [1  1  5];
%x2 = [7  6  5];
%x3 = [1  8 14];
%x4 = [7 13 14];
%col = {'-*b' '-*k' '*r'};
%col = char(col);
%leg = {'no grad','grad'};
%xticks = [1,2,3,4,5,6,7];

%%
%% Sim C known sd:
%%
%x1 = [1  1 ];
%x2 = [6  6 ];
%x3 = [1  7 ];
%x4 = [6 12 ];
%col = {'-*b' '-*k' '*r'};
%col = char(col);
%leg = {'no grad','grad'};
%xticks = [1,2,3,4,5,6];

%%
%% Site 01:
%%
x1 = [0];
x2 = [3];
x3 = [1];
x4 = [4];
col = {'-*b'};
col = char(col);
leg = {'no grad'};
xticks = [0,1,2,3];

%%
%% Site 02:
%%
%x1 = [3];
%x2 = [7];
%x3 = [4];
%x4 = [8];
%col = {'-*b' '-*b'};
%col = char(col);
%leg = {'no grad'};
%xticks = [3,4,5,6,7];

%%
%% Site 02 brute:
%%
%x1 = [0];
%x2 = [5];
%x3 = [1];
%x4 = [6];
%col = {'-*b'};
%col = char(col);
%leg = {''};
%xticks = [0,1,2,3,4,5];

%%
%% Site 13:
%%
%x1 = [ 4];
%x2 = [ 8];
%x3 = [10];
%x4 = [14];
%col = {'-*b' '--ob' '-*k' '--ok'};
%col = char(col);
%leg = {'case a','case b'};
%xticks = [4,5,6,7,8];

%%
%% Site 13:
%%
%x1 = [2  3];
%x2 = [8  7];
%x3 = [1  8];
%x4 = [7 12];
%col = {'-*b' '--ok'};
%col = char(col);
%leg = {'case a','case b'};
%xticks = [1,2,3,4,5,6];

%
% ---------------------------------------------------------------------------

nx = 1;
ny = 2;
nsubfig = nx*ny;
xim = 0.03;
yim = 0.1/ny;
xymarg = [0.14 0.04 0.04 0.10];
opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',12,'Height',8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);


gc1 = figure(1);
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;

%
% ---------------------------------------------------------------------------

subplot(h1)
hold on;
box on;

for i = 1:length(x1)
   plot([x1(i):x2(i)],log10EN( x3(i):x4(i)),col(i,:))
end;

ylabel('log10(E(H))','FontSize',14)
set(gca,'XLim',[min(x1)-1 max(x2)+1])
set(gca,'XTick',[xticks],'FontSize',14);
set(gca,'XTickLabel',[],'FontSize',14);
set(gca,'YGrid','on')
%legend(leg)
%set(gca,'YLim',[-530 -300])
%set(gca,'YLim',[-530 -300])

%
% ---------------------------------------------------------------------------

subplot(h2)
hold on;
box on;

for i = 1:length(x1)
   plot([x1(i):x2(i)],minL( x3(i):x4(i)),col(i,:))
end;

%for i = 1:length(x1)
%   plot([x1(i):x2(i)],L( x3(i):x4(i)),col(i,:))
%end;

ylabel('Misfit','FontSize',14)
xlabel('No. Layers','FontSize',14)
set(gca,'XLim',[min(x1)-1 max(x2)+1])
set(gca,'XTick',[xticks],'FontSize',14);
set(gca,'YGrid','on')
%legend(leg)

nx = 1;
ny = 3;
nsubfig = nx*ny;
xim = 0.03;
yim = 0.1/ny;
xymarg = [0.14 0.04 0.04 0.10];
opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',12,'Height',8,'Color','cmyk',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

gc2 = figure(2);
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
hold on; box off;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
hold on; box off;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
hold on; box off;

%
% ---------------------------------------------------------------------------

subplot(h1)
hold on;
box on;

for i = 1:length(x1)
   plot([x1(i):x2(i)],occn( x3(i):x4(i)),col(i,:))
end;

ylabel('Linear approx.','FontSize',14)
set(gca,'XLim',[min(x1)-1 max(x2)+1])
set(gca,'XTick',[xticks],'FontSize',14);
set(gca,'XTickLabel',[],'FontSize',14);
set(gca,'YGrid','on')
%legend(leg)

%
% ---------------------------------------------------------------------------

subplot(h2)
hold on;
box on;

for i = 1:length(x1)
   plot([x1(i):x2(i)],BICN( x3(i):x4(i)),col(i,:))
end;

ylabel('BIC','FontSize',14)
set(gca,'XLim',[min(x1)-1 max(x2)+1])
set(gca,'XTick',[xticks],'FontSize',14);
set(gca,'XTickLabel',[],'FontSize',14);
set(gca,'YGrid','on')
%legend(leg)

%
% ---------------------------------------------------------------------------

subplot(h3)
hold on;
box on;

for i = 1:length(x1)
   plot([x1(i):x2(i)],AICN( x3(i):x4(i)),col(i,:))
end;

ylabel('AIC','FontSize',14)
xlabel('No. Layers','FontSize',14)
set(gca,'XLim',[min(x1)-1 max(x2)+1])
set(gca,'XTick',[xticks],'FontSize',14);
set(gca,'YGrid','on')
%legend(leg)

if(i_save == 1)
   exportfig(gc1,plotfile1,opts);
   exportfig(gc2,plotfile2,opts);
%   saveas(gc1,plotfile3,'fig');
%   saveas(gc2,plotfile4,'fig');
end

return;
