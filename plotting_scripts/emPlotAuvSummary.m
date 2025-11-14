function [] = emPlotAuvSummary(baseFileName,depthMax,plotLogAlpha,trueFile,plotDataFit,showOnScreen,trueLogl,nfreq,iping)
close all
set(0, 'DefaultFigurePaperPosition', [0 0 16 12]);
mycolormap = [ones(1,3) ; jet(128)];

mPixels = 1000;
nPixels = 300;

%track_env = dlmread('../profile/track_environment.dat');
track_env = dlmread(trueFile);
isave = 1;
A = 0;

%% subaxis package values
%figw = 18;
%figh = 8;
%nx = 4;
%ny = 4;
ML = .05;
MR = .03;
%MB = .08;
MT = .02;
SP = .02;
PAD = 0;
%FNT = 12;
%inter = 20;
%%

zSample = zeros(1,250);

freqStr = ['1000 Hz';'1200 Hz';'2000 Hz';'2400 Hz'];

grd = linspace(0,depthMax,mPixels);
grd = (grd(2)-grd(1)/2)+grd;
kmx = 30;



clear ztru ctru atru rtru;
ktru = track_env(iping,1);
ztru = (track_env(iping,2:4:4*ktru+1));
ctru = [track_env(iping,3:4:4*ktru+2) track_env(iping,4*ktru+2)];
rtru = [track_env(iping,4:4:4*ktru+3) track_env(iping,4*ktru+3)];
atru = [track_env(iping,5:4:4*ktru+4) track_env(iping,4*ktru+4)];

%sampleFile = strrep(baseFileName,'.txt','_sample.mat');
sampleFile = strrep(baseFileName,'.txt','_particles.mat');
datTmp = dlmread(baseFileName);
saveFile = strrep(baseFileName,'.txt','_summary');
profileFile = strrep(baseFileName,'.txt','_profile');


%
%     paperUnits = get(fig1, 'PaperUnits');
% set(fig1,'PaperUnits','inches');
% paperSize = get(fig1,'PaperSize');
% paperPosition = [.5 .5 paperSize - .5];
% set(fig1,'PaperPosition', paperPosition);
% set(fig1,'PaperUnits',paperUnits);

%% LOAD SAMPLE
load(sampleFile)
disp (['working on ' sampleFile])

logL = A(:,1);

%% FIGURE 1
if(showOnScreen == true)
    fig1 = figure('visible','on');
else
    fig1 = figure('visible','off');
end
%% PLOT LOGL CHAIN
%subplot(4,4,1:3);hold on; box off;

subaxis(1,2,1,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);hold on;
set(gca,'FontSize',14);
plot(logL,'k-');
ylabel('log(L)');
xlabel('sample index');
if(exist('trueLogl','var')==1)
    plot([0,length(logL)],[trueLogl trueLogl],'w-','Linewidth',3);
    plot([0,length(logL)],[trueLogl trueLogl],'k--','Linewidth',2);
end;
ylim([250 450]);
hold off;

k = A(:,4);

%% PLOT K HISTOGRAM
%subplot(4,4,4);hold on;box on;

subaxis(1,2,2,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);hold on;
set(gca,'FontSize',14);
[n,lim]=hist(k,0:kmx);n = [0, n, 0];lim = [lim(1) lim lim(end)+1];
n = n/sum(n);
lim = lim-0.5;

[xx,yy]=stairs(lim,n,'k');
patch(xx,yy,[0.8,0.8,0.8]);
plot([ktru,ktru],[0,1],'w-','Linewidth',3);
plot([ktru,ktru],[0,1],'k--','Linewidth',2);

stairs(lim,n,'k');
clear n lim;
xlabel('No. interfaces in partition');
ylabel('Probability');
set(gca,'XLim',[0 kmx+0.5],'YLim',[0. 1.0]);


%% FIGURE 2
%% COMPUT MARGINALS
cSample=zeros(mPixels,size(A,1));
rSample=zeros(mPixels,size(A,1));
aSample=zeros(mPixels,size(A,1));
for iParticles = 1:size(A,1)
    zTmp = A(iParticles,5:4:4*(A(iParticles,4)+1));
    cTmp = A(iParticles,6:4:4*(A(iParticles,4)+1));
    rTmp = A(iParticles,7:4:4*(A(iParticles,4)+1));
    aTmp = A(iParticles,8:4:4*(A(iParticles,4)+1));
    chs = A(iParticles,4+4*A(iParticles,4)+1);
    rhs = A(iParticles,4+4*A(iParticles,4)+2);
    ahs = A(iParticles,4+4*A(iParticles,4)+3);
    
    
    mTmp = ones(mPixels,3);
    mTmp(:,1) = mTmp(:,1).*chs;
    mTmp(:,2) = mTmp(:,2).*rhs;
    mTmp(:,3) = mTmp(:,3).*ahs;
    for ik = k(iParticles):-1:1
        index = find(grd<zTmp(ik));
        mTmp(index,1) = cTmp(ik);
        mTmp(index,2) = rTmp(ik);
        mTmp(index,3) = aTmp(ik);
    end
    
    zSample = zSample+histc(zTmp,0:depthMax/249:depthMax);
    cSample(:,iParticles) = mTmp(:,1);
    rSample(:,iParticles) = mTmp(:,2);
    aSample(:,iParticles) = mTmp(:,3);
end


if(showOnScreen == true)
    fig2 = figure('visible','on');
else
    fig2 = figure('visible','off');
end
%% velocity plot
hc = histc(cSample,linspace(1450,1750,nPixels),2);
%subplot(4,4,[6 10]);set(gca,'FontSize',14);hold on;
subaxis(1,4,2,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);
set(gca,'FontSize',14);hold on;
imagesc(linspace(1450,1750,nPixels),0:depthMax/999:depthMax,hc);set(gca,'Ydir','reverse');
stairs([ctru(1) ctru],[0 ztru depthMax],'w-', 'LineWidth',3);stairs([ctru(1) ctru],[0 ztru depthMax],'k--', 'LineWidth',2);%stairs(ctru,[0 ztru],'w--')
xlabel('Velocity (m/s)');
%set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
set(gca,'YTickLabel',[],'TickDir','out','Layer','top');box on;
xlim([1450 1750]);
ylim([0 depthMax]);
box on;hold off;

hc = histc(rSample,linspace(1.2,2.2,nPixels),2);


%[hc,bins] = hist(rSample,binr);

%% density plot
%subplot(4,4,[7 11]);
subaxis(1,4,3,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);
set(gca,'FontSize',14);hold on;
imagesc(linspace(1.2,2.2,nPixels),0:depthMax/999:depthMax,hc);set(gca,'Ydir','reverse');
stairs([rtru(1) rtru],[0 ztru depthMax],'w-', 'LineWidth',3);stairs([rtru(1) rtru],[0 ztru depthMax],'k--', 'LineWidth',2);
xlabel('Density (g/ccm)');
set(gca,'YTickLabel',[],'TickDir','out','Layer','top');box on;
%set(gca,'XTickLabel',[1.4 1.6 1.8 2.0]);
set(gca,'XTick',1.:.2:3.0);
xlim([1.0 2.2]);
ylim([0 depthMax]);
colormap(mycolormap);
hold off;

%% attenuation plot
if(plotLogAlpha)
    atru = log10(atru); %aSample = 10.^aSample;
    hc = histc(aSample,linspace(-3,0,nPixels),2);
    %subplot(4,4,[8 12]);
    subaxis(1,4,4,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);
    set(gca,'FontSize',14);hold on;
    
    imagesc(linspace(-3,0,nPixels),0:depthMax/999:depthMax,hc);set(gca,'Ydir','reverse');
    
    stairs([atru(1) atru],[0 ztru depthMax],'w-', 'LineWidth',3);stairs([atru(1) atru],[0 ztru depthMax],'k--', 'LineWidth',2);
    xlabel('Attenuation (dB/m/kHz)');
    set(gca,'YTickLabel',[],'TickDir','out','Layer','top');
    xlim([-3 0]);
    ylim([0 depthMax]);
    colormap(mycolormap);
    box on;hold off;
else
    hc = histc(aSample,linspace(0.001,1,nPixels),2);
    %subplot(4,4,[8 12]);
    subaxis(1,4,4,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);
    set(gca,'FontSize',14);hold on;
    
    imagesc(linspace(0.001,1,nPixels),0:depthMax/999:depthMax,hc);set(gca,'Ydir','reverse');
    
    stairs([atru(1) atru],[0 ztru depthMax],'w-', 'LineWidth',3);stairs([atru(1) atru],[0 ztru depthMax],'k--', 'LineWidth',2);
    xlabel('Attenuation (dB/m/kHz)');
    set(gca,'YTickLabel',[],'TickDir','out','Layer','top');
    xlim([0.001 1]);
    ylim([0 depthMax]);
    colormap(mycolormap);
    box on;hold off;
end

%% z depth histogram
%subplot(4,4,[5 9]);set(gca,'FontSize',14);hold on;

subaxis(1,4,1,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);
set(gca,'FontSize',14);hold on;
yaxisDisc = 0:depthMax/249:depthMax;
halfBin = depthMax/249*0.5;
yaxisDisc = yaxisDisc+2*halfBin;

%[xx,yy] = stairs(zSample/sum(zSample),0:depthMax/249:depthMax);
[xx,yy] = stairs(zSample/sum(zSample),yaxisDisc);
%sum(zSample/sum(zSample))
patch([0;xx],[0;yy],[0.8,0.8,0.8]);set(gca,'Ydir','reverse');
xlim([0 0.1]);
ylim([0 depthMax]);
ylabel('Depth(m)');
xlabel('Probability');
%set(gca,'YTickLabel',0:0.5:32);
%set(gca,'YTickLabel',[0 0.1]);
for jtmp = 1:length(ztru)
    plot([0,1],[ztru(1,jtmp),ztru(1,jtmp)],'w-', 'LineWidth',3);plot([0,1],[ztru(1,jtmp),ztru(1,jtmp)],'k--', 'LineWidth',2);%    plot([0,1],[ztru(ifile,1),ztru(ifile,1)],'w--');
end
colormap(mycolormap);
set(gca,'TickDir','out','Layer','top');
box on;hold off;

if (plotDataFit == true)
    %% FIGURE 3
    m = ceil(nfreq/3);
    if(showOnScreen == true)
        fig3 = figure('visible','on');
    else
        fig3 = figure('visible','off');
    end
    % PLOT DATA FIT
    rep = load(['replica_' strrep(baseFileName,'.txt','_particles.txt')]);
    n_size = size(rep,1) - 1 - nfreq;
    for i = 1:nfreq
        subaxis(m,3,i,'Spacing',SP,'Padding',PAD,'ML',ML,'MR',MR,'MT',MT);
        set(gca,'FontSize',14);
        plot(rep(n_size+1,:),rep(i:nfreq:n_size,:),':','Color',[0.8,0.8,0.8]);hold on;
        plot(rep(n_size+1,:),rep(n_size+1+i,:),'k+');hold off;
        ylim([0 1]);
        %xlim([0 90]);
        xlim([29 67]);
        if(mod(i-1,3) == 0)
            set(gca,'YTick',0.2:0.2:1,'TickDir','out','Layer','top');
            ylabel('R');
            box on;
        else
            set(gca,'YTickLabel',[]);
            set(gca,'YTick',0.2:0.2:1,'TickDir','out','Layer','top');box on;
        end
        if(i == nfreq || i == nfreq -1 || i == nfreq - 2)
            set(gca,'XTick',30:30:90,'TickDir','out','Layer','top');
            %set(gca,'XTick',10:5:25,'TickDir','out','Layer','top');
            xlabel('angle (\circ)');
            box on;
        else
            set(gca,'XTick',30:30:90,'TickDir','out','Layer','top');
            
            %set(gca,'XTick',15:5:25,'TickDir','out','Layer','top');
            set(gca,'XTickLabel',[]);
        end
    end
    
end

if (isave == 1)
    disp(['saving file ' saveFile '.png'])
    saveas(fig1,saveFile,'png');
    saveas(fig2,profileFile,'png');
    if(plotDataFit)
    saveas(fig3,strrep(baseFileName,'.txt','_datafit'),'png');
    end
end
clear datTmp repTmp