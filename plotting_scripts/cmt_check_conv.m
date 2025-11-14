%function []=cmt_check_conv(filename);

%filename = 'maule_data_sample.mat';
%filename = 'bo10_data_sample.mat';
%filename = 'ec10_data_sample.mat';
filename = 'ma10_data_sample.mat';
%filename = 'mi10_data_sample.mat';
%filename = 'am16_data_sample.mat';
%filename = 'cr17_data_sample.mat';
%filename = 'hi17_data_sample.mat';
%filename = 'kk17_data_sample.mat';
%filename = 'td17_data_sample.mat';
%filename = 'ak17_6.3_sample.mat';
%filename = 'ak17_5.7_sample.mat';
%filename = 'an17_data_sample.mat';

isave   = 1;
ieps    = 1;
idatfit = 0;
isyn    = 0;
inonstat= 1;
FNT = 21;

%% Was sampling done with double log prior on M0?
i_dlogprior = 1;

filebase = strrep(filename,'_sample.mat','');
parfile  = strcat(filebase,'_parameter.dat');
M0file   = strcat(filebase,'_M0.dat');
datfile  = strcat(filebase,'.hdf5');
[IMAP,ICOV,NVMX,ILOC,IAR,IEXCHANGE,IDBLCPL,...
 NPTCHAINS1,dTlog,ICHAINTHIN,NKEEP,IADAPT,NBUF,...
 minlat,maxlat,minlon,maxlon,mindepth,maxdepth,...
 mindelay,maxdelay,minMw,minstr,mindip,minrak,...
 maxMw,maxstr,maxdip,maxrak]=cmt_read_parfile(parfile);
[NMRF,NMETA,NSTN,NDAT,dobs,NTSMP]=cmt_read_datafiles();

load(filename);
NSMP = length(A);
N = round(NSMP/3);
A1 = A(1:N,:);
A3 = A(end-N:end,:);

figure();
subplot(1,2,1);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[n1,lim]=hist(A1(:,1),80);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
n1 = n1/trapz(lim,n1);
[xx,yy]=stairs(lim,n1,'b');
patch(xx,yy,'b');
stairs(lim,n1,'b');
clear n1 lim;
ylabel('Probability');

set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
[n1,lim]=hist(A3(:,1),80);n1 = [0, n1, 0];lim = [lim(1) lim lim(end)];
n1 = n1/trapz(lim,n1);
[xx,yy]=stairs(lim,n1,'r');
patch(xx,yy,'r');
stairs(lim,n1,'r');
clear n1 lim;
xlabel('Probability');
alpha(.3)

subplot(1,2,2);hold on;box on;
set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
stairs(A1(:,1),'b');
clear n1 lim;
xlabel('Probability');

set(gca,'FontSize',FNT,'layer','top','LineWidth',1)
stairs(A3(:,1),'r');
clear n1 lim;
xlabel('Probability');
alpha(.3)











