function [] = auv_compute_buckingham_disp_track(istart);
%clear all;
close all;
%istart = 1;
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

imarg   = 1; %% Plot depth-marginal distributions?
inorm   = 1; %% Normalize profile marginals line by line to unit area
isyn    = 0;
i_vref  = 1;
ibucking= 3;
ivgsnosh= 0;
ishear  = 0;
ivoro   = 0;
imap    = 0;    
imead   = 0;
imean   = 0;
imax    = 0;      %%
isave   = 1;      %%
idatfit = 0;      %% plot data misfits?
isd     = 1;      %% sampled over sigma?
iar     = 1;
ivpvs   = 0;
icore   = 0;
itrackmat = 0;    %% 1 -> Already pre-computes track.mat
FNT = 14;

%files=dir('*_particles.mat');
files=dir('*2513_sample.mat');

pinglist = load('ping_list.txt');
NPING = pinglist(1);
dping = 300;      %% Interval for ping axis label
pinglist(1)=[];

NPING = length(files);
filename = files(istart).name
[fr, z, nfc, nfa, nfr, meac, mear, meaa]=compute_buckingham_disp(filename);
[Nfr,NZ,tmp] = size(nfc);
nfct = zeros(NPING,Nfr,NZ,2);
nfrt = zeros(NPING,Nfr,NZ,2);
nfat = zeros(NPING,Nfr,NZ,2);
meact = zeros(NPING,Nfr,NZ);
meart = zeros(NPING,Nfr,NZ);
meaat = zeros(NPING,Nfr,NZ);
nfct(1,:,:,:)  = nfc;
nfrt(1,:,:,:)  = nfr;
nfat(1,:,:,:)  = nfa;
meact(1,:,:) = meac;
meart(1,:,:) = mear;
meaat(1,:,:) = meaa;

if(itrackmat == 0)
  ifile2 = istart;
  for ifile=1:NPING;

       ifile2 = ifile2 + 1;
       if(ifile2 > pinglist(end));exit;end;
       filename = files(ifile2).name
       [fr, z, nfc, nfa, nfr, meac, mear, meaa]=compute_buckingham_disp(filename);
       
       nfct(ifile,:,:,:)  = nfc;
       nfrt(ifile,:,:,:)  = nfr;
       nfat(ifile,:,:,:)  = nfa;
       meact(ifile,:,:) = meac;
       meart(ifile,:,:) = mear;
       meaat(ifile,:,:) = meaa;
  end;
  save(strcat('track_cra',num2str(istart),'.mat'),'-v7.3');
  ifile2 = istart - 1;
else
  ifile2tmp = istart-1;
  load track_cra.mat
  ifile2 = ifile2tmp;
  hmax2 = 7.5;
end;
