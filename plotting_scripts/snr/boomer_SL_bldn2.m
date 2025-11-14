% boomer_sl_blde4; SS98 Site 4
% Extract Boomer SL using received level and TL model
% plot results of piston model
% INPUTS: TL modeled and RL measured
% eliminates outliers more than 3 standard deviations

%select third octave band numbers for plotting
fr_otonum=[22:39]; szfr=(size(fr_otonum)); numfrq=szfr(2);
snr=6; snrx=100;
cs=1514.44 ;  cb=1510.8;
dirc=''
id='id4_neg'; phone='2'; poly_num=6;
%id='id6'; phone='2'; poly_num=6;

%maxrange=1000; minrange=0; colr2='b';  
maxrange=0; minrange=-800; colr2='r';

load otofc.dat;  load otofedg.dat;
frq=otofc(fr_otonum)
foff=4 - 20; frcl=foff+fr_otonum;

dirnm='d:\experiment_data\Sicily Scarab98\bl_vla\Site2\';
dirnmtl=[dirnm 'Proc Deep Boomer\no surf\'];

ymin=30; ymax=100; iplot=1; ylab=1:3:24;
 
% LOAD DIRECT PATH RL DATA;
fnbl=['bldn2_' dirc phone]
rld=load ([ dirnm 'bl data\hp_results\' fnbl '_rld_' id '.asc']);
tldsnr=load ([ dirnm 'bl data\hp_results\' fnbl '_rldsnr_' id '.asc']);
sztld=size(rld);
sld=-999*ones(size(rld)); sld(:,1:3)=rld(:,1:3);
rldscr=-999*ones(size(rld)); rldscr(:,1:3)=rld(:,1:3);

boomer_SL

%  save(fsnm, 'sld', '-ascii')
