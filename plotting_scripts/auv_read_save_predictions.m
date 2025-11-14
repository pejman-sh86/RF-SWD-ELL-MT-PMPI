%function auv_read_save_predictions;

istart = 1;
files=dir('*_particles.mat');
%files=dir('*2700_sample.mat');

pinglist = load('ping_list.txt');
NPING = pinglist(1);
%NPING = 5;
pinglist(1)=[];
ipst = pinglist(1);

%% Malta 
bands = [988. 1113. 1288. 1913. 2263. 2513.];

NBAND = length(bands);
NANG = 32;
NDAVE = 1200;
NSD = NBAND;
ifile2 = istart-1;

ref2 = zeros(NANG,NBAND,NDAVE,NPING);
dobs = zeros(NANG,NBAND,NPING);
rex = zeros(NANG,NBAND,NPING);
rex = 1;
for ifile=1:NPING;
    
  ifile2 = ifile2 + 1;
  filename = files(ifile2).name;
  if(rem(ifile,100)==0);disp(filename);end;
  filebase    = strrep(filename,'_particles.mat','');
  repfile     = strcat(filebase,'_rep_ens.mat');
  datafile     = strcat(filebase,'.txt');

  load(repfile);
  ntmp1 = size(ref,1);
  ntmp2 = size(ref,2);
  ntmp3 = size(ref,3);
  ref2(1:ntmp1,1:ntmp2,1:ntmp3,ifile2) = ref;
  
  tmp = dlmread(datafile);
  z_t    = tmp(1,1);
  cw     = tmp(2,1);
  rw     = tmp(3,1);
  hmax   = tmp(4,1)+tmp(4,1)/10.;
  hmax2  = hmax;
  dobs(:,:,ifile2)   = tmp(5:length(bands)+4,:)';
  angobs = tmp(length(bands)+5,:);
  %rex(:,:,ifile2)   = tmp(6+length(bands):length(bands)+length(bands)+5,:)';
end;
save('ref2.mat','ref2','dobs','rex','angobs','-v7.3');
%return;