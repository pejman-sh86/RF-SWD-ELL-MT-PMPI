function [] = convert_replicas(istart);

%% Malta Plateau track:
NBAND = 6;
files=dir('*_rep_ens.dat');
%% Simulation track:
%NBAND = 4;
%files=dir('*_replicas.txt');
NFILE = length(files)
for ifile = istart:NFILE;
  disp(['converting  ',files(ifile).name]);
  repfile = files(ifile).name;
  reptmp=dlmread(repfile);
  %reptmp(end,:)=[];
  NANG = size(reptmp,2);
  NDAVE = length(reptmp)/NBAND;
  ref = zeros(NANG,NBAND,NDAVE);
  for j=1:NDAVE;
    for iband=1:NBAND;
      ref(:,iband,j) = reptmp(NBAND*(j-1)+iband,:);
    end;
  end;
  ref = ref(:,:,1:2:end);
  NDAVE = size(ref,3);
  
  if(repfile(end-2:end) == 'dat');
    matfile = strrep(repfile,'.dat','.mat');
  else;
    matfile = strrep(repfile,'.txt','.mat');
  end;
  matfile
  save(matfile,'ref');

end;

return;
