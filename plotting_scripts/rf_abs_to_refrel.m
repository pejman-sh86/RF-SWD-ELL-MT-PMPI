function [] = rf_plot_hist_varpar(filename);

filebase    = strrep(filename,'sample.mat','');
velreffile  = strrep(filename,'_sample.mat','_vel_ref.txt');
mapfile     = strrep(filename,'_sample.mat','_map_voro.dat');
mapfiletru  = strrep(filename,'_sample.mat','_map_voro_true.dat');
parfile     = strrep(filename,'_sample.mat','_parameter.dat');

NPL = 3;
[IMAP,ICOV,iar,i_varpar,irv,itr,iswd,ivref,ivpvs,ISMPPRIOR,IEXCHANGE,idip,...
 NDAT_SWD,NMODE,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,NKEEP,NPTCHAINS1,...
 dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,sampling_dt,dVs,dVpVs,...
 ntr,baz]=rf_read_parfile(parfile);
ktru = 4;
voroabs(1,:) = [  0.0, 3.5, 1.80 ];
voroabs(2,:) = [ 10.0, 3.7, 1.80 ];
voroabs(3,:) = [ 37.0, 4.7, 1.75 ];
voroabs(4,:) = [180.0, 4.3, 1.85 ];

map = load(mapfile);
map(1)=ktru;
map(2:NVMX*NPL) = 0;

%%
%% Ref velocity model:
%%
tmp = dlmread(velreffile);
NREF = tmp(1,1);
ntmp = tmp(1,2);
NPREM = ntmp-NREF;
vel_ref = tmp(2:NREF+1,:);
vel_prem = tmp(NREF+1:NREF+1+NPREM,:);

voropert = voroabs;
for ivo = 1:ktru;
  [vref,vpvsref]=rf_getref(voroabs(ivo,1),vel_ref);
  voropert(ivo,2) = - vref + voroabs(ivo,2);
  voropert(ivo,3) = - vpvsref + voroabs(ivo,3);
end;

for ivo = 1:ktru;
  map(2+(ivo-1)*NPL:2+(ivo*NPL)-1) = squeeze(voropert(ivo,:));
end;

save(mapfiletru,'map','-ascii');

return;