function [] = resamp(filename);

filebase = strrep(filename,'.mat','');

i_inter = 0;
nang_re = 31;
load(filename);
x.ang=fliplr(x.ang);
x.Rsmooth=fliplr(x.Rsmooth);

%% Frequencies:
%%1012 1212 2362
idx = [8, 16, 64];
%idx2= [9, 17, 63];
%%1000 1200 2000 2400 2800 3200
%idx = [8, 16, 48, 64, 80, 96];
%%1000 1200 2000 2200 2400 2600 2800 3000 3200
%idx = [8, 16, 48, 56, 64, 72, 80, 88, 96];
%% Frequencies:
%%1000 1200 2000 2200 2400 2600 2800 3000 3200
%idx = [48, 56, 64, 72, 80, 88, 96];

nf = length(idx);
theta = x.ang*360./(2.*pi);
R = x.Rsmooth(idx,:);
%R = (x.R(idx,:)+x.R(idx2,:))/2.;
z_t = 2.*x.altitude;
cw = x.c;
rw = 1.029;
hmx = x.tb_win*1600./2.;

outfile  = strcat(filebase,'_',num2str(x.freq(idx(1))),'_',num2str(x.freq(idx(end))),'.mat')
outfilea  = strcat(filebase,'_',num2str(x.freq(idx(1))),'_',num2str(x.freq(idx(end))),'.txt')

thdiff = diff(theta);
if(i_inter == 1);
   theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];
else
   theta_even = theta;
end;

theta_even
Rex = ones(nf,length(theta_even));
for ifr = 1:nf

   R_even(ifr,:)   = interp1(theta,R(ifr,:),theta_even);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = cw
F(1).rw = rw;

F(1).freq = x.freq(idx)
for ifr = 1:nf;
    F(ifr).dat = R_even(ifr,:);
    F(ifr).ang = theta_even;
end
F(1).Rex = Rex;


save(outfile,'F');

filebase = strrep(outfile,'.mat','');
fid = fopen('filebase.txt','w');
filebase
fprintf(fid,'%i\n',length(filebase));
fprintf(fid,'%s',filebase);
fclose(fid);

save(outfilea,'z_t','-ascii');
save(outfilea,'cw','-ascii','-append');
save(outfilea,'rw','-ascii','-append');
save(outfilea,'hmx','-ascii','-append');
for i = 1:nf

   tmp = F(i).dat;
   save(outfilea,'tmp','-ascii','-append');

end
save(outfilea,'theta_even','-ascii','-append');

for i = 1:nf

    tmp = F(1).Rex(i,:);
    save(outfilea,'tmp','-ascii','-append');

end

return;
