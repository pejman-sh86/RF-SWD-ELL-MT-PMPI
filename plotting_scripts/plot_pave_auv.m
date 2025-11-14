
idx = [5, 6, 7];
ipskip = 5;
filebase1 = 'p';
filebase2 = '_pave_10_1000_2400';
fileext  = '.txt';
N=3421

j=0;
for i=1:ipskip:N;
   j = j + 1;
   infile = strcat(filebase1,num2str(i,'%04i'),filebase2,fileext);
   fp(:,:, j) = dlmread(infile);
end;
figure(3);box on;hold on;
subplot(3,1,1);box on;
size(squeeze(fp(idx(1),:,:)))
size([1:ipskip:N])
size(squeeze(fp(8,:,1)))
pcolor([1:ipskip:N],squeeze(fp(8,:,1)),squeeze(fp(idx(1),:,:)));
shading flat;
subplot(3,1,2);box on;
pcolor([1:ipskip:N],squeeze(fp(8,:,1)),squeeze(fp(idx(2),:,:)));
shading flat;
subplot(3,1,3);box on;
pcolor([1:ipskip:N],squeeze(fp(8,:,1)),squeeze(fp(idx(3),:,:)));
shading flat;
clear;

