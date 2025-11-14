function []=reformat_seacable();

NFREQ    = 25;
NPING    = 50;
dat      = load('dispersion_PH1-80Greentw4p5secfp68-6.dat');
filebase = 'seacable_p';
fileext  = '.txt';

for iping = 1:NPING;

   idx1 = find(dat(1:NFREQ,iping) ~= 0);
   idx2 = find(dat(NFREQ+1:2*NFREQ,iping) ~= 0);
   NDAT = length(idx1);
   freq = dat(idx1,iping);
   vsgr = dat(idx2+NFREQ,iping)

   savedat = [NDAT,NDAT;freq,vsgr];
   filename = strcat(filebase,num2str(iping,'%03i'),fileext);
   save(filename,'savedat','-ascii');

   clear idx1 idx2 freq vsgr savedat;
end;

return;
