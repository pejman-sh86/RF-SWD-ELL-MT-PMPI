function []=test_resample();

A=load('resampling.txt');

j = 0;
for i=1:5;
   figure();

   subplot(3,2,1);
   j = j + 1;
   hist(A((j-1)*1000+1:j*1000,2),50);
   set(gca,'XLim',[150 270]);

   subplot(3,2,3);
   hist(A((j-1)*1000+1:j*1000,4),50);
   set(gca,'XLim',[150 270]);

   subplot(3,2,5);
   hist(A((j-1)*1000+1:j*1000,6),50);
   set(gca,'XLim',[150 270]);


   subplot(3,2,2);
   j = j + 1;
   hist(A((j-1)*1000+1:j*1000,2),50);
   set(gca,'XLim',[150 270]);

   subplot(3,2,4);
   hist(A((j-1)*1000+1:j*1000,4),50);
   set(gca,'XLim',[150 270]);

   subplot(3,2,6);
   hist(A((j-1)*1000+1:j*1000,6),50);
   set(gca,'XLim',[150 270]);

end;

return;
