function []=test_TDSA();

d(:,:,1) = load('fort.71');
d(:,:,2) = load('fort.72');
d(:,:,3) = load('fort.73');
d(:,:,4) = load('fort.74');
d(:,:,5) = load('fort.75');
d(:,:,6) = load('fort.76');
d(:,:,7) = load('fort.77');

col = {'r','b','k','--r','--b','--k',':r'};
col = char(col);

for i=1:7;

   figure(1);hold on;box on;
   subplot(3,1,1);hold on;box on;
   plot(d(:,1,i),col(i,:));
   subplot(3,1,2);hold on;box on;
   plot(d(:,3,i),col(i,:));
   subplot(3,1,3);hold on;box on;
   plot(d(:,2,i),col(i,:));

end;

return;
