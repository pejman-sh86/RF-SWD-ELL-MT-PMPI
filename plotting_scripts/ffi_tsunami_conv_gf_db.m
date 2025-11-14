function []=ffi_tsunami_conv_gf_db(datafile);

newfile = 'new_db.h5';

NPATCH = 576;
NSTN = 22;

stn = {'1002','1006','202','203','205','21401','21413','21418',...
       '2672','2673','5741','5742','5861','5862','602','613',...
       '801','802','803','804','806','807'};

b1=[0;boxcar(59);0];
b1=b1/sum(b1);
for ip=1:NPATCH;
  for ist=1:NSTN;
    path = strcat('/computed/gf',num2str(ip-1),'/',stn(ist),'/data');
    d(:,:,ip,ist) = h5read(datafile,char(path));
    d2(:,:,ip,ist) = d(:,:,ip,ist);
    x = conv(squeeze(d(2,:,ip,ist)),b1);
    d2(2,:,ip,ist) = x(16:9001+15);
    
    h5create(newfile,char(path),size(d2(:,:,ip,ist)));
    h5write(newfile,char(path),d2(:,:,ip,ist));
    
    %figure();hold on;
    %plot(d(1,:,ip,ist),d(2,:,ip,ist));
    %plot(d2(1,15:end,ip,ist),d2(2,15:end,ip,ist));
    %save tmp.mat
    %stop
  end;
end;


return;



