function debug_spher;

%
% Load data averaged in Fortran code
%
sphera = load('spher_refldBa.dat');
BLa = (sphera(1:end-1,:));
anga = (sphera(end,:));
spherb = load('spher_refldBb.dat');
BLb = (spherb(1:end-1,:));
angb = (spherb(end,:));
%figure();
%hold on;box on;
%plot(angd,BLd,'-r')

%load pack1.txt
sim = load('dag3_02_1b.txt');
BLe = sim(1:end-1,:);
ange = sim(end,:);
% plot(angd,pack1,'-k')
% BLe = pack1(1:5,:);
    figure();
    hold on;box on;
for ifrq=1:9

    ifrq
    subplot(3,3,ifrq);
    hold on;box on;
    plot(anga,BLa(ifrq,:),'-b')
    plot(angb,BLb(ifrq,:),'-r')
    plot(ange,BLe(ifrq,:),'--k')
    legend('oastl true','oastl MAP','oasp (timeseries)')

end

%ange = angd;
%save angs ange

return;
