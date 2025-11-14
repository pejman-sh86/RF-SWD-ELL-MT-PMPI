load tohoku_data_sample.mat

Nsub=50000;

nx = 2;
ny = 1;
ML = .03;
MR = .03;
MB = .1;
MT = .03;
SP = .02;
PAD = 0;
FNT = 14;
inter = 20;

figure();hold on;box on;
subaxis(ny,nx,2,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);
hold on; box on;
for isub = 1:4;
    [a(:,isub),b]=hist(A((isub-1)*Nsub+1:isub*Nsub,2),[100:150]);
    a(:,isub)=a(:,isub)/trapz(b,a(:,isub));
    [xx,yy]=stairs(b,a(:,isub));
    stairs(b,a(:,isub));
end;
xlabel('No. coefficients')
set(gca,'YTickLabel',[]);
set(gca,'XLim',[100,135]);

subaxis(ny,nx,1,'Spacing',SP,'Padding',PAD,'ML', ML,'MR',MR,'MB',MB,'MT',MT);
set(gca,'FontSize',FNT);
hold on; box on;
clear a b;
for isub = 1:4;
    [a(:,isub),b]=hist(A((isub-1)*Nsub+1:isub*Nsub,1),100);
    a(:,isub)=a(:,isub)/trapz(b,a(:,isub));
    [xx,yy]=stairs(b,a(:,isub));
    stairs(b,a(:,isub));
end;
xlabel('log likelihood')
ylabel('Probability')
set(gca,'YTickLabel',[]);
set(gca,'XLim',[-1280,-1155]);

