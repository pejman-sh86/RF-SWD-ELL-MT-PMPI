function []=compute_prior_vol();

clb = 1.500;
chb = 2.000;
rlb = 0.99;
rhb = 2.1;

r1=1.2;
r2=2.1;
c1=1.500;
c2=2.000;

N=2000.;

r = [r1:(r2-r1)/N:r2];
cl = (1.54-0.907.*r+0.3700.*r.^1.88)*1.4904;
ch=(1.49-0.907*r+0.4895*r.^1.85)*1.4914;

cl(find(cl<clb))=clb;
cl(find(cl>chb))=chb;
ch(find(ch<clb))=clb;
ch(find(ch>chb))=chb;

idx1=find(r<rlb);
r(idx1)=[];
cl(idx1)=[];
ch(idx1)=[];


idx2=find(r>rhb);
r(idx2)=[];
cl(idx2)=[];
ch(idx2)=[];

idx3=find(ch-cl==0);
r(idx3)=[];
cl(idx3)=[];
ch(idx3)=[];

vol = sum(ch-cl)*(r2-r1)/N
vol2 = (chb-clb)*(rhb-rlb)

figure();hold on;
plot(r,cl)
plot(r,ch)
xlabel('Desity');
ylabel('Velocity');
plot(1.5,1.6953000,'xk');
set(gca,'XLim',[rlb rhb],'YLim',[clb chb])
%figure();hold on;
%plot(r,ch-cl)
