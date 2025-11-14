obs=load('input_obsdat.dat');
pred=dlmread('maule_data_replica.dat');
scale = 0.10;
rawnoise = randn(size(obs));
noise = zeros(size(obs));
covdat=dlmread('maule_data_icovmat.dat');
nrow = sum(pred(2,1:pred(1,1)))
for istn=1:pred(1,1);
  ie=sum(pred(2,1:istn));
  is=ie-pred(2,istn)+1;
  for irow=1:pred(2,istn);
    F(istn).Cdi(irow,:) = covdat(sum(pred(2,1:istn))-pred(2,istn)+irow,1:pred(2,istn));
    F(istn).Cd(irow,:)  = covdat(sum(pred(2,1:istn))-pred(2,istn)+irow+nrow,1:pred(2,istn));
  end;
  F(istn).L = chol(F(istn).Cd,'lower');
  disp([is,ie])
  %sd(istn)=scale*mean(abs(pred(3,is:ie)));
  sd(istn)=sqrt(F(istn).Cd(1,1));
  noise(is:ie) = sd(istn)*rawnoise(is:ie);
  cornoise(is:ie) = F(istn).L*rawnoise(is:ie);
end;

%sim = pred(3,:)+noise';
sim = pred(3,:)+cornoise;
sim = sim';

figure();hold on;box on;
plot(noise,'-b');
plot(cornoise,'--r');
figure();hold on;box on;
plot(pred(3,:),'-b');
plot(sim,'.k');
save('input_obsdat.dat','sim','-ascii');
