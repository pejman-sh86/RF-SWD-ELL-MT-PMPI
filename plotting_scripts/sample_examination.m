function [] = sample_examination();

load('x_s13_2_10_50_3lay_sample.mat');
A1 = A(1:2:end,:);
A2 = A(2:2:end,:);
ncpu = 8;
mkeep = 20;
nmod = length(A1(:,1))

j = 1;
for i = 1:80:nmod

   A11(j:j+mkeep-1,:) = A1([i:i+mkeep-1],:);
   A12(j:j+mkeep-1,:) = A1([i:i+mkeep-1]+20,:);
   A13(j:j+mkeep-1,:) = A1([i:i+mkeep-1]+40,:);
   A14(j:j+mkeep-1,:) = A1([i:i+mkeep-1]+60,:);
   A21(j:j+mkeep-1,:) = A2([i:i+mkeep-1],:);
   A22(j:j+mkeep-1,:) = A2([i:i+mkeep-1]+20,:);
   A23(j:j+mkeep-1,:) = A2([i:i+mkeep-1]+40,:);
   A24(j:j+mkeep-1,:) = A2([i:i+mkeep-1]+60,:);

   j = j+mkeep;

end;

nx = 5;
ny = 5;
nsubfig = nx*ny;
xim = 0.03;
yim = 0.1/ny;
xymarg = [0.04 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

nbin=20;
nfp = size(A,2);
for i=1:nfp;[n1(1,:,i),x(1,:,i)] = hist(A11(:,i),nbin);
   n1(1,:,i) = n1(1,:,i)/(sum(n1(1,:,i))*(max(A11(:,i))-min(A11(:,i)))/nbin);
end;
for i=1:nfp;[n1(2,:,i),x(2,:,i)] = hist(A12(:,i),nbin);
   n1(2,:,i) = n1(2,:,i)/(sum(n1(2,:,i))*(max(A12(:,i))-min(A12(:,i)))/nbin);
end;
for i=1:nfp;[n1(3,:,i),x(3,:,i)] = hist(A13(:,i),nbin);
   n1(3,:,i) = n1(3,:,i)/(sum(n1(3,:,i))*(max(A13(:,i))-min(A13(:,i)))/nbin);
end;
for i=1:nfp;[n1(4,:,i),x(4,:,i)] = hist(A14(:,i),nbin);
   n1(4,:,i) = n1(4,:,i)/(sum(n1(4,:,i))*(max(A14(:,i))-min(A14(:,i)))/nbin);
end;
for i=1:nfp;[n1(5,:,i),x(5,:,i)] = hist(A21(:,i),nbin);
   n1(5,:,i) = n1(5,:,i)/(sum(n1(5,:,i))*(max(A21(:,i))-min(A21(:,i)))/nbin);
end;
for i=1:nfp;[n1(6,:,i),x(6,:,i)] = hist(A22(:,i),nbin);
   n1(6,:,i) = n1(6,:,i)/(sum(n1(6,:,i))*(max(A22(:,i))-min(A22(:,i)))/nbin);
end;
for i=1:nfp;[n1(7,:,i),x(7,:,i)] = hist(A23(:,i),nbin);
   n1(7,:,i) = n1(7,:,i)/(sum(n1(7,:,i))*(max(A23(:,i))-min(A23(:,i)))/nbin);
end;
for i=1:nfp;[n1(8,:,i),x(8,:,i)] = hist(A24(:,i),nbin);
   n1(8,:,i) = n1(8,:,i)/(sum(n1(8,:,i))*(max(A24(:,i))-min(A24(:,i)))/nbin);
end;

for i=1:nfp;[n1(1,:,i),x2(1,:,i)] = hist(A1(:,i),nbin);
   n2(1,:,i) = n1(1,:,i)/(sum(n1(1,:,i))*(max(A1(:,i))-min(A1(:,i)))/nbin);
end;
for i=1:nfp;[n1(2,:,i),x2(2,:,i)] = hist(A2(:,i),nbin);
   n2(2,:,i) = n1(2,:,i)/(sum(n1(2,:,i))*(max(A2(:,i))-min(A2(:,i)))/nbin);
end;

for i=1:nfp;[n1(3,:,i),x2(3,:,i)] = hist(A(:,i),nbin);
   n2(3,:,i) = n1(3,:,i)/(sum(n1(3,:,i))*(max(A(:,i))-min(A(:,i)))/nbin);
end;

figure(1);
for i=1:nfp;

   subplot(5,6,i);
   hold on; box on;
   stairs(x(1,:,i),n1(1,:,i),'k');
   stairs(x(2,:,i),n1(2,:,i),'b');
   stairs(x(3,:,i),n1(3,:,i),'r');
   stairs(x(4,:,i),n1(4,:,i),'--k');
   stairs(x(5,:,i),n1(5,:,i),'--b');
   stairs(x(6,:,i),n1(6,:,i),'--r');
   stairs(x(7,:,i),n1(7,:,i),':k');
   stairs(x(8,:,i),n1(8,:,i),':b');

end;

figure(2);
isubfig = 1;
for i=1:nfp;

   if(i == nfp-2) 
      isubfig = isubfig + 1;
   end;
   subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
   hold on;box on;

   stairs(x2(1,:,i),n2(1,:,i),'b');
   stairs(x2(2,:,i),n2(2,:,i),'r');
   stairs(x2(3,:,i),n2(3,:,i),':k');
   
   isubfig = isubfig +1;

end;

save bla;
return;
