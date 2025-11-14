function [] = test_residual_sample(filename);

% filename = tmp_residuals.dat
A=load(filename);
A = A(:,:);
meanA=mean(A,2);
meanAA=repmat(meanA,1,size(A,2));
B=A-meanAA;

NSMP = size(A,1);
NDAT = size(A,2);
NDATB = 35;
NDAT2 = 4000;

NDAT
NDATB

BA = B(:,1:NDATB-1);
BB = B(:,NDATB:NDAT);

size(BA)
size(BB)

ipassa = 0;ifaila = 0;ikpassa = 0;ikfaila = 0;iadpassa = 0;iadfaila = 0;
ilpassa = 0;ilfaila = 0;il2passa = 0;il2faila = 0;iswpassa = 0;iswfaila = 0;
ipassb = 0;ifailb = 0;ikpassb = 0;ikfailb = 0;iadpassb = 0;iadfailb = 0;
ilpassb = 0;ilfailb = 0;il2passb = 0;il2failb = 0;iswpassb = 0;iswfailb = 0;

xlim = 5.;
dx = 0.2;
edges = [-xlim-dx/2.:dx:xlim+dx/2.];
ctrs = edges(1:end-1) + dx./2.;
for i=1:NSMP;
 
   if(rem(i,1000) == 0);disp(i);end;

   [ha(i),pa(i)]=runstest(BA(i,:),0.);
   if(ha(i) == 1);ifaila = ifaila+1;end;
   if(ha(i) == 0);ipassa = ipassa+1;end;
   if(size(BB,2)>10);
     [hb(i),pb(i)]=runstest(BB(i,:),0.);
     if(hb(i) == 1);ifailb = ifailb+1;end;
     if(hb(i) == 0);ipassb = ipassb+1;end;
   end;

   [hka(i),pka(i)]=kstest(BA(i,:),[],0.05,'smaller');
   if(hka(i) == 1);ikfaila = ikfaila+1;end;
   if(hka(i) == 0);ikpassa = ikpassa+1;end;
   if(size(BB,2)>10);
     [hkb(i),pkb(i)]=kstest(BB(i,:),[],0.05,'smaller');
     if(hkb(i) == 1);ikfailb = ikfailb+1;end;
     if(hkb(i) == 0);ikpassb = ikpassb+1;end;
   end;

   [hada(i),pada(i)]=AnDartest(BA(i,:));
   if(hada(i) == 1);iadfaila = iadfaila+1;end;
   if(hada(i) == 0);iadpassa = iadpassa+1;end;
   if(size(BB,2)>10);
     [hadb(i),padb(i)]=AnDartest(BB(i,:));
     if(hadb(i) == 1);iadfailb = iadfailb+1;end;
     if(hadb(i) == 0);iadpassb = iadpassb+1;end;
   end;

   [hswa(i),pswa(i)]=swtest(BA(i,:));
   if(hswa(i) == 1);iswfaila = iswfaila+1;end;
   if(hswa(i) == 0);iswpassa = iswpassa+1;end;
   if(size(BB,2)>10);
     [hswb(i),pswb(i)]=swtest(BB(i,:));
     if(hswb(i) == 1);iswfailb = iswfailb+1;end;
     if(hswb(i) == 0);iswpassb = iswpassb+1;end;
   end;

   C = B(i,:)/std(B(i,:));
   D = randn(1,NDAT2);
   n1 = histc(C,edges);
   n2 = histc(D,edges);
   nn1(i,:) = n1(1:end-1);
   nn2(i,:) = n2(1:end-1);
   n_norm1(i,:) = nn1(i,:)/(NDAT*dx);
   n_norm2(i,:) = nn2(i,:)/(NDAT2*dx);

end;

mxprob = max(max(n_norm1))+.05;
prob = [0:mxprob/100:mxprob];

nn = zeros(length(prob),length(ctrs));

for i=1:NSMP;
   for j=1:length(ctrs);

      if(n_norm1(i,j) > 0); 
         idx = find(n_norm1(i,j)>=prob);
         nn(idx,j) = nn(idx,j) + 1;
      end;

   end;
end;

figure(1);hold on;box on;
nn = nn/NSMP;
pcolor(ctrs,prob,nn);shading flat;

xx = [-6:.1:6];
nd = 1/sqrt(2*pi)*exp(-(xx.^2)/2);
plot(xx,nd,'-k')
plot(xx,nd,'--w')
set(gca,'XLim',[-xlim xlim],'YLim',[0 1.0],'CLim',[0 1])
colorbar;

%figure();hold on;box on;
%stem(h);

%figure();hold on;box on;
%hist(p,100);

fid = fopen('residual_statistics.dat','w');
fprintf(fid,'Runstest A:      %12.4f %% passed\n',ipassa/(ipassa+ifaila));
fprintf(fid,'Runstest B:      %12.4f %% passed\n',ipassb/(ipassb+ifailb));
fprintf(fid,'KS test A:       %12.4f %% passed\n',ikpassa/(ikpassa+ikfaila));
fprintf(fid,'KS test B:       %12.4f %% passed\n',ikpassb/(ikpassb+ikfailb));
fprintf(fid,'AD test A:       %12.4f %% passed\n',iadpassa/(iadpassa+iadfaila));
fprintf(fid,'AD test B:       %12.4f %% passed\n',iadpassb/(iadpassb+iadfailb));
fprintf(fid,'SW test A:       %12.4f %% passed\n',iswpassa/(iswpassa+iswfaila));
fprintf(fid,'SW test B:       %12.4f %% passed\n',iswpassb/(iswpassb+iswfailb));
fclose(fid);

return;
