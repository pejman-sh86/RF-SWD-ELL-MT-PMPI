function [] = test_residual_sample(filename);

% filename = tmp_residuals.dat
A=load(filename);
A = A(:,:);
meanA=mean(A,2);
meanAA=repmat(meanA,1,size(A,2));
B=A-meanAA;

NSMP = size(A,1);
NDAT = size(A,2);
NDAT2 = 4000;

ipass = 0;ifail = 0;ikpass = 0;ikfail = 0;iadpass = 0;iadfail = 0;
ilpass = 0;ilfail = 0;il2pass = 0;il2fail = 0;iswpass = 0;iswfail = 0;

xlim = 5.;
dx = 0.2;
edges = [-xlim-dx/2.:dx:xlim+dx/2.];
ctrs = edges(1:end-1) + dx./2.;
for i=1:NSMP;
 
   if(rem(i,1000) == 0);disp(i);end;

%   D = randn(1,NDAT);
   [h(i),p(i)]=runstest(B(i,:),0.);
%   [h(i),p(i)]=runstest(D,0.);
   if(h(i) == 1);ifail = ifail+1;end;
   if(h(i) == 0);ipass = ipass+1;end;

   [hk(i),pk(i)]=kstest(B(i,:),[],0.05,'smaller');
%   [hk(i),pk(i)]=kstest(D,[],0.05,'smaller');
   if(hk(i) == 1);ikfail = ikfail+1;end;
   if(hk(i) == 0);ikpass = ikpass+1;end;

   [had(i),pad(i)]=AnDartest(B(i,:));
%   [had(i),pad(i)]=AnDartest(D);
   if(had(i) == 1);iadfail = iadfail+1;end;
   if(had(i) == 0);iadpass = iadpass+1;end;

   [hl(i),pl(i)]=lillietest(B(i,:));
%   [hl(i),pl(i)]=lillietest(D);
%   [hl(i),pl(i)]=lillietest(B(i,:),0.05,'norm',1e-3);
   if(hl(i) == 1);ilfail = ilfail+1;end;
   if(hl(i) == 0);ilpass = ilpass+1;end;

   [hl2(i),pl2(i)]=lillietest(A(i,:));
%   [hl2(i),pl2(i)]=lillietest(D);
%   [hl2(i),pl2(i)]=lillietest(A(i,:),0.05,'norm',1e-3);
   if(hl2(i) == 1);il2fail = il2fail+1;end;
   if(hl2(i) == 0);il2pass = il2pass+1;end;

   [hsw(i),psw(i)]=swtest(B(i,:));
%   [hsw(i),psw(i)]=swtest(D);
   if(hsw(i) == 1);iswfail = iswfail+1;end;
   if(hsw(i) == 0);iswpass = iswpass+1;end;

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
fprintf(fid,'Runstest:      %12.4f %% passed\n',ipass/(ipass+ifail));
fprintf(fid,'KS test:       %12.4f %% passed\n',ikpass/(ikpass+ikfail));
fprintf(fid,'AD test:       %12.4f %% passed\n',iadpass/(iadpass+iadfail));
fprintf(fid,'Lillie test A: %12.4f %% passed\n',ilpass/(ilpass+ilfail));
fprintf(fid,'Lillie test B: %12.4f %% passed\n',il2pass/(il2pass+il2fail));
fprintf(fid,'SW test:       %12.4f %% passed\n',iswpass/(iswpass+iswfail));
fclose(fid);

return;
