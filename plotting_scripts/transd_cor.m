%
% Compute and plot correlation matrix
%

function [] = transd_cor(samplefile);

load(samplefile);
NPL = 6;

k = A(:,4);
kmin = min(k)
kmax = max(k)

%minlim = [ 0.15,  1450., 1.20, 0.001 ];
%maxlim = [ 7.50,  1700., 2.10, 1.000 ];
%minlimh = [1450., 1.20, 0.001 ];
%maxlimh = [1700., 2.10, 1.000 ];
minlim = [0.1 0.10 log(1.e6)  log(0.001) log(4.53e-3) log(1.e6)];
maxlim = [8.0 0.91 log(1.e11) log(0.05) log(0.4)   log(2.e10)];
minlimh = minlim(2:end);
maxlimh = maxlim(2:end);
ylabels = {'h1','c1','r1','a1','h2','c2','r2','a2',...
                  'h3','c3','r3','a3','h4','c4','r4','a4',...
                  'h5','c5','r5','a5','h6','c6','r6','a6',...
                  'h7','c7','r7','a7','h8','c8','r8','a8',...
                  'h9','c9','r9','a9','h9','c9','r9','a9'};
xlabels = {'h1','c1','r1','a1','h2','c2','r2','a2',...
                  'h3','c3','r3','a3','h4','c4','r4','a4',...
                  'h5','c5','r5','a5','h6','c6','r6','a6',...
                  'h7','c7','r7','a7','h8','c8','r8','a8',...
                  'h9','c9','r9','a9','h9','c9','r9','a9'};

for il = 1:kmax-kmin+1;
   ik = kmin+il-1
   idx = find(A(:,4) == ik);
   NFP = ((ik+1)*NPL)-1;
   idxh = [1:NPL:NFP-(NPL-1)]
   m  = A(idx,5:4+NFP);
   mh = m;mh(:,idxh)=[];

   ylab = fliplr(ylabels(1:NFP));
   xlab = xlabels(1:NFP);
   ylabh = ylabels(1:NFP);
   xlabh = xlabels(1:NFP);
   ylabh(idxh) = [];
   xlabh(idxh) = [];
   ylabh = fliplr(ylabh);

   minlim2 = repmat(minlim,1,ik);
   maxlim2 = repmat(maxlim,1,ik);
   minlim2 = [minlim2,minlim(2:end)];
   maxlim2 = [maxlim2,maxlim(2:end)];
   minlim2h = repmat(minlimh,1,ik+1);
   maxlim2h = repmat(maxlimh,1,ik+1);

   dlim = maxlim2-minlim2;
   dlimh = maxlim2h-minlim2h;
   npar = length(m(1,:));
   ndat = length(m);
   nparh = length(mh(1,:));
   ndath = length(mh);

   for i = 1:npar
     m(:,i) = (m(:,i) - minlim2(i))/dlim(i);
   end 
   for i = 1:nparh
     mh(:,i) = (mh(:,i) - minlim2h(i))/dlimh(i);
   end 

   mcov = zeros(npar,npar);
   mcor = zeros(size(mcov));
   mcorp = zeros(size(mcov));
   mcovh = zeros(nparh,nparh);
   mcorh = zeros(size(mcovh));
   mcorph = zeros(size(mcovh));

   mcov = cov(m);
   mcovh = cov(mh);

   for i = 1:npar
     for j = 1:npar
       mcor(i,j) = mcov(i,j)/sqrt(mcov(i,i)*mcov(j,j));
       if(mcor(i,j) < 0.7);mcor(i,j) = 0.;end;
     end
   end
   for i = 1:nparh
     for j = 1:nparh
       mcorh(i,j) = mcovh(i,j)/sqrt(mcovh(i,i)*mcovh(j,j));
       if(mcorh(i,j) < 0.7);mcorh(i,j) = 0.;end;
     end
   end
   kk = 1;
   for i = 1:npar
     for j = i:npar
       if(abs(mcor(i,j)) > 0.55 & abs(mcor(i,j)) < 1.0)
         a(kk,:) = [i j];
         mcorplot(kk)=mcor(i,j);
         kk = kk + 1;
       end
     end
   end
   kkh = 1;
   for i = 1:nparh
     for j = i:nparh
       if(abs(mcorh(i,j)) > 0.55 & abs(mcorh(i,j)) < 1.0)
         ah(kkh,:) = [i j];
         mcorploth(kk)=mcorh(i,j);
         kkh = kkh + 1;
       end
     end
   end
   mcorp = flipud(mcor);
   mcorph = flipud(mcorh);

   gc1=figure();
   hold on;box on;grid on;
   for i = 1:npar
%     plot(mcorp(i,:));
     offset = ((i-1)*2)+1;
     area([(npar-i)+1:npar],mcorp(i,(npar-i)+1:end)+offset,offset,'FaceColor',[0 0 0]);
     set(gca,'YLim',[0 2*npar+1],'YTick',[1:2:npar*2]);
     set(gca,'XLim',[1 npar],'XTick',[1:1:npar]);
     set(gca,'YTickLabel',ylab);
     set(gca,'XTickLabel',xlab);

   end
   gc2=figure();
   hold on;box on;grid on;
   for i = 1:nparh
     offset = ((i-1)*2)+1;
     area([(nparh-i)+1:nparh],mcorph(i,(nparh-i)+1:end)+offset,offset,'FaceColor',[0 0 0]);
     set(gca,'YLim',[0 2*nparh+1],'YTick',[1:2:nparh*2]);
     set(gca,'XLim',[1 nparh],'XTick',[1:1:nparh]);
     set(gca,'YTickLabel',ylabh);
     set(gca,'XTickLabel',xlabh);

   end
   clear m idx minlim2 maxlim2;
end

return;
