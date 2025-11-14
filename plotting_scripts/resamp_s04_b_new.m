function [] = resamp();

outfile  = strcat('x_s04_1_3_32_0layb.mat');
outfilea = strcat('x_s04_1_3_32_0layb.txt');

nang_re = 130;
load blde4_2_ping_ida_inv.mat;

%idx = [6:13];
idx = [6:16];
nf = length(idx);
theta = xde4.ang(1:132);
R = xde4.bl(idx,1:132);
R = R';

%minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 0.001];
%maxlim = [2.0, 1550, 1550, 1.6, 1.8, 0.500, 0.500];
%mtrue  = [1.0, 1500, 1500, 1.4, 1.55,0.100, 0.1005];

minlim = [0.5, 1450, 1450, 1.2, 1.3, 0.001, 1450, 1.3, 0.001];
maxlim = [2.0, 1550, 1550, 1.6, 1.8, 1.000, 1550, 1.8, 1.000];
mtrue  = [1.0, 1500, 1500, 1.4, 1.55,0.500, 1500, 1.5, 0.5005];

%minlim = [];
%maxlim = [];
%mtrue  = [];

%thdiff = diff(theta);
%idx2 = find(thdiff == 0);
%theta(idx2) = [];
%R(idx2,:) = [];

theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];
te_diff = 3*mean(diff(theta_even));

%Rex = -1*(isnan(R)-1);
%
%for ifr = 1:nf
%   R(find(Rex(:,ifr) == 0),ifr) = 0;
%end
%clear Rex;
Rex = ones(length(theta_even),nf);
for ifr = 1:nf

   for j=1:length(theta_even)

      idx_ave = find((theta<theta_even(j)+(te_diff/2.)) & ...
                     (theta >= theta_even(j)-(te_diff/2.)) & ...
                     (isnan(R(:,ifr)') == 0));
      R_even(j,ifr) = mean(R(idx_ave,ifr));

   end;


   Rex(:,ifr) = -1*(isnan(R_even(:,ifr))-1);
   R_even(find(Rex(:,ifr) == 0),ifr) = 0;

end

Rnan = (isnan(R));
for ifr = 1:nf
   R_tmp = R(:,ifr);
   theta_tmp = theta;
   idx_int = find(Rnan(:,ifr) == 1)
   R_tmp(idx_int) = [];
   theta_tmp(idx_int) = [];
   R_int(:,ifr)   = interp1(theta_tmp,R_tmp,theta_even);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = 1511.37;
F(1).rw = 1.029;

F(1).freq = xde4.pref(idx,1)
for ifr = 1:nf;
    F(ifr).dat = R_even(:,ifr);
    F(ifr).ang = theta_even;
end


F(1).minlim = minlim;
F(1).maxlim = maxlim;
F(1).mtrue  = mtrue;
F(1).Rex    = Rex;
F(1).nmod   = length(F(1).maxlim);

save(outfile,'F');
save(outfilea,'minlim','-ascii');
save(outfilea,'maxlim','-ascii','-append');
save(outfilea,'mtrue','-ascii','-append');
for i = 1:nf

    tmp = F(i).dat';
    save(outfilea,'tmp','-ascii','-append');

end
ang = F(1).ang(1:size(Rex,1));
save(outfilea,'ang','-ascii','-append');

for i = 1:nf

    tmp = F(1).Rex(:,i)';
    save(outfilea,'tmp','-ascii','-append');

end

for i=1:nf;
   subplot(3,4,i);hold on;box on;
   plot(F(1).ang,F(i).dat(:),'+k');
   text(70,30,num2str(F(1).freq(i)));
   set(gca,'XLim',[0 90],'YLim',[0 40]);
end;

%for i=1:8;
%   subplot(3,3,i);hold on;box on;
%   plot(theta_even,R_int(:,i),'.b');
%   text(70,30,num2str(F(1).freq(i)));
%   set(gca,'XLim',[0 90],'YLim',[0 40]);
%end;


%for i=1:8;
%   subplot(3,3,i);hold on;box on;
%   plot(xde4.ang,xde4.bl(i+5,:),'or');
%   text(70,30,num2str(xde4.pref(i+5,1)));
%   set(gca,'XLim',[0 90],'YLim',[0 40]);
%end;

return;
