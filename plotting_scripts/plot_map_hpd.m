function []=plot_models();
%
% Loads: map.mat (from plot_hist)
%

set(0, 'DefaultFigurePaperPosition', [0 0 7 4]); 
dex = 0;
paper = 1;
showcore = 1;
errbar = 1;	%  0 == print cores only
		%  1 == print errbars on every 10th datum
barstep = 10;
max_depth = 2;
npar = 6;
col = {'+r' '+g'};
col = char(col);
cole = {'r' 'g'};
cole = char(cole);

lay_thick = 0.02 ;
znorm=lay_thick/2:lay_thick:1-lay_thick/2;
nlayer = round(max_depth/lay_thick);

nmod    = 10000;
nrbin   = 100;
site = 201;
if site == 4
    h1 = 1.5;
    nc = 2;
    data    = 'hol_site04_b.mat';
    sample  = 'site04_sample.mat';
    core = {'cores/site04_core5.mat' 'cores/site04_core5.mat'};
elseif site == 19
    h1 = 1.;
    nc = 1;
    data    = 'hol_site19.mat';
    sample  = 'site19_sample.mat';
    core = {'cores/site19_core5b.mat'};
elseif site == 191
    h1 = 1.;
    nc = 1;
    data    = 'hol_site19_b.mat';
    sample  = 'site19b_sample.mat';
    core = {'cores/site19_core5b.mat'};
elseif site == 20
    h1 = 1.;
    nc = 2;
    data    = 'hol_site20.mat';
    sample  = 'site20_sample.mat';
    core = {'cores/site20_core4b.mat' 'cores/site20_core4bq.mat'};
elseif site == 201
    h1 = 1.;
    nc = 2;
    data    = 'hol_site20_b.mat';
    sample  = 'site20_b_sample.mat';
    core = {'cores/site20_core4b.mat' 'cores/site20_core4bq.mat'};
end

plotfile1 = strrep(sample,'sample.mat','models_rand.eps');
mapfile = strrep(sample,'sample.mat','map.mat');

if errbar == 0
  plotfile2 = strrep(sample,'sample.mat','models_hpd.eps');
elseif errbar == 1
  plotfile2 = strrep(sample,'sample.mat','models_hpd_10.eps');
end

load(mapfile);
load(sample);
load(data);

%
% Use ML model
%
%xtrue = [1.9 1.32 1.50 0.80 1472.83  1465.87 0.30]';
%xtrue = [1.12 1.1989922  1.4990106  0.21987834  1471.061   1460.0173 0.3]';
%xml = xtrue;
%
% Use MAP model
%
xmap = [h1 xmap];
xml = xmap;

%
% Correct cores to in-situ measurements:
%
if(showcore == 1)
  for i = 1:nc
    clear tmp; clear tmp2; clear tmp3;
    tmp = load(core{1,i});
    tmp2=cell2mat(fieldnames(tmp));
    tmp3 = getfield(tmp,tmp2);
    C(i).depth = tmp3.depth;
    C(i).density = tmp3.density;
    C(i).density_ref = tmp3.density_ref;
    C(i).sound_speed = tmp3.sound_speed;
    C(i).sound_speed_ref = tmp3.sound_speed_ref;
    sc(i) = F(1).cw/C(i).sound_speed_ref;
    sr(i) = F(1).rw/C(i).density_ref;
    C(i).sound_speed = C(i).sound_speed/sc(i);
    C(i).density = C(i).density/sr(i);
    if(errbar == 1)
      C(i).Err_r = 2/100 * C(i).density;
      C(i).Err_c = 6 * ones(length(C(i).sound_speed),1);
    end
  end
end

m = A(:,2:end-1);
m = [ones(length(A),1)*h1 m];
disp('read sample');
if(nmod > length(m))
  nmod = length(m);
end

length(m)
%rand('state',265456);
Y = rand(nmod,1);
Y = floor(length(m) * Y);
if(min(Y) == 0)
  [YM I] = min(Y);
  Y(I) = 1;
end

X = m(Y,:);

h = X(:,1);
rhot = X(:,2);
rhob = X(:,3);
nu = X(:,4);
ct = X(:,5);
cb = X(:,6);
alpha = X(:,7);

c = zeros(1,length(znorm)+1);
rho = zeros(1,length(znorm)+1);

%figure(1);
%hold on;
for i = 1:nmod
  z = [0 znorm*h(i)];
  rho = rhot(i) + (sin(znorm .* pi/2)).^nu(i) .* (rhob(i) - rhot(i));
  rho = [rhot(i) rho];
  c = ct(i) + (cb(i) - ct(i)) .* znorm;
  c = [ct(i) c];
  a = alpha(i) * ones(size(znorm));
  a = [alpha(i) a];
end

%
% ML model:
%
h_ml = xml(1);
rhot_ml = xml(2);
rhob_ml = xml(3);
nu_ml = xml(4);
ct_ml = xml(5);
cb_ml = xml(6);
alpha_ml = xml(7);
%
% MAP model:
%
h_map = xmap(1);
rhot_map = xmap(2);
rhob_map = xmap(3);
nu_map = xmap(4);
ct_map = xmap(5);
cb_map = xmap(6);
alpha_map = xmap(7);

z_ml = [0 znorm] *h_ml;
rho_ml = rhot_ml + (sin([0 znorm] .* pi/2)).^nu_ml .* (rhob_ml - rhot_ml);
c_ml = ct_ml + (cb_ml - ct_ml) .* [0 znorm];

z_map = [0 znorm] *h_map;
rho_map = rhot_map + (sin([0 znorm] .* pi/2)).^nu_map .* (rhob_map - rhot_map);
c_map = ct_map + (cb_map - ct_map) .* [0 znorm];

%z_ml = znorm *h_ml;
%rho_ml(1) = rhot_ml;
%rho_ml(2:end) = rhot_ml + (sin(znorm .* pi/2)).^nu_ml .* (rhob_ml - rhot_ml);
%c_ml(1) = ct_ml;
%c_ml(2:end) = ct_ml + (cb_ml - ct_ml) .* znorm;
%a_ml = alpha_ml * ones(size(znorm));
%subplot('Position',[.12 .15 .4 .8]);
%plot(rho_ml,z_ml,'-+k','MarkerSize',8);
%plot(rho_map,z_map,'-+k','MarkerSize',8);
%legend('sample','core 5','core 6','ML estimate');

%subplot('Position',[.56 .15 .4 .8]);
%plot(c_ml,z_ml,'-+k','MarkerSize',8);
%plot(c_map,z_map,'-+k','MarkerSize',8);
%legend('sample','core 5','core 6','ML estimate');

%saveas(gca,plotfile1,'epsc2');

%
% Plot bounds for ML model:
%
z_fixed = [0.0:lay_thick:1.8];
r_prof = zeros(length(z_fixed),nmod);
c_prof = zeros(length(z_fixed),nmod);
for i = 1:nmod
  z(1) = 0;
  z(2:end) = znorm *h(i);
  rho(1) = rhot(i);
  rho(2:end) = rhot(i) + sin(znorm .* pi/2).^nu(i) .* (rhob(i) - rhot(i));
  c(1) = ct(i);
  c(2:end) = ct(i) + (cb(i) - ct(i)) .* znorm;
  
  k = 1;
  for j = 1:length(znorm)		% Interpolation loop
    
    r_m = rho(j); r_p = rho(j+1);
    c_m = c(j); c_p = c(j+1);
    z_m = z(j);   z_p = z(j+1);
    dr = r_p - r_m;
    dc = c_p - c_m;
    dz = z_p - z_m;
    r_grad = dr/dz;
    c_grad = dc/dz;
    if(k <= length(z_fixed))
      while(z_fixed(k) < z_p)
        r_prof(k,i) = r_m + r_grad*(z_fixed(k)-z_m);
        c_prof(k,i) = c_m + c_grad*(z_fixed(k)-z_m);
        k = k+1;
        if(k > length(z_fixed)) 
          break; 
        end
      end
    end
  end
end

%
% Calc 95% HPDs for each layer:
%
  nf_int_r = ninetyfive(r_prof,nmod,nrbin,z_fixed);
  nf_int_c = ninetyfive(c_prof,nmod,nrbin,z_fixed);

  figure(2);
  hold on;
  subplot('Position',[.12 .15 .4 .8]);
  box on;
  hold on;
  if(showcore == 1)
    for i = 1:nc
      if(errbar == 1)
        errorbarxy(C(i).density(1:barstep:end),C(i).depth(1:barstep:end),...
                   C(i).Err_r(1:barstep:end),...
                   zeros(size(C(i).Err_r(1:barstep:end))),cole(i,:),cole(i,:));
        plot(C(i).density,C(i).depth,col(i,:));
      else
        plot(C(i).density,C(i).depth,col(i,:));
      end
    end
  end
  plot(nf_int_r(:,1),z_fixed,'--k');
  plot(nf_int_r(:,2),z_fixed,'--k');
%  plot(rho_ml,z_ml,'-k.','MarkerSize',8);
  plot(rho_map,z_map,'-k.','MarkerSize',8);
  plot([rho_map(end) rho_map(end)],[z_map(end) 1.8],'-k.','MarkerSize',8);
  if(paper == 1)
    set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[0 0.5 1.0 1.5],...
        'YTickLabel',[0 0.5 1.0 1.5],'FontSize',16);
  else
    set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[0 0.5 1.0 1.5],...
        'YTickLabel',[0 0.5 1.0 1.5],'FontSize',16);
  end
  set(gca,'XLim',[1.15 1.65]);
  set(gca,'XTickLabel',[1.2 1.4 1.6],'XTick',[1.2 1.4 1.6]);
  xlabel('\rho [g/cm^3]');
  ylabel('depth [m]');
%return;
  subplot('Position',[.56 .15 .4 .8]);
  box on;
  hold on;
  if(showcore == 1)
    for i = 1:nc
      if(errbar == 1)
        errorbarxy(C(i).sound_speed(1:barstep:end), ... 
                   C(i).depth(1:barstep:end),...
                   C(i).Err_c(1:barstep:end),...
                   zeros(size(C(i).Err_c(1:barstep:end))),cole(i,:),cole(i,:));
        plot(C(i).sound_speed,C(i).depth,col(i,:));
      else
        plot(C(i).sound_speed,C(i).depth,col(i,:));
      end
    end
  end
  plot(nf_int_c(:,1),z_fixed,'--k');
  plot(nf_int_c(:,2),z_fixed,'--k');
%  plot(c_ml,z_ml,'-k.','MarkerSize',8);
  plot(c_map,z_map,'-k.','MarkerSize',8);
  plot([c_map(end) c_map(end)],[z_map(end) 1.8],'-k.','MarkerSize',8);
  set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[0 0.5 1.0 1.5],...
      'YTickLabel',[],'FontSize',16);
  set(gca,'XLim',[1450 1490]);
  set(gca,'XTickLabel',[1460 1470 1480],'XTick',[1460 1470 1480]);
  xlabel('c [m/s]');
%  ylabel('depth [m]');
  saveas(gca,plotfile2,'epsc2');
  hgsave('fig16.fig');
%  saveas(gca,plotfile2,'png');
return;

function [nf_int] = ninetyfive(r_prof,nmod,nrbin,z_fixed)
  for k = 1:length(z_fixed)
    j = 1;
    for i = 1:nmod
      if(r_prof(k,i) ~= 0)
        r_tmp(j) = r_prof(k,i);
        j = j + 1;
      end
    end
    xmin = min(r_tmp);
    xmax = max(r_tmp);
    xdiff = xmax - xmin;
    for i =0:nrbin
      edges(i+1) = xmin + xdiff/nrbin * i;
    end

    n1 = histc(r_tmp,edges);
    nintyfive = sum(n1)-5*sum(n1)/100;
    lb = 0;
    rb = 0;
 
    p = 1;
    for i = 1:length(n1)
      s = 0;
      stopj = 0;
      for j = i:length(n1)
        if(stopj == 0)
          s = s + n1(j);
          if(s >= nintyfive)
            lb(p) = edges(i);
            rb(p) = edges(j);
            p = p+1;
            stopj = 1;
          end
        end
      end
    end
    [foo,nf_index] = min(abs(lb-rb));
    nf_int(k,1) = lb(nf_index);
    nf_int(k,2) = rb(nf_index);
  end
return;
