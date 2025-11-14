function plot_models;
%
% Loads: map.mat (from plot_hist)
%

set(0, 'DefaultFigurePaperPosition', [0 0 7 4]); 

nmod = 10000;
nrbin = 100;
sample1 = 'sample1.dat';
sample2 = 'sample2.dat';
plotfile1 = strrep(sample1,'sample1.dat','models_rand.eps');
plotfile2 = strrep(sample1,'sample1.dat','models_hpd.eps');
A = load(sample1);
B = load(sample2);


m = [A(:,2:8) ; B(:,2:8)];

C = load('core5.dat');
D = load('core6.dat');
load('map.mat');
xtrue = [2.0930214      1.3508944      1.4841878     0.73839097 ...
         1472.8123 1467.5688     0.25152184]';
xmap = xtrue;

%
% Correct cores to in-situ measurements:
%
sc_5 = 1.003551;
sc_6 = 1.010277;
sr_5 = 1.006025;
sr_6 = 0.996599;

Y = rand(nmod,1);
Y = floor(length(m) * Y);
X = m(Y,:);

h = X(:,1);
rhot = X(:,2);
rhob = X(:,3);
nu = X(:,4);
ct = X(:,5);
cb = X(:,6);
alpha = X(:,7);

znorm = zeros(10,1);
lay_thick = 0.1;
znorm=lay_thick/2:lay_thick:1-lay_thick/2;

c = zeros(size(znorm));
rho = zeros(size(znorm));

figure(1);
hold on;
for i = 1:nmod
  z = [0 znorm *h(i)];
  rho = rhot(i) + (sin(znorm .* pi/2)).^nu(i) .* (rhob(i) - rhot(i));
  rho = [rhot(i) rho];
  c = ct(i) + (cb(i) - ct(i)) .* znorm;
  c = [ct(i) c];
  a = alpha(i) * ones(size(znorm));
  a = [alpha(i) a];
  subplot('Position',[.12 .15 .4 .8]);
  box on;
  hold on;
  plot(rho,z);
  set(gca,'YDir','reverse');
  subplot('Position',[.56 .15 .4 .8]);
  box on;
  hold on;
  plot(c,z);
  set(gca,'YDir','reverse');
end

%
% Core data:
%
subplot('Position',[.12 .15 .4 .8]);
plot(C(:,3)/sr_5,C(:,1),'+g');
plot(D(:,3)/sr_6,D(:,1),'+r');
set(gca,'YDir','reverse','YLim',[0 1.5],'FontSize',16);
set(gca,'XLim',[1.2 1.6]);
set(gca,'XTickLabel',[1.2 1.4 1.6],'XTick',[1.2 1.4 1.6]);
xlabel('\rho [g/cm^3]');
ylabel('depth [m]');

subplot('Position',[.56 .15 .4 .8]);
plot(C(:,2)/sc_5,C(:,1),'+g');
plot(D(:,2)/sc_6,D(:,1),'+r');
set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[],'FontSize',16);
set(gca,'XLim',[1450 1480]);
set(gca,'XTickLabel',[1460 1470],'XTick',[1460 1470]);
xlabel('c [m/s]');

%
% ML model:
%
h_ml = xmap(1);
rhot_ml = xmap(2);
rhob_ml = xmap(3);
nu_ml = xmap(4);
ct_ml = xmap(5);
cb_ml = xmap(6);
alpha_ml = xmap(7);

z_ml = znorm *h_ml;
rho_ml = rhot_ml + (sin(znorm .* pi/2)).^nu_ml .* (rhob_ml - rhot_ml);
c_ml = ct_ml + (cb_ml - ct_ml) .* znorm;
a_ml = alpha_ml * ones(size(znorm));
subplot('Position',[.12 .15 .4 .8]);
plot(rho_ml,z_ml,'-+k','MarkerSize',8);
%legend('sample','core 5','core 6','ML estimate');

subplot('Position',[.56 .15 .4 .8]);
plot(c_ml,z_ml,'-+k','MarkerSize',8);
%legend('sample','core 5','core 6','ML estimate');

%subplot(2,2,3)
%plot(a_ml,z_ml,'-+k','MarkerSize',8);
%legend('sample','core 5','core 6','ML estimate');

saveas(gca,plotfile1,'epsc2');

%
% Plot bounds for ML model:
%
z_fixed = [0.0:0.1:1.8];
r_prof = zeros(length(z_fixed),nmod);
c_prof = zeros(length(z_fixed),nmod);
for i = 1:nmod
  z = znorm *h(i);
  rho = rhot(i) + sin(znorm .* pi/2).^nu(i) .* (rhob(i) - rhot(i));
  c = ct(i) + (cb(i) - ct(i)) .* znorm;
  
  k = 1;
  for j = 1:9		% Interpolation loop
    
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
  plot(nf_int_r(:,1),z_fixed,'--k');
  plot(nf_int_r(:,2),z_fixed,'--k');
  plot(rho_ml,z_ml,'-k.','MarkerSize',8);
  plot(C(:,3)/sr_5,C(:,1),'+g');
  plot(D(:,3)/sr_6,D(:,1),'+r');
  set(gca,'YDir','reverse','YLim',[0 1.5],'FontSize',16);
  set(gca,'XLim',[1.2 1.6]);
  set(gca,'XTickLabel',[1.2 1.4 1.6],'XTick',[1.2 1.4 1.6]);
  xlabel('\rho [g/cm^3]');
  ylabel('depth [m]');
%return;
  subplot('Position',[.56 .15 .4 .8]);
  box on;
  hold on;
  plot(nf_int_c(:,1),z_fixed,'--k');
  plot(nf_int_c(:,2),z_fixed,'--k');
  plot(c_ml,z_ml,'-k.','MarkerSize',8);
  plot(C(:,2)/sc_5,C(:,1),'+g');
  plot(D(:,2)/sc_6,D(:,1),'+r');
  set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[],'FontSize',16);
  set(gca,'XLim',[1450 1480]);
  set(gca,'XTickLabel',[1460 1470],'XTick',[1460 1470]);
  xlabel('c [m/s]');
%  ylabel('depth [m]');
  saveas(gca,plotfile2,'epsc2');
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
