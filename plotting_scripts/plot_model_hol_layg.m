function [] = plot_model_hol_lay();

set(0, 'DefaultFigurePaperPosition', [0 0 10 6]);

hmax = 10;

hpdfile   = 'tmp_ci_hpds.txt';
mapfile   = 'tmp_ci_map.txt';
corefile  = 'core15.mat';
plotfile1 = 'tmp_results.eps';

hpds = load(hpdfile);

ishade = 1;
icore = 1;

if(icore == 1)
    load(corefile);
end;

M(1).m = load(mapfile);

npl = 4;

col = {'-k' '-r' '-r' '--k' '--r' ':k'};
col = char(col);

nx = 3
ny = 1
xim = 0.06;
yim = 0.06;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

hf1=figure(1);hold on; box on;
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
set(gca,'Fontsize',14,'XLim',[1450 1700]);
set(gca,'XTickLabel',[1500 1600 1700]);
set(gca,'XTick',[1500 1600 1700]);
hold on; box on;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
set(gca,'Fontsize',14,'XLim',[1.4 2.4]);
set(gca,'XTickLabel',[1.4 1.6 1.8 2.0 2.2 2.4]);
set(gca,'XTick',[1.4 1.6 1.8 2.0 2.2 2.4]);
set(gca,'YTickLabel',[]);
hold on; box on;

h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
set(gca,'Fontsize',14);
set(gca,'Fontsize',14,'XLim',[0 1.0]);
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[0.0 0.5 1.0]);
set(gca,'XTick',[0.0 0.5 1.0]);
%set(gca,'XTickLabel',[0.0 0.25 0.5 0.75 1.0]);
%set(gca,'XTick',[0.0 0.25 0.5 0.75 1.0]);
hold on; box on;

hsum1 = 0;
hsum2 = 0;
hsum3 = 0;
nl = (length(M(1).m)-3)/npl
for i = 1:nl

    c(i,1) = hpds(((i-1)*npl)+2,1);
    c(i,2) = hpds(((i-1)*npl)+2,2);
    rho(i,1) = hpds(((i-1)*npl)+3,1);
    rho(i,2) = hpds(((i-1)*npl)+3,2);
    alpha(i,1) = hpds(((i-1)*npl)+4,1);
    alpha(i,2) = hpds(((i-1)*npl)+4,2);
    h(i,1) = hpds(((i-1)*npl)+1,1)+hsum1;
    h(i,2) = hpds(((i-1)*npl)+1,2)+hsum2;
    hsum1 = h(i,1);
    hsum2 = h(i,2);
    hh(i,1) = M(1).m(((i-1)*npl)+1);
    hh(i,2) = M(1).m(((i-1)*npl)+1)+hsum3;
    hsum3 = hh(i,2);
    hh(i,3) = hpds(((i-1)*npl)+1,1);
    hh(i,4) = hpds(((i-1)*npl)+1,2);

end;

c(nl+1,1) = hpds((nl*npl)+1,1);
c(nl+1,2) = hpds((nl*npl)+1,2);
rho(nl+1,1) = hpds((nl*npl)+2,1);
rho(nl+1,2) = hpds((nl*npl)+2,2);
alpha(nl+1,1) = hpds((nl*npl)+3,1);
alpha(nl+1,2) = hpds((nl*npl)+3,2);

if(c(2,1)-c(1,1) > 0)
    hc(1,1) = hh(1,4);
else
    hc(1,1) = hh(1,3);
end
if(c(2,2)-c(1,2) > 0)
    hc(1,2) = hh(1,3);
else
    hc(1,2) = hh(1,4);
end

for i = 2:length(c)-1
    if(c(i+1,1)-c(i,1) > 0)
        hc(i,1) = hh(i,2)+(hh(i,4)-hh(i,1));
    else
        hc(i,1) = hh(i,2)-(hh(i,1)-hh(i,3));
    end
    if(c(i+1,2)-c(i,2) > 0)
        hc(i,2) = hh(i,2)-(hh(i,1)-hh(i,3));
    else
        hc(i,2) = hh(i,2)+(hh(i,4)-hh(i,1));
    end
end
hc(length(c),1) = hmax;
hc(length(c),2) = hmax;

if(rho(2,1)-rho(1,1) > 0)
    hr(1,1) = hh(1,4);
else
    hr(1,1) = hh(1,3);
end
if(rho(2,2)-rho(1,2) > 0)
    hr(1,2) = hh(1,3);
else
    hr(1,2) = hh(1,4);
end

for i = 2:length(c)-1
    if(rho(i+1,1)-rho(i,1) > 0)
        hr(i,1) = hh(i,2)+(hh(i,4)-hh(i,1));
    else
        hr(i,1) = hh(i,2)-(hh(i,1)-hh(i,3));
    end
    if(rho(i+1,2)-rho(i,2) > 0)
        hr(i,2) = hh(i,2)-(hh(i,1)-hh(i,3));
    else
        hr(i,2) = hh(i,2)+(hh(i,4)-hh(i,1));
    end
end
hr(length(c),1) = hmax;
hr(length(c),2) = hmax;

if(ishade == 0)
  subplot(h1);hold on; box on;
  h = 0;h_old = 0;
  v = 0;v_old = 0;
  for il = 1:nl

      h = hc(il,1);
      v = c(il,1);
      if(il > 1)
          plot([v_old v],[h_old h_old],col(3,:));
      end
      plot([v v],[h_old h],col(3,:));
      h_old = h;
      v_old = v;

  end
  v = c(end,1);
  plot([v v],[h_old hc(end,1)],col(3,:));
  plot([v_old v],[h_old h_old],col(3,:));
  h = 0;h_old = 0;
  v = 0;v_old = 0;
  for il = 1:nl

      h = hc(il,2);
      v = c(il,2);
      if(il > 1)
          plot([v_old v],[h_old h_old],col(3,:));
      end
      plot([v v],[h_old h],col(3,:));
      h_old = h;
      v_old = v;

  end
  v = c(end,2);
  plot([v v],[h_old hc(end,2)],col(3,:));
  plot([v_old v],[h_old h_old],col(3,:));

  subplot(h2);hold on; box on;
  h = 0;h_old = 0;
  r = 0;r_old = 0;
  for il = 1:nl

      h = hr(il,1);
      r = rho(il,1);
      if(il > 1)
          plot([r_old r],[h_old h_old],col(3,:));
      end
      plot([r r],[h_old h],col(3,:));
      h_old = h;
      r_old = r;

  end
  r = rho(end,1);
  plot([r r],[h_old hr(end,1)],col(3,:));
  plot([r_old r],[h_old h_old],col(3,:));
  h = 0;h_old = 0;
  r = 0;r_old = 0;
  for il = 1:nl

      h = hr(il,2);
      r = rho(il,2);
      if(il > 1)
          plot([r_old r],[h_old h_old],col(3,:));
      end
      plot([r r],[h_old h],col(3,:));
      h_old = h;
      r_old = r;

  end
  r = rho(end,2);
  plot([r r],[h_old hr(end,2)],col(3,:));
  plot([r_old r],[h_old h_old],col(3,:));

else

poly(1,1) = c(1,1);
poly(1,2) = 0.001;
poly(2,1) = c (1,1);
poly(2,2) = hc(1,1);
j = 3;
for i = 2:nl

    poly(j,1) = c(i,1);
    poly(j,2) = hc(i-1,1);
    poly(j+1,1) = c(i,1);
    poly(j+1,2) = hc(i,1);
    j = j + 2;

end;
poly(j,1) = c(end,1);
poly(j,2) = hc(end-1,1);
poly(j+1,1) = c (end,1);
poly(j+1,2) = hc(end,1)-0.01;


poly2(1,1) = c (1,2);
poly2(1,2) = 0.001;
poly2(2,1) = c (1,2);
poly2(2,2) = hc(1,2);

j = 3;
for i = 2:nl

    poly2(j,1) = c(i,2);
    poly2(j,2) = hc(i-1,2);
    poly2(j+1,1) = c(i,2);
    poly2(j+1,2) = hc(i,2);
    j = j + 2;

end;
poly2(j,1) = c(end,2);
poly2(j,2) = hc(end-1,2);
poly2(j+1,1) = c (end,2);
poly2(j+1,2) = hc(end,2)-0.01;

poly = [poly; flipud(poly2)];
subplot(h1);hold on; box on;
fill(poly(:,1),poly(:,2),[.8 .8 .8],'LineStyle','none');

clear poly;
clear poly2;

poly(1,1) = rho(1,1);
poly(1,2) = 0.001;
poly(2,1) = rho(1,1);
poly(2,2) = hr(1,1);
j = 3;
for i = 2:nl

    poly(j,1) = rho(i,1);
    poly(j,2) = hr(i-1,1);
    poly(j+1,1) = rho(i,1);
    poly(j+1,2) = hr(i,1);
    j = j + 2;

end;
poly(j,1) = rho(end,1);
poly(j,2) = hr(end-1,1);
poly(j+1,1) = rho(end,1);
poly(j+1,2) = hr(end,1)-0.01;


poly2(1,1) = rho(1,2);
poly2(1,2) = 0.001;
poly2(2,1) = rho(1,2);
poly2(2,2) = hr(1,2);

j = 3;
for i = 2:nl

    poly2(j,1) = rho(i,2);
    poly2(j,2) = hr(i-1,2);
    poly2(j+1,1) = rho(i,2);
    poly2(j+1,2) = hr(i,2);
    j = j + 2;

end;
poly2(j,1) = rho(end,2);
poly2(j,2) = hr(end-1,2);
poly2(j+1,1) = rho(end,2);
poly2(j+1,2) = hr(end,2)-0.01;

poly = [poly; flipud(poly2)];
subplot(h2);hold on; box on;
fill(poly(:,1),poly(:,2),[.8 .8 .8],'LineStyle','none');

clear poly;
clear poly2;

% Alpha
if(alpha(2,1)-alpha(1,1) > 0)
    ha(1,1) = hh(1,4);
else
    ha(1,1) = hh(1,3);
end
if(alpha(2,2)-alpha(1,2) > 0)
    ha(1,2) = hh(1,3);
else
    ha(1,2) = hh(1,4);
end

for i = 2:length(c)-1
    if(alpha(i+1,1)-alpha(i,1) > 0)
        ha(i,1) = hh(i,2)+(hh(i,4)-hh(i,1));
    else
        ha(i,1) = hh(i,2)-(hh(i,1)-hh(i,3));
    end
    if(alpha(i+1,2)-alpha(i,2) > 0)
        ha(i,2) = hh(i,2)-(hh(i,1)-hh(i,3));
    else
        ha(i,2) = hh(i,2)+(hh(i,4)-hh(i,1));
    end
end
ha(length(c),1) = hmax;
ha(length(c),2) = hmax;

poly(1,1) = alpha(1,1);
poly(1,2) = 0.001;
poly(2,1) = alpha(1,1);
poly(2,2) = ha(1,1);
j = 3;
for i = 2:nl

    poly(j,1) = alpha(i,1);
    poly(j,2) = ha(i-1,1);
    poly(j+1,1) = alpha(i,1);
    poly(j+1,2) = ha(i,1);
    j = j + 2;

end;
poly(j,1) = alpha(end,1);
poly(j,2) = ha(end-1,1);
poly(j+1,1) = alpha(end,1);
poly(j+1,2) = ha(end,1)-0.01;


poly2(1,1) = alpha(1,2);
poly2(1,2) = 0.001;
poly2(2,1) = alpha(1,2);
poly2(2,2) = ha(1,2);
j = 3;
for i = 2:nl

    poly2(j,1) = alpha(i,2);
    poly2(j,2) = ha(i-1,2);
    poly2(j+1,1) = alpha(i,2);
    poly2(j+1,2) = ha(i,2);
    j = j + 2;

end;
poly2(j,1) = alpha(end,2);
poly2(j,2) = ha(end-1,2);
poly2(j+1,1) = alpha(end,2);
poly2(j+1,2) = ha(end,2)-0.01;

poly = [poly; flipud(poly2)];
subplot(h3);
hold on; box on;
fill(poly(:,1),poly(:,2),[.8 .8 .8],'LineStyle','none');
figure(1);hold on; box on;

save bla poly poly2 alpha ha hh;

end;
for i = 1:size(M,2)

  nl = (length(M(i).m)-3)/npl
% Velocity
  subplot(h1);box on;
  set(gca,'YDir','reverse','YLim',[0 hmax]);
  h = 0;h_old = 0;
  v = 0;v_old = 0;
  for il = 1:nl

      h = h + M(i).m(((il - 1) * npl) + 1);
      v = M(i).m(((il - 1) * npl) + 2);
      if(il > 1)
          plot([v_old v],[h_old h_old],col(i,:),'Linewidth',2);
      end
      plot([v v],[h_old h],col(i,:),'Linewidth',2);
      h_old = h;
      v_old = v;

  end
  v = M(i).m(end - 2);
  plot([v v],[h_old hmax],col(i,:),'Linewidth',2);
  handle(i) = plot([v_old v],[h_old h_old],col(i,:),'Linewidth',2);
  xlabel('Velocity (m/s)');
  ylabel('Depth (m)');
  ylims = get(gca,'Ylim');



% Density
  subplot(h2);box on;
  set(gca,'YDir','reverse','YLim',[0 hmax]);

  h = 0;h_old = 0;
  r = 0;r_old = 0;
  for il = 1:nl

      h = h + M(i).m(((il - 1) * npl) + 1);
      r = M(i).m(((il - 1) * npl) + 3);
      if(il > 1)
          plot([r_old r],[h_old h_old],col(i,:),'Linewidth',2);
      end
      plot([r r],[h_old h],col(i,:),'Linewidth',2);
      h_old = h;
      r_old = r;

  end
  r = M(i).m(end - 1);
  plot([r r],[h_old hmax],col(i,:),'Linewidth',2);
  plot([r_old r],[h_old h_old],col(i,:),'Linewidth',2);
  xlabel('Density (g/ccm)');
%  subplot(h1);box on;
%  xlims = get(gca,'XLim');
%  plot([xlims(1) xlims(2)],[3.5 3.5],':k');
%  plot([xlims(1) xlims(2)],[10.0 10.0],':k');
%  subplot(h2);box on;
%  xlims = get(gca,'XLim');
%  plot([xlims(1) xlims(2)],[3.5 3.5],':k');
%  plot([xlims(1) xlims(2)],[10.0 10.0],':k');

% Attenuation
  subplot(h3);hold on; box on;
  set(gca,'YDir','reverse','YLim',[0 hmax]);

  h = 0;h_old = 0;
  a = 0;a_old = 0;
  for il = 1:nl

      h = h + M(i).m(((il - 1) * npl) + 1);
      a = M(i).m(((il - 1) * npl) + 4);
      if(il > 1)
          plot([a_old a],[h_old h_old],col(i,:),'Linewidth',2);
      end
      plot([a a],[h_old h],col(i,:),'Linewidth',2);
      h_old = h;
      a_old = a;

  end
  a = M(i).m(end);
  plot([a a],[h_old hmax],col(i,:),'Linewidth',2);
  plot([a_old a],[h_old h_old],col(i,:),'Linewidth',2);
  xlabel('Attenuation (dB/L)');

end;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'--k','Linewidth',2);
    plot(c2(:,2),c2(:,1),'--k','Linewidth',2);

    subplot(h2);
    plot(r1(:,2),r1(:,1),'--k','Linewidth',2);
    plot(r2(:,2),r2(:,1),'--k','Linewidth',2);
end;

subplot(h1); box on;
subplot(h2); box on;
subplot(h3); box on;
saveas(hf1,plotfile1,'epsc2');

%subplot(h1);hold on; box on;
%legend([handle(1) handle(2) handle(3)],'true model','sph. inv.','pl. inv.')

return;
