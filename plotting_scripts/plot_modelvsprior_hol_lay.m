function [] = plot_modelvsprior_hol_lay();

set(0, 'DefaultFigurePaperPosition', [0 0 7 4]);
datafile1 = 'x_jan_1.mat';
datafile2 = 'x_jan_2.mat';
datafile4 = 'x_jan_4.mat';

ishade = 1;

%
% MAP model
%
M(1).m = [1.6756, 1473.8019, 1.6601, 0.1500, ...
          1.8600, 1555.2124, 1.8313, 0.0100, ...
          1.1815, 1540.6675, 1.7023, 0.1444, ...
          3.84065413,1621.49719,1.95841837,0.149963692,...
          1.33724737,1681.96155,1.82700408,0.149790883,...
          3.15582967,1600.0011, 1.70010376,0.148903981,...
                     1676.42761,1.99389005,1.050061919E-2];

hpd = zeros(length(M(1).m),2);
load(datafile1);
hpd(1:4,1) = F(1).minlim(1:end-3);
hpd(1:4,2) = F(1).maxlim(1:end-3);

load(datafile2);
hpd(5:12,1) = F(1).minlim(1:end-3);
hpd(5:12,2) = F(1).maxlim(1:end-3);

load(datafile4);
hpd(13:27,1) = F(1).minlim;
hpd(13:27,2) = F(1).maxlim;

hpd

npl = 4;
length(M(1).m)
nl = (length(M(1).m)-3)/npl

col = {'-k' '--k' ':k' '--k' '--r' ':k'};
col = char(col);

figure(1);hold on; box on;
h1 = subplot('Position',[.08 .15 .44 .8]);
set(gca,'Fontsize',14,'XLim',[1440 1720]);
set(gca,'XTickLabel',[1500 1550 1600 1650 1700]);
set(gca,'XTick',[1500 1550 1600 1650 1700]);
hold on; box on;
h2 = subplot('Position',[.54 .15 .44 .8]);
set(gca,'Fontsize',14,'XLim',[1.3 2.45]);
set(gca,'YTickLabel',[]);
hold on; box on;
hsum1 = 0;
hsum2 = 0;
for i = 1:nl

    c(i,1) = hpd(((i-1)*npl)+2,1);
    c(i,2) = hpd(((i-1)*npl)+2,2);
    rho(i,1) = hpd(((i-1)*npl)+3,1);
    rho(i,2) = hpd(((i-1)*npl)+3,2);
    h(i,1) = hpd(((i-1)*npl)+1,1)+hsum1;
    h(i,2) = hpd(((i-1)*npl)+1,2)+hsum2;
    hsum1 = h(i,1);
    hsum2 = h(i,2);

end;

hmax = h(end,2) + h(end,2)/10
c(nl+1,1) = hpd((nl*npl)+1,1);
c(nl+1,2) = hpd((nl*npl)+1,2);
rho(nl+1,1) = hpd((nl*npl)+2,1);
rho(nl+1,2) = hpd((nl*npl)+2,2);

for i = 1:length(c)-1
    if(c(i+1)-c(i) > 0)
        hc(i,1) = h(i,2);
        hc(i,2) = h(i,1);
    else
        hc(i,1) = h(i,1);
        hc(i,2) = h(i,2);
    end
    if(rho(i+1)-rho(i) > 0)
        hr(i,1) = h(i,2);
        hr(i,2) = h(i,1);
    else
        hr(i,1) = h(i,1);
        hr(i,2) = h(i,2);
    end
end

hr(length(c),1) = hmax;
hr(length(c),2) = hmax;
hc(length(c),1) = hmax;
hc(length(c),2) = hmax;

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
poly(1,2) = 0.03;
poly(2,1) = c (1,1);
poly(2,2) = hc(1,1);
j = 3
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
poly2(1,2) = 0.03;
poly2(2,1) = c (1,2);
poly2(2,2) = hc(1,2);

j = 3
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
poly(1,2) = 0.03;
poly(2,1) = rho(1,1);
poly(2,2) = hr(1,1);
j = 3
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
poly2(1,2) = 0.03;
poly2(2,1) = rho(1,2);
poly2(2,2) = hr(1,2);

j = 3
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



end;
for i = 1:size(M,2)

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

end;

return;
