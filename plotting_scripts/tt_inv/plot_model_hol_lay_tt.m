function [] = plot_model_hol_lay();

set(0, 'DefaultFigurePaperPosition', [0 0 5 7]);
%set(0, 'DefaultFigurePaperPosition', [0 0 9 5.14]);
%set(0, 'DefaultFigurePaperPosition', [0 0 14 8]);
hpdfile = 'picked_tt_int_covstat_hpds.txt';
mapfile = 'picked_tt_int_covstat_map.txt';
meanfile = 'picked_tt_int_covstat_mean.txt';
stdfile = 'picked_tt_int_covstat_std.txt';
plotfile1 = 'picked_tt_int_covstat_results.eps';
hpd = load(hpdfile);

load core;

ishade = 1;
i_nf = 0;

%M(1).m = load(mapfile);
M(1).m = load(meanfile);
M(2).m = load(stdfile);
%
% true model
%

npl = 2;
nl = (length(M(1).m)-7)/npl

%M(2).m = hpd(:,1);
%M(3).m = hpd(:,2);
%M(3).m([1 5 9 13 17 21]) = M(2).m([1 5 9 13 17 21]);
%M(4).m([1 5 9 13 17 21]) = M(2).m([1 5 9 13 17 21]);

col = {'-k' '--k' ':k' '--k' '--r' ':k'};
col = char(col);

figure(1);hold on; box on;
set(gca,'Fontsize',14,'XLim',[1450 1750]);
set(gca,'XTickLabel',[1500 1550 1600 1650 1700]);
set(gca,'XTick',[1500 1550 1600 1650 1700]);
hold on; box on;
hsum1 = 0;
hsum2 = 0;
for i = 1:nl

    if(i_nf == 1)
        c(i,1) = hpd(((i-1)*npl)+2,1);
        c(i,2) = hpd(((i-1)*npl)+2,2);
        h(i,1) = hpd(((i-1)*npl)+1,1)+hsum1;
        h(i,2) = hpd(((i-1)*npl)+1,2)+hsum2;
        hsum1 = h(i,1);
        hsum2 = h(i,2);
    else
       c(i,1) = M(1).m(((i-1)*npl)+2)-M(2).m(((i-1)*npl)+2);
       c(i,2) = M(1).m(((i-1)*npl)+2)+M(2).m(((i-1)*npl)+2);
       h(i,1) = M(1).m(((i-1)*npl)+1)-M(2).m(((i-1)*npl)+1)+hsum1;
       h(i,2) = M(1).m(((i-1)*npl)+1)+M(2).m(((i-1)*npl)+1)+hsum2;
       hsum1 = h(i,1);
       hsum2 = h(i,2);
    end

end;
%hmax = h(end,2) + h(end,2)/10
hmax = 23;

%c(nl+1,1) = hpd((nl*npl)+1,1);
%c(nl+1,2) = hpd((nl*npl)+1,2);

for i = 1:length(c)-1
    if(c(i+1)-c(i) > 0)
        hc(i,1) = h(i,2);
        hc(i,2) = h(i,1);
    else
        hc(i,1) = h(i,1);
        hc(i,2) = h(i,2);
    end
end

hc(length(c),1) = hmax;
hc(length(c),2) = hmax;

if(ishade == 0)
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
fill(poly(:,1),poly(:,2),[.8 .8 .8],'LineStyle','none');

clear poly;
clear poly2;

end;
for i = 1:size(M,2)

% Velocity
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
      h_old = h
      v_old = v

  end
  v = M(i).m(end - 2);
  plot([v],[h_old],col(i,:),'Linewidth',2);
%  handle(i) = plot([v_old v],[h_old h_old],col(i,:),'Linewidth',2);
  xlabel('Velocity (m/s)');
  ylabel('Depth (m)');
  ylims = get(gca,'Ylim');

end;
%plot(c1(:,2),c1(:,1),'r','Linewidth',2);
%plot(c2(:,2),c2(:,1),'r','Linewidth',2);

saveas(gca,plotfile1,'epsc2');

%legend([handle(1) handle(2) handle(3)],'true model','sph. inv.','pl. inv.')

return;
