function [] = plot_model_hol_lay();

set(0, 'DefaultFigurePaperPosition', [0 0 11 8]);
%set(0, 'DefaultFigurePaperPosition', [0 0 14 8]);

corefile  = 'core13.mat';
plotfile1 = 'tmp.eps';

icore = 1;

if(icore == 1)
    load(corefile);
end;

hmax = 3;

%M(1).m = load(mapfile);
M(1).m = [0.200000000001456        1500.19377054163        1.56512831418132       1.000000001521314E-002  0.453462672231634        1560.78279238104        1.60000000000370       1.000000292928564E-002   1500.00007320315        1.60000000001249       1.000000017811677E-002];
M(2).m = [0.102550920377054        1450.00000458916        1.50000000002236       0.339721215185228       0.302552797229459        1552.95836415918        1.66929471457838       1.000000047683751E-002  0.537287856134042        1576.83556887956        1.79561666282983       0.459446562467888        1565.87444588445        1.98771616929437       1.000000085905595E-002];
         
npl = 4;

col = {'-k' '-r' '-b' '--k' '--r' ':k'};
col = char(col);

nx = 3
ny = 1

xim = 0.03;
yim = 0.28/ny;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 20;
nbin2 = 10;

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',8,'Height',2*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

hf1=figure(1);hold on; box on;
h1 = subplot('Position',[loc(1,1) loc(2,1) spw sph]);
set(gca,'Fontsize',14,'XLim',[1450 1700]);
set(gca,'XTickLabel',[1500 1550 1600 1650 1700]);
set(gca,'XTick',[1500 1550 1600 1650 1700]);
hold on; box on;
h2 = subplot('Position',[loc(1,2) loc(2,2) spw sph]);
set(gca,'Fontsize',14,'XLim',[1.3 2.2]);
set(gca,'YTickLabel',[]);
hold on; box on;
h3 = subplot('Position',[loc(1,3) loc(2,3) spw sph]);
hold on; box on;

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
  set(gca,'Fontsize',14);
  set(gca,'Fontsize',14,'XLim',[0 1.0]);
  set(gca,'YTickLabel',[]);
  set(gca,'XTickLabel',[0.0 0.25 0.5 0.75 1.0]);
  set(gca,'XTick',[0.0 0.25 0.5 0.75 1.0]);

end;
if(icore == 1)
    subplot(h1);
    plot(c1(:,2),c1(:,1),'r','Linewidth',2);
%    plot(c2(:,2),c2(:,1),'r','Linewidth',2);

    subplot(h2);
    plot(r1(:,2),r1(:,1),'r','Linewidth',2);
%    plot(r2(:,2),r2(:,1),'r','Linewidth',2);
end;

subplot(h1); box on;
subplot(h2); box on;
subplot(h3); box on;
%saveas(hf1,plotfile1,'epsc2');
exportfig(hf1,plotfile1,opts);

%subplot(h1);hold on; box on;
%legend([handle(1) handle(2) handle(3)],'true model','sph. inv.','pl. inv.')

return;
