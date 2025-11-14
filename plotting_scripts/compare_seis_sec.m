function [] = compare_seis_sec();

trial = 1;

set(0, 'DefaultFigurePaperPosition', [0 0 7 9.5]);

load('pro_jan_real.mat');
if (min(x.r) < 0)
    x.r  = -1 * x.r;      %  shot ranges
end
xv = x;
load('x.mat');
if (min(x.r) < 0)
    x.r  = -1 * x.r;      %  shot ranges
end

tsq_lim=[0.00  1.];
t_lim=[0.05 0.55];
col = {'b' 'k' ':k' '--k' '--r' ':k'};
col = char(col);
offset = 0;
offset2 = 0;

%-----------------------------------------------------------------------
%  Plot data
%-----------------------------------------------------------------------

fig = 1;
[figh_tr] = raw_plot(x, t_lim, 4e2,fig,col(1,:),0);
%fig = fig+1;
[figh_tr] = raw_plot(xv, t_lim, 4e2,fig,col(2,:),offset);
%saveas(figh_tr,tr_file,'epsc2');
fig = fig+1;

if(trial == 1)
    [figh_red] = reduce_sqrt(x,tsq_lim,2e2,fig,col(1,:),0);
    [figh_red] = reduce_sqrt(xv,tsq_lim,2e2,fig,col(2,:),offset2);
%    saveas(figh_red,red_file,'epsc2');
    fig = fig+1;
    return;
end

%=============================================================================

function [figh] = raw_plot(x, t_lim, fact, fig,col,offset)

col
figh = figure(fig);
hold on;box on;
i = find(t_lim(1) <= x.t_ax & x.t_ax <= t_lim(2));
x1 = x.t_s(:, i);
for i1 = 1:length(x.r)
     x1(i1, :) = fact*x1(i1, :)+x.r(i1);
end
plot(x1, x.t_ax(i)+offset, col)
axis([0, max(x.r), t_lim(1), t_lim(2)])
%axis([10, 140, 0.128, 0.19])
xlabel('range (m)')
ylabel('TWT (s)')
set(gca,'YDir','reverse')
title ('Traces');

saveas(gca,'plot.eps','epsc2');

return;

%=============================================================================

function [figh] = reduce_sqrt(x, tsq_lim, fact, fig,col,offset,stretch);

col
x.t_ax = x.t_ax+offset;
figh = figure(fig);
for i = 1:length(x.r)
     tsq = sqrt(x.t_ax.^2-x.r(i)^2/x.c^2);
     i1 = find(tsq_lim(1) <= tsq & tsq <= tsq_lim(2));
     line('XData', fact*x.t_s(i, i1)+x.r(i), ...
          'YData', tsq(i1), 'LineWidth',.1,'Color',col)
end
xlabel('Range (m)'); ylabel('Reduced Time (s)');
axis([min(x.r) max(x.r) tsq_lim]); box on;
set(gca,'YDir','reverse'); title ('Traces');
hold on;

return;

%=============================================================================
%  This is the end my fiend...
%  EOF
%
