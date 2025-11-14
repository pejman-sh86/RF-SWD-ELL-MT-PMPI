function [] = pick_tt();

format long g;
i_old = 1;
ntrace = 70;

set(0, 'DefaultFigurePaperPosition', [0 0 8 4]);

%data = 'pro_sim_A.mat';
%plotfile = 'pro_sim_A.eps';
data = 'x.mat';
plotfile = 'x.eps';
load(data);
if(i_old == 1)
%    replica = 'replica_true.mat';
%    load(replica);
    old_picked = 'tt_rep.mat';
%    old_picked = 'tmp.mat';
    load(old_picked);
end
if (min(x.r) < 0); x.r  = -1 * x.r;end; %  shot ranges
%x.r = x.r(1:32)

nlay = 7;
tsq_lim=[0.05  0.4];
%t_lim=[0.128 0.175];
t_lim=[0.05 0.4];
col = {'b' 'k' ':k' '--k' '--r' ':k'};
col = char(col);
offset = 0;
offset2 = 0;

% -----------------------------------------------------------------------
%  Hyperbolas
% -----------------------------------------------------------------------

% Site 13:
itbb = [5200 5440 5500 6100 6330 6990];
%itbb = [2841 2891 2979 3015 3219 3258 3368];
%itbb = [3080 3133 3219 3274 3503 3565 3767];
%itbb = [3944 4000 4400 5000 6330];

t0 = x.t_ax(itbb)
cnmo = 1513.26*ones(1,length(itbb));
cnmo = [1513.26 1513.26 1513.26 1513.26 1513.26 1513.26]

% -----------------------------------------------------------------------
%  Plot data
% -----------------------------------------------------------------------
fig = 1;

%fact = 4e3
fact = 1e2

figh = figure(fig+1);
hold on;box on;
i = find(t_lim(1) <= x.t_ax & x.t_ax <= t_lim(2));
x1 = x.t_s(:, i);
for i1 = 1:length(x.r)
     x1(i1, :) = fact*x1(i1, :)+x.r(i1);
end
plot(x1, x.t_ax(i), col(1,:));

axis([x.r(1), x.r(70), t_lim(1), t_lim(2)])
%axis([10, 140, 0.128, 0.19])
xlabel('range (m)','FontSize',14)
ylabel('TWT (s)','FontSize',14)
set(gca,'YDir','reverse')
set(gca,'FontSize',14)
title ('Traces');

%-----------------------------------------------------------------------
%
%-----------------------------------------------------------------------
%t0 = 0.15667;
%v = 1900.;

%tt(length(itbb)+1,:) = x.r(1:ntrace);

if(i_old == 1)

  for hyper = 1:length(itbb)
     plot(r(1:ntrace), tt(hyper,1:ntrace), '--k')
  end

end

%tsq_lim=[0.07 0.14];
%tsq2 = sqrt(rep2.^2-repmat(r.^2/1511.^2,nlay,1));
%[figh_red] = reduce_sqrt(x,r,tsq_lim,4e2,tsq2,itbb,3);

%saveas(figh,plotfile,'epsc2');
return;

itlen = 40;
icorlen = 120;

for hyper = 1:length(itbb)

    load(data);
    if (min(x.r) < 0); x.r  = -1 * x.r;end; %  shot ranges
    itb = itbb(hyper)
    if(i_old == 0)
    dt(1) = x.t_ax(itb);
    itbdiff = 0;
    for ir = 2:ntrace
  
        if(x.r(ir)-x.r(ir-1) < 10)
            itb2 = itb + floor(3/4*itbdiff);
        else
            itb2 = itb + floor(2*itbdiff);
        end
%        [xcy,lags] = xcorr(x.t_s(ir-1,itb:itb+itlen),...
%                     x.t_s(ir,itb2:itb2+icorlen));
%        [val,idx] = max(xcy);
%        lags(idx);
%        itbold = itb2;
%        itb = itb2-lags(idx);
%        dt(ir) = x.t_ax(itb);
        dt(ir) = dt(ir-1);

    end
    else
       dt = tt(hyper,:);
    end

dt = t0(hyper)*(1 + (abs(x.r(1:ntrace))/(cnmo(hyper)*t0(hyper))).^2).^0.5;

handle = plot(x.r(1:length(dt)), dt, 'k')

n=0;
for n=1:ntrace*10
    zoom on;
    while ~waitforbuttonpress
    end
    [x,y,b]=ginput(1);
    L = handle;
    X=get(L,'xdata');
    Y=get(L,'ydata');
    [v,j]=min((X-x).^2+(Y-y).^2);% find closest point to clicked location
    if (b == 3)
        break;
    else
        [x,y,b]=ginput(1);
        R = [x y]; % keep the point location
        R
    end
    Xf=[X(1:j-1) R(1) X(j+1:end)]; % insert the new location
    Yf=[Y(1:j-1) R(2) Y(j+1:end)];
    set(L,'xdata',Xf,'ydata',Yf);
end
Yf = get(L,'ydata');

tt(hyper,:) = Yf;
%save tt_tmp tt;
end % End of interfaces loop

%save picked_tt2 tt;
return;
%=============================================================================

function [figh] = reduce_sqrt(x,r, tsq_lim, fact, tsq2, cnmo, fig);

figh = figure(fig);
for i = 1:length(x.r)
     tsq = sqrt(x.t_ax.^2-x.r(i)^2/x.c^2);
     i1 = find(tsq_lim(1) <= tsq & tsq <= tsq_lim(2));
     line('Color', 'b', ...
          'XData', fact*x.t_s(i, i1)+x.r(i), ...
          'YData', tsq(i1), 'LineWidth',.1)
end
xlabel('Range (m)'); ylabel('Reduced Time (s)');
axis([min(x.r) 350 tsq_lim]); box on;
set(gca,'YDir','reverse'); title ('Traces');
hold on;

for iw = 1:length(cnmo)

  plot(r,tsq2(iw,:),'k');

end

return;

% =============================================================================
%  This is the end my fiend...
%  EOF
