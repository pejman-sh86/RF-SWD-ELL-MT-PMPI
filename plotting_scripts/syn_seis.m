function [] = syn_seis();
%
% Simple synthetic seismograms for 2 halfspaces with
% constant velocity
%
load pro_sim_E;
load source.dat;

x.r = -1 * x.r;
nr = length(x.r);
nt = length(x.t_ax);
zs   =   0.35;
zr   = 122.00;
zw   = 150.00;
c    = 1510;
r1   = 1.1;
r2   = 1.5;

h  = zr - zs + 2*(zw-zr)
h2 = zr - zs
rho = (r2-r1)/(r2+r1)

x.t_s = x.t_s.*0;

for ir = 1:nr

    len_r(ir) = (h^2+x.r(ir)^2)^(.5);
    len_d(ir) = (h2^2 + x.r(ir)^2)^(.5);
    td(ir) = len_d(ir)/c;
    tr(ir) = len_r(ir)/c;
    ad(ir) = 1/len_d(ir);
    ar(ir) = 1 * rho * 1/len_r(ir);
%    ad(ir) = 1;
%    ar(ir) = rho;

    td(ir)-x.t_ax(2);
    td(ir)+x.t_ax(2);

    i_d = find(td(ir)-x.t_ax(2) < x.t_ax & x.t_ax <= td(ir)+x.t_ax(2));
    i_r = find(tr(ir)-x.t_ax(2) < x.t_ax & x.t_ax <= tr(ir)+x.t_ax(2));
    x.t_s(ir,i_d(1)) = ad(ir);
    x.t_s(ir,i_r(1)) = ar(ir);

end

for ir = 1:nr

    x.tmp(ir,:) = conv(x.t_s(ir,:),source(:,2));
    x.t_s(ir,:) = x.tmp(ir,1:length(x.t_ax));

end

tnmo = [td;tr];
t_lim=[0.05 0.55];
[figh_tr] = raw_plot(x, t_lim, 2e2,tnmo,[1510 1510],1);

save pro_sim_F x;

return;


%=============================================================================

function [figh] = raw_plot(x, t_lim, fact,tnmo,cnmo, fig)

figh = figure(fig);
hold on;box on;
i = find(t_lim(1) <= x.t_ax & x.t_ax <= t_lim(2));
x1 = x.t_s(:, i);
for i1 = 1:length(x.r)
     x1(i1, :) = fact*x1(i1, :)+x.r(i1);
end
plot(x1, x.t_ax(i), 'b')
axis([0, max(x.r), x.t_ax(i(1)), x.t_ax(i(end))])
xlabel('range (m)')
ylabel('Time re trigger (s)')
set(gca,'YDir','reverse')
title ('Traces');


for iw = 1:length(cnmo)

  plot(x.r,tnmo(iw,:),'k');

end

return;

