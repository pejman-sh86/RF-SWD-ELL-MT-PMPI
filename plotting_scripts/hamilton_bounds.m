%function [log_area]=hamilton_bounds();

minlim = [1450. 1.2];
maxlim = [2500. 2.2];

rho = [minlim(2):(maxlim(2)-minlim(2))/1000.:maxlim(2)];

cl=(1.515-0.89*rho+0.3695*rho.^1.87)*1.5004*1000.;
ch=(1.58-0.907*rho+0.3695*rho.^2.05)*1.5014*1000.;

%cl(find(cl<minlim(1))) = minlim(1);
%ch(find(ch>maxlim(1))) = maxlim(1);

log_area = sum((ch-cl).*(rho(2)-rho(1)))
tot_area = prod(maxlim-minlim)
ratio = log_area/tot_area

figure();box on;hold on;
plot(rho,cl,'LineWidth',2);
plot(rho,ch,'LineWidth',2);
set(gca,'XLim',[minlim(2) maxlim(2)]);
set(gca,'YLim',[minlim(1) maxlim(1)]);
xlabel('Density (g/ccm)');
ylabel('Velocity (m/s)');
%return;
