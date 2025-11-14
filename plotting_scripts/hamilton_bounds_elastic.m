function [log_area]=hamilton_bounds();

minlim = [1450. 1.2];
maxlim = [6000. 2.6];

rho = [minlim(2):(maxlim(2)-minlim(2))/5000.:maxlim(2)];

cl=(1.54-0.907*rho+0.3695*rho.^1.88)*1.5004*1000.;
ch=(1.62-0.907*rho+0.3695*rho.^2.05)*1.5014*1000.;

cl(find(cl<minlim(1))) = minlim(1);
ch(find(ch>maxlim(1))) = maxlim(1);
ch(find(rho>2.000)) = maxlim(1);

log_area = sum((ch-cl).*(rho(2)-rho(1)))
tot_area = prod(maxlim-minlim)
ratio = log_area/tot_area

figure();box on;hold on;
plot(rho,cl);plot(rho,ch);
set(gca,'XLim',[minlim(2) maxlim(2)]);
set(gca,'YLim',[minlim(1) maxlim(1)]);

return;
