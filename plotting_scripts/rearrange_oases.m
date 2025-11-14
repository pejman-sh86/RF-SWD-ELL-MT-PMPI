%
% Rearrange Holland data simulated by OASES
%
function [] = rearrange_oases();

dataload  = 'trex_sim_boomer.mat';
datawrite = 'x_trex_sim_boomer.mat';
load(dataload);

x.r = double(RANGE');
%
% OASES outputs negative pressure!
%
x.t_s = -1 * double(DATA');
for i_step = 1:length(x.t_s)

    x.t_ax(i_step) = (i_step-1) * double(DELTAT);

end
x.c = 1525.;

save(datawrite,'x');

return;

