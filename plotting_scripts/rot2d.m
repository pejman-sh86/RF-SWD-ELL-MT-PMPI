function [xr,yr] = rot2d(x,y,angle);
%%
%% Counter-clockwise rotation of pair (x,y)
%%

angr = -pi/180.*angle;
xr = x * cos(angr) - y * sin(angr);
yr = x * sin(angr) + y * cos(angr);

return;
