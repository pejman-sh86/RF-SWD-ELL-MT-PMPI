function []=cmt_make_bounds();

dip  = 20.; % dip angle of subducting plate in degree
zmin = 10.; % min depth due to GF data base
latitude = 52.55;
width = 1.0;

offsetkm = zmin/tan(dip*(pi/180.));
offsetdeg = offsetkm/(cos(latitude*(pi/180.))*111.);

%% For Haida Gwaii, 
offsetdeg = -0.5;

%% load plate:
pa = load('pacific_plate_near_haida.txt');
lat_up  = [pa(1,2):(pa(end,2)-pa(1,2))/(size(pa,1)*10):pa(end,2)];
pa_up(:,2) = lat_up
pa_up(:,1) = interp1(pa(:,2),pa(:,1),lat_up);

pa_prior = [pa_up,pa_up];
pa_prior(:,1) = pa_prior(:,1)+offsetdeg;
pa_prior(:,3) = pa_prior(:,3)+offsetdeg+width;

figure();hold on;
plot(pa(:,1),pa(:,2),'-b');
plot(pa_prior(:,1),pa_prior(:,2),'-r');
plot(pa_prior(:,3),pa_prior(:,4),'-r');

save('bounds_based_on_pacific_plate.txt','pa_prior','-ascii');
return;
