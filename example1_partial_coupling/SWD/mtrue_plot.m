
z = -[-200.0 -170. -35. -8. -1. 0.];
Vs = [4.7 4.7 4.7 3.8 3.2 2.];
VpVs = [2. 2. 1.7 1.7 1.7 1.8];
logsigma = [-2. -2. -3. -3. -3. -1.];

fig = figure(1000);
hold on; box on;

left = 0.2;
bottom= 0.14;
w = 0.1;
h = 0.8;
spw = 0.2;

loc1 = [left; bottom];
loc2 = [left+w+spw; bottom];
loc3 = [left+2*w+2*spw; bottom];

m1 = subplot('Position', [loc1(1,1), loc1(2,1), w, h]);
hold on; box on;
subplot(m1)
hold on; box on;
set(gca, 'Fontsize', 21)
stairs(-logsigma, z);
xlabel('log_{10}\rho(\Omega m)', 'Fontsize', 21)
ylabel('Depth (km)', 'Fontsize', 21)
xlim([0.0 5.])
%ylim([-400 0])
set(gca, 'YDir', 'reverse')

m2 = subplot('Position', [loc2(1,1), loc2(2,1), w, h]);
hold on; box on;
subplot(m2)
hold on; box on;
set(gca, 'Fontsize', 21)
stairs(Vs, z);
xlabel('V_s(Km/s)', 'Fontsize', 21)
ylabel('Depth (km)', 'Fontsize', 21)
xlim([1.5 5.])
%ylim([-400 0])
set(gca, 'YDir', 'reverse')

m3 = subplot('Position', [loc3(1,1), loc3(2,1), w, h]);
hold on; box on;
subplot(m3)
hold on; box on;
set(gca, 'Fontsize', 21)
stairs(VpVs, z);
xlabel('V_p/V_s ratio', 'Fontsize', 21)
ylabel('Depth (km)', 'Fontsize', 21)
xlim([1.65 2.1])
%ylim([-400 0])
set(gca, 'YDir', 'reverse')