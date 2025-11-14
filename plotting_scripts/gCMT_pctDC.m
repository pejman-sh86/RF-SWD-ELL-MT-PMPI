%% Some gCMT solutions

%% Bohol 2013:
mt=[4.970 -2.910 -2.060 -1.020 -0.297 -3.530];
MT=[mt(1),mt(4),mt(5);mt(4),mt(2),mt(6);mt(5),mt(6),mt(3)]
e=abs(eig(MT))
pDC=100*(1-2*(min(e)/max(e)))

%% Amberley 2016
mt=[3.520 1.630 -5.150 -1.380 4.390 -1.930];
MT=[mt(1),mt(4),mt(5);mt(4),mt(2),mt(6);mt(5),mt(6),mt(3)]
e=abs(eig(MT))
pDC=100*(1-2*(min(e)/max(e)))
