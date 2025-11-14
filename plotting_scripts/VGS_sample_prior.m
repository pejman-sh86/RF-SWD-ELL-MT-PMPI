
f = [1000];
pmin = [0.22 6.9 -1.4 -4.5]';
pmax = [0.95 8.6 -0.4 -1.3]';
maxpert = pmax - pmin;


for i=1,100000;
    m = pmin
    [cp(i),alfp(i),cs(i),alfs(i),rho(i)] = ...   
    VGSlambda_Mod(mmisc(1),mmisc(2),mmisc(3),mmisc(4),...
          m(ipar),10.^(m(ipar+1)),10.^(m(ipar+2)),10.^(m(ipar+3)),fr);
end;

