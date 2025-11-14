function [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
xlm = xymarg(1);
xrm = xymarg(2);
yum = xymarg(3);
ybm = xymarg(4);

spw = (1 - xlm - xrm - (nx-1)*xim)/nx;
sph = (1 - ybm - yum - (ny-1)*yim)/ny;
ii = 1;
for j = 1:ny;
for i = 1:nx;

    loc(1,ii) = xlm + (i-1)*spw + (i-1)*xim;
    loc(2,ii) = 1 - yum - j*sph - (j-1)*yim;
    ii = ii + 1;

end; end;
return;
