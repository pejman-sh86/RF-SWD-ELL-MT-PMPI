%function test_rayfast;

paths = [1 0 0 0 0 0];
z_in  = [0 10 20 30 140];
c_in  = [1510 1500 1505 1520 1525];
zs_in = 0.5;
zr_in = 122.;
r_tol = 0.01;
r = 7000.;

rayfast(paths,z_in,c_in,zs_in,zr_in,r,r_tol)


%return;
