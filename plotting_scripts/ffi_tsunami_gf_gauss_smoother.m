function [dep,ran,image2] = ffi_tsunami_gf_gauss_smoother(deltr,image);

[NDEP,NRAN]=size(image);

res = 20.;
kmpdeg = 111.;
sec = 1./3600.;
kmps = kmpdeg*sec;
resm = 20*kmps;
L = deltr;

%% First estimate
xgrd = [0:resm:deltr];
nx=length(xgrd);
%% Update estimate
xgrd = [0:deltr/(nx-1):deltr];
nx=length(xgrd);
over = 3*nx;

ltape = over*2; % taper length (should be two times the overlap)
nsrc = ltape + nx;
deltr2 = deltr + 2*over; % length of source patch (overlapping grid, so is larger than L by the amount of taper length)
ndtot = NDEP*nx+ltape;
nrtot = NRAN*nx+ltape;
image2 = zeros(ndtot,nrtot);

dep = [-resm*over:resm:(-resm*over)+(ndtot-1)*resm];
ran = [-resm*over:resm:(-resm*over)+(nrtot-1)*resm];

ss = ones(nsrc,nsrc); % source array

%% Now design a cos taper:
tx=[-floor(nsrc/2):floor(nsrc/2)]*deltr/(nx-1);
C=[L^2,0;0,L^2];

%% Apply taper to source array ss:
for i=1:nsrc;
  for j=1:nsrc;
    xx=[tx(i),tx(j)];xxc=[0,0];
    ss(i,j)=exp(-0.5*(xx-xxc)*inv(C)*(xx-xxc)');
end;end;
%figure,imagesc(tx,tx,ss);

%% Now stack a bunch of sources to make sure they sum to unity.
js = 1;
for j=1:NRAN;
is = 1;
for i=1:NDEP;
  image2(is:is+nsrc-1,js:js+nsrc-1) = image2(is:is+nsrc-1,js:js+nsrc-1) + ss*image(i,j);
  %image2(is:is+nsrc-1,js:js+nsrc-1) = image2(is:is+nsrc-1,js:js+nsrc-1) + ss;
  is = is + nx;
end;
js = js + nx;
end;

%figure();
%imagesc(image2);
%axis('equal');
%colormap(darkb2r(-2,6));colorbar;axis equal;set(gca,'YDir','reverse')
%set(gca,'XLim',[0,500]);
%set(gca,'YLim',[0,950]);

return;
