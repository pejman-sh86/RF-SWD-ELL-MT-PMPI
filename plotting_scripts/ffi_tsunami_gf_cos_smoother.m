function [dep,ran,image2] = ffi_tsunami_gf_cos_smoother(deltr,image);


[NDEP,NRAN]=size(image);

res = 60.;
kmpdeg = 111.;
sec = 1./3600.;
kmps = kmpdeg*sec;
resm = res*kmps;

%% First estimate
xgrd = [0:resm:deltr];
nx=length(xgrd);
%% Update estimate
xgrd = [0:deltr/(nx-1):deltr];
nx=length(xgrd);
over = floor(nx/2);

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
tx=[-pi/2:pi/(ltape-1):pi/2];
ty=[0:pi/(ltape-1):pi];
ta = (sin(tx)+1)/2;
tb = (cos(ty)+1)/2;
tape = ones(nsrc,1);
tape(1:ltape) = ta;
tape(end-ltape+1:end) = tb;

%% Apply taper to source array ss:
for i=1:nsrc;
  ss(:,i) = ss(:,i).*tape;
  ss(i,:) = ss(i,:).*tape';
end;
%figure,surf(ss);

%% Up-sample image
%for j=1:NRAN;
%    for i=1:NDEP;
%        image2((i-1)*nx+1:i*nx,(j-1)*nx+1:j*nx) = image(i,j);
%    end;
%end;
%figure();imagesc(image2);

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
%imagesc(image3);
%axis('equal');

return;
