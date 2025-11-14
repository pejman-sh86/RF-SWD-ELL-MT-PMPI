%%
%% This of course will have to be properly adapted to lat lon etc. to make a source .grd file which Jakir uses.
%%


N  = 30; % Number of source patches in lat and lon
l  = 0.5; % fine grid spacing (km)
L  = 30; % coarse grid spacing between source patches (grid points)
over = 5; % amount of overlap (grid points)

ltape = over*2; % taper length (should be two times the overlap)
L2 = L + 2*over; % length of source patch (overlapping grid, so is larger than L by the amount of taper length)

X = zeros(N*(L+1)); % total source region

ss = ones(L2); % source array

%% Now design a cos taper:
tx=[-pi/2:pi/(ltape-1):pi/2];
ty=[0:pi/(ltape-1):pi];
ta = (sin(tx)+1)/2;
tb = (cos(ty)+1)/2;
tape = ones(L2,1);
tape(1:ltape) = ta;
tape(end-ltape+1:end) = tb;

%% Apply taper to source array ss:
for i=1:L2;
  ss(:,i) = ss(:,i).*tape;
  ss(i,:) = ss(i,:).*tape';
end;
figure,surf(ss);

%% Now stack a bunch of sources to make sure they sum to unity.
is = 1;
for i=1:N;
js = 1;
for j=1:N;
  X(is:is+L2-1,js:js+L2-1) = X(is:is+L2-1,js:js+L2-1) + ss;
  js = js + L;
end;
is = is + L;
end;

%figure,surf(X);
figure,imagesc(X);
grid off;

