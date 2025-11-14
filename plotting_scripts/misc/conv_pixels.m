function [] = conv_pixels();

corefile = 'core4b-dns.raw';
savefile = 'core4b-dns.mat';
%corefile = 'co4bquin-dns.raw';
%savefile = 'co4bquin-dns.mat';
A = load(corefile);

px1 = A(1,1);
px2 = A(2,1);
pz1 = A(3,1);
pz2 = A(4,1);
x1  = A(1,2);
x2  = A(2,2);
z1  = A(3,2);
z2  = A(4,2);

dx = (x2-x1)/(px2-px1)
dz = (z2-z1)/(pz2-pz1);

B = A(5:end,:);

for iz = 1:length(B)

  z(iz) = (B(iz,1)-B(1,1)) * dz;
  x(iz) = x1 + (B(iz,2)-px1) * dx;
  (B(iz,2)-px1)*x1

end
c = [z',x']
save(savefile,'c')

return;

