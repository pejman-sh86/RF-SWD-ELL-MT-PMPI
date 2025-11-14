

N=100;
L = 15;
x=[0:1:100];
y=[0:1:100];
X=[0:15:101];
Y=[0:15:101];

xi=50;
yi=50;
C=[L^2,0;0,L^2];

for i=1:101;
  for j=1:101;
    xx=[x(i),y(j)];xxc=[xi,yi];
    G(i,j)=exp(-0.5*(xx-xxc)*inv(C)*(xx-xxc)');
end;end;
figure();
imagesc(G);colorbar;

GG=zeros(2000,2000);
for i=1:15:1500;
  for j=1:15:1500;
    GG(i:i+N,j:j+N) = GG(i:i+N,j:j+N)+G;
end;end;
figure();
imagesc(GG)
