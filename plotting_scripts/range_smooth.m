function [] = range_smooth();

idx = 5

tmp = load('spher_noave.dat');
tmp2 = load('sph_gauss_av.txt');
tmp3 = load('spher_refldBb.dat');

Rp = tmp2(idx,:);
ang_Rp = tmp2(end,:);
Rint = tmp3(idx,:);
ang_Rint = tmp3(end,:);
N = fliplr(tmp(idx,:));
nran = length(N);
nfrq = size(N,1);
r = tmp(end-1,:);
th = tmp(end,:);

alpha = 0.022;
for j=1:2
disp([j alpha])
dr = mean(diff(r));
for i = 1:nran

    r0 = r(i);
    NR(i) = sum(N .* exp(-(r-r0).^2/(alpha*r0)^2)*dr)/...
            sum(exp(-(r-r0).^2/(alpha*r0)^2)*dr);

end;

figure(j);hold on;box on;
plot(ang_Rp,Rp,'-b')
plot(th,NR,'--k')
plot(ang_Rint,Rint,'--r')
legend('processr','matlab','jantest')
alpha = alpha + 0.003;
end

return;
