function debug_spher;

load field.dat
load fieldb.dat
nfrq = 9;
nbnd = 9;
nline = 4;

for ifrq = 1:nfrq

    a(ifrq,:) = field(nline*(ifrq-1)+1,:) + i*field(nline*(ifrq-1)+2,:);
    kr1(ifrq,:) = field(nline*(ifrq-1)+3,:);
    kr2(ifrq,:) = field(nline*(ifrq-1)+4,:);
    
    b(ifrq,:) = fieldb(nline*(ifrq-1)+1,:) + i*fieldb(nline*(ifrq-1)+2,:);
    kr1b(ifrq,:) = fieldb(nline*(ifrq-1)+3,:);
    kr2b(ifrq,:) = fieldb(nline*(ifrq-1)+4,:);

end

len1 = field(4*nfrq+1,:);
len2 = field(4*nfrq+2,:);
ra = field(4*nfrq+3,:)*1000;
theta = field(4*nfrq+4,:);
theta = fliplr(theta);

len1b = fieldb(4*nfrq+1,:);
len2b = fieldb(4*nfrq+2,:);
rb = fieldb(4*nfrq+3,:);
thetab = fieldb(4*nfrq+4,:);
thetab = fliplr(thetab);

mean(diff(ra))
mean(diff(rb))

rb
size(b)

figure();
hold on;box on;
plot(ra,len1,'b');
plot(ra,len2,'b');
plot(rb,len1b,'k');
plot(rb,len2b,'k');


figure();
hold on;box on;
plot(ra,kr1,'b');
plot(ra,kr2,'b');
plot(rb,kr1b,'k');
plot(rb,kr2b,'k');


figure();
title('Real Part')
plot(ra,real(a(2,:)))
hold on;
plot(rb,real(b(2,:)),'--k')


figure();
title('Imaginary Part')
plot(ra,imag(a(2,:)))
hold on;
plot(rb,imag(b(2,:)),'--k')

x = unwrap(angle(a(2,:)));
y = unwrap(angle(b(2,:)));
figure();
plot(x)
hold on;
plot(y,'--k')

figure();
title('Abs Field')
plot(ra,abs(a(2,:)))
hold on;
plot(rb,abs(b(2,:)),'--k')

for jr=1:length(rb)
    for ifrq=1:nfrq

        Rscb(ifrq,jr) = (b(ifrq,jr) - ...
                         exp(i*kr1b(ifrq,jr))/len1b(jr)) / ...
                        (exp(i*kr2b(ifrq,jr))/len2b(jr));

    end

end

Rsb = abs(fliplr(Rscb));
BLb = -20*log10(Rsb);

Rsbb = Rsb.*Rsb;
Rsbav = sqrt(sum(Rsbb,1)./(nfrq));
BLbav = -20*log10(Rsbav);

%
% Load data averaged in Fortran code
%
load spher_refldBb.dat;
BLd = (spher_refldBb(1:end-1,:));
angd = (spher_refldBb(end,:));
%figure();
%hold on;box on;
%plot(angd,BLd,'-r')

%load pack1.txt
load sim_A_1.txt
BLe = sim_A_1(1:end-1,:);
ange = sim_A_1(end,:);
% plot(angd,pack1,'-k')
% BLe = pack1(1:5,:);
for ifrq=1:5

    ifrq
    figure();
    hold on;box on;
    plot(angd,BLd(ifrq,:),'-r')
    plot(ange,BLe(ifrq,:),'--k')
    legend('jan test','ref')

end

%ange = angd;
%save angs ange

return;
