function compare_profile;

set(0, 'DefaultFigurePaperPosition', [0 0 7 4]); 
showcore = 0;
site20 = 0;
errbar = 0;
barstep = 10;

%
% Use ML model
%
%x(:,1) = [1.6665  1.3563  1.5461  1.0035  1472.7412  1465.8553  0.3799]';
%x(:,1) = [1 1.3142633      1.4938605     0.14312698      1484.1188      1488.0772     0.34090275]';
x(:,1) = [1. 1.4550371      1.6010057      0.8      1483.5776      1453.2873     0.26679658]';
x(:,2) = [1. 1.4550371      1.6010057      1.0      1483.5776      1453.2873     0.26679658]';
x(:,3) = [1. 1.4550371      1.6010057      1.2      1483.5776      1453.2873     0.26679658]';
x(:,4) = [1. 1.4550371      1.6010057      1.4      1483.5776      1453.2873     0.26679658]';
x(:,5) = [1. 1.4550371      1.6010057      1.6      1483.5776      1453.2873     0.26679658]';
%x(:,4) = [1 1.1      1.4976439    0.057720791      1484.4985      1488.6057     0.35936887]';
%x(:,2) = [1.8 1.3461897      1.4986794      0.7795354      1473.4571      1465.3035     0.31035113]';
%x(:,3) = [1.8 1.3585412      1.5054058     0.95254483      1473.4314      1467.3606     0.31351455]';
%x(:,4) = [1.8 1.355651      1.5332277     0.93675168      1473.7369      1464.9685     0.34219434]';
%x(:,5) = [1.8 1.3561251      1.4621229     0.76159934      1472.6783      1472.0898     0.26903145]';
%x(:,6) = [1.8 1.3626426      1.4962751     0.99955341      1474.7739       1472.986     0.28968489]';
%x(:,7) = [1.8 1.1000028      1.4482859     0.15235421      1477.4642      1479.0265     0.41189471]';
%x(:,8) = [1.8 1.1039379      1.4541665     0.15258941      1475.3046      1478.5997     0.42047152]';


nmod = size(x,2)
col = {'b' 'r' 'k' :'k' '-b' '--b' ':b' 'o-b'};
%col = {r --r b --b g --g y --y k}
col = char(col);
%save bla col

sc_5 = 1.003551;
sc_6 = 1.010277;
sr_5 = 1.006025;
sr_6 = 0.996599;

if(showcore == 1)
  C = load('core5.dat');
  D = load('core6.dat');
  C(:,2) = C(:,2)/sc_5;
  C(:,3) = C(:,3)/sr_5;
  D(:,2) = D(:,2)/sc_6;
  D(:,3) = D(:,3)/sr_6;
  if(errbar == 1)
    ErrC_r = 2/100 * C(:,3);
    ErrC_c = 6 * ones(length(C),1);
    ErrD_r = 2/100 * D(:,3);
    ErrD_c = 6 * ones(length(D),1);
  end
  figure(1);
  hold on;box on;
  subplot('Position',[.12 .15 .4 .8]);
  hold on;box on;
  plot(C(:,3),C(:,1),'+g');
  plot(D(:,3),D(:,1),'+r');
  if(errbar == 1)
    errorbarxy(C(1:barstep:end,3),C(1:barstep:end,1),ErrC_r(1:barstep:end),...
    zeros(size(ErrC_r(1:barstep:end))),'b','g');
    errorbarxy(D(3:barstep:end,3),D(3:barstep:end,1),ErrD_r(3:barstep:end),...
    zeros(size(ErrD_r(3:barstep:end))),'b','r');
  end
  subplot('Position',[.56 .15 .4 .8]);
  hold on;box on;
  plot(C(:,2),C(:,1),'+g');
  plot(D(:,2),D(:,1),'+r');
  if(errbar == 1)
    errorbarxy(C(1:barstep:end,2),C(1:barstep:end,1),...
    ErrC_c(1:barstep:end),...
    zeros(size(ErrC_c(1:barstep:end))),'b','g');
    errorbarxy(D(3:barstep:end,2),D(3:barstep:end,1),ErrD_c(3:barstep:end),...
    zeros(size(ErrD_c(3:barstep:end))),'b','r');
  end

end


znorm = zeros(10,1);
lay_thick = 0.01;
znorm=lay_thick/2:lay_thick:1-lay_thick/2;

for imod=1:nmod
%
% model:
%
  
  h_ref = x(1,imod);
  rhot_ref = x(2,imod);
  rhob_ref = x(3,imod);
  nu_ref = x(4,imod);
  ct_ref = x(5,imod);
  cb_ref = x(6,imod);
  alpha_ref = x(7,imod);

  z_ref = [0 znorm] *h_ref;
  rho_ref = rhot_ref + (sin([0 znorm] .* pi/2)).^nu_ref...
                  .* (rhob_ref - rhot_ref);
  c_ref = ct_ref + (cb_ref - ct_ref) .* [0 znorm];

  sub1 = subplot('Position',[.12 .15 .4 .8]);hold on;box on;
  plot(rho_ref,z_ref,col(imod,:),'MarkerSize',8);

  sub2 = subplot('Position',[.56 .15 .4 .8]);hold on;box on;
  plot(c_ref,z_ref,col(imod,:),'MarkerSize',8);

end %imod loop

if site20 == 1
 
  subplot(sub1);
  load core4b-dns.mat;
  
  plot(c(:,2),c(:,1),'-k');
  
end

subplot(sub1);
legend('04','19','20','20 with cov');
set(gca,'XLim',[1.2 1.7]);
set(gca,'XTickLabel',[1.2 1.4 1.6],'XTick',[1.2 1.4 1.6]);
set(gca,'YDir','reverse','YLim',[0 1.5],'FontSize',12);

subplot(sub2);
legend('04','19','20','20 with cov');
set(gca,'XLim',[1460 1490]);
set(gca,'XTickLabel',[1460 1470 1480],'XTick',[1460 1470 1480]);
set(gca,'YDir','reverse','YLim',[0 1.5],'YTick',[],'FontSize',12);

return;
