function []=compare_var_par_models(filename,n1,n2);
NPL = 6;
zmx = 8.;
load(filename);
B=A;
filename2 = strrep(filename,'voro_sample.mat','sample.mat');
load(filename2);
m1 = A(n1,5:A(n1,4)*NPL+NPL-1+4);
m2 = A(n2,5:A(n2,4)*NPL+NPL-1+4);
voro1 = B(n1,5:B(n1,4)*NPL+4);
voro2 = B(n2,5:B(n2,4)*NPL+4);

minlim = [0.300, log(5.e7), log(0.001),log(1.5e-3),log(1.e6)];
maxlim = [0.910, log(5.e9), log(0.050),log(0.37),  log(8.e8)];

disp([A(n1,4),B(n1,4),A(n2,4),B(n2,4)]);

mtru =[0.1450,  0.8, 18.0, -6.6036, -4.1730, 14.0,...
       0.6000,  0.7, 19.0, -6.6036, -4.1730, 14.0,...
       0.9450,  0.6, 20.0, -6.6036, -4.1730, 14.0,...
       1.4000,  0.6, 20.0, -6.6036, -5.2100, 18.6,...
       4.4450,  0.5, 21.0, -6.6036, -5.2100, 18.6,...
                0.4, 22.0, -6.6036, -5.2100, 20.0];

fig1=figure();hold on;box on;
figw = 12;
figh = 6;
set(fig1,'PaperUnits','inches','PaperPosition',[0 0 figw figh]);
nx = NPL-1;
ny = 1;
xim = 0.01;
yim = 0.00;
xymarg = [0.1 0.04 0.04 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
for ipar = 1:5;
  sub1 = subplot('Position',[loc(1,ipar) loc(2,ipar) spw sph]);hold on; box on;
  h1=plprof(mtru,zmx,NPL,'k',ipar,0);
  set(h1,'LineWidth',4,'Color',[.6,.6,.6]);
  h2=plprof(m1,zmx,NPL,'b',ipar,0);
  set(h2,'LineWidth',2);
%  h3=plprof(m2,zmx,NPL,'--r',ipar,0);
%  set(h3,'LineWidth',2);
  set(gca,'YDir','reverse');
  for ivo = 1:B(n1,4);
    if(voro1((ivo-1)*NPL+1+ipar) > -99.);
    plot(voro1((ivo-1)*NPL+1+ipar),voro1((ivo-1)*NPL+1),'ob',...
         'LineWidth',2,'MarkerSize',10);end;
  end;
%  for ivo = 1:B(n2,4);
%    if(voro2((ivo-1)*NPL+1+ipar) > -99.);
%    plot(voro2((ivo-1)*NPL+1+ipar),voro2((ivo-1)*NPL+1),'or',...
%         'LineWidth',2,'MarkerSize',10);end;
%  end;
  set(gca,'XLim',[minlim(ipar) maxlim(ipar)]);
  if(ipar > 1);set(gca,'YTickLabel',[]);end;
end;

print(fig1,'coplexity.png','-dpng');
return;
