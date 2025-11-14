%----------------------------------------------------------
% Plot nested sampling results
%----------------------------------------------------------
function [] = plot_nest();

i_save = 0;
ippd2 = 1;

Tstar2 = 1;
i_scaley = 0;
itrue = 0;
nx = 2;
ny = 2;

% INPUT:
%ppd     = 'sample.mat';
%ppd2    = '../../ais/ais_mh/sample.mat';
ppd     = 's03_sample.mat';
ppd2    = '../../ais/ais_tomo/sample.mat';

% OUTPUT:
if(i_save == 1)
   mapfile  = 'nestmap.txt';
   meanfile = 'mean.txt';
   nestfile = 'nest.eps';
   plotfile = 'marg.eps';
end

nsubfig = nx*ny;
xim = 0.03;
yim = 0.1/ny;
xymarg = [0.04 0.04 0.04 0.14];
nbin1 = 100;
nbin2 = 100;

opts = struct('bounds','tight','linestylemap','bw','LockAxes',1, ...
              'Width',12,'Height',1*ny,'Color','bw',...
              'Renderer','painters',...
              'FontMode','fixed','FontSize',12,'FontEncoding','adobe');

%% Toy lighthouse
F(1).xticks = [0 0.3 0.6 1.0];
F(2).xticks = [0 0.3 0.6 1.0];
F(3).xticks = [0 0.3 0.6 1.0];
F(4).xticks = [0 0.3 0.6 1.0];

%% Toy tomography
%F(1).xticks = [0 0.3 0.6 1.0];
%F(2).xticks = [0 0.3 0.6 1.0];
%F(3).xticks = [0 0.3 0.6 1.0];
%F(4).xticks = [0 0.3 0.6 1.0];

load(ppd);
m        = A(:,6:end);
logL     = A(:,1);
logwidth = A(:,2);
logwt    = A(:,3);
H        = A(:,4);
logZ     = A(:,5);
w = exp(logwt-logZ(end));

if(ippd2 == 1)
   load(ppd2);
   m2 = A(:,3:end);
   logwt2 = A(:,2);
   logZ2  = A(1,5);
   w2 = exp(logwt2-logZ2);
end

%
% Plot nested sampling convergence plot
%
gc1=figure(1);
plot(exp(w));
xlabel('No. iteration');
ylabel('Importance weight');
set(gca,'XLim',[1 length(w)])

[logL_max i_max] = max(logL);
xmap = m(i_max,:);
if(i_save == 1)
    save(mapfile,'xmap','-ASCII');
end

npar = size(m,2);
nlay = (npar-3)/4;

%% Toy lighthouse
plo = [-2.0, 0.0];
phi = [ 2.0, 2.0];

%% Toy tomography
%plo = [0.0, 0.0, 0.0, 0.0,...
%       0.0, 0.0, 0.0, 0.0,...
%       0.0, 0.0, 0.0, 0.0];
%phi = [1.0, 1.0, 1.0, 1.0,...
%       1.0, 1.0, 1.0, 1.0,...
%       1.0, 1.0, 1.0, 1.0];

pplo = plo;
pphi = phi;

nmods = size(m,1);
nfig = ceil(npar/nsubfig);

%----------------------------------------------------------
%  Calc histograms:
%----------------------------------------------------------

for ipar = 1:npar

   [n1(:,ipar),lim(:,ipar)] = hist_wt(m(:,ipar),w,nbin1);
   if(ippd2 ==1)
       [n2(:,ipar),lim2(:,ipar)] = hist_wt(m2(:,ipar),w2,nbin2);
%      [n2(:,ipar),lim2(:,ipar)] = hist_tstar(m2(:,ipar),E2,nbin2,Tstar2);
   end

end

%----------------------------------------------------------
%  Plot histograms:
%----------------------------------------------------------

[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

gc2=figure(2);
isubfig = 1;
ipar = 1;
for iy = 1:ny
for ix = 1:nx

   if(ipar <= npar)
      if(ipar == npar-2) 
         isubfig = isubfig + 1;
      end;
      if(ipar > npar-3) 
         ylims_tmp(iy,ix+1) = max(n1(:,ipar));
      else
         ylims_tmp(iy,ix) = max(n1(:,ipar));
      end
      subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
      hold on;box on;

      lim_2 = [min(m(:,ipar));min(m(:,ipar));
               lim(:,ipar); max(m(:,ipar)); max(m(:,ipar))];
      n1_2 = [0; n1(1,ipar); n1(:,ipar); n1(end,ipar); 0];
      fill(lim_2,n1_2,[0.8,0.8,0.8]);

      if(ippd2 ==1)
      lim2_2 = [min(m2(:,ipar));min(m2(:,ipar));
               lim2(:,ipar); max(m2(:,ipar)); max(m2(:,ipar))];
      n2_2 = [0; n2(1,ipar); n2(:,ipar); n2(end,ipar); 0];
      stairs(lim2_2,n2_2);
      end

      set(gca,'layer','top')
      ipar = ipar +1;
      isubfig = isubfig + 1;
   end
end
end

%ylims = max(max(ylims_tmp));
ylims = max(ylims_tmp);

ipar = 1;
isubfig = 1;
for iy = 1:ny
for ix = 1:nx
   if(ipar <= npar)
      if(ipar == npar-2) 
         isubfig = isubfig + 1;
      end;
      subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
      set(gca,'YTick',[],'FontSize',14);
      set(gca,'XLim',[pplo(ipar) pphi(ipar)],'FontSize',14);
%      plot([xmap(ipar) xmap(ipar)],[0 ylims(2)],'--k')

      if(i_scaley == 1)
         if(ipar > npar-3) 
            set(gca,'YLim',[0 ylims(ix+1)+.05*ylims(ix+1)]);
            if(itrue == 1)
               plot([mtrue(ipar) mtrue(ipar)],[0 ylims(ix+1)+.05*ylims(ix+1)],'-k')
            end
         else
            set(gca,'YLim',[0 ylims(ix)+.05*ylims(ix)]);
            if(itrue == 1)
               plot([mtrue(ipar) mtrue(ipar)],[0 ylims(ix)+.05*ylims(ix)],'-k')
            end
         end 
      else
        if(iy > (npar-3)/4) 
            set(gca,'YLim',[0 ylims_tmp(iy,ix+1)+.05*ylims_tmp(iy,ix+1)]);
            if(itrue == 1)
               plot([mtrue(ipar) mtrue(ipar)],[0 ylims_tmp(iy,ix+1)+.05*ylims_tmp(iy,ix+1)],'-k')
            end
         else
            set(gca,'YLim',[0 ylims_tmp(iy,ix)+.05*ylims_tmp(iy,ix)]);
            if(itrue == 1)
               plot([mtrue(ipar) mtrue(ipar)],[0 ylims_tmp(iy,ix)+.05*ylims_tmp(iy,ix)],'-k')
            end
         end 
      end
%      set(gca,'YLim',[0 30]);
      ipar = ipar +1;
      isubfig = isubfig + 1;
   end
end
end

if(i_save == 1)
   exportfig(gc1,nestfile,opts);
   exportfig(gc2,plotfile,opts);
end

return;
