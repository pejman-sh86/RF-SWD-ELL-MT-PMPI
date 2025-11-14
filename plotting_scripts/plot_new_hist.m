function [] = plot_new_hist();

npar = 16;
ipar = 1;

nx = 4;
ny = 4;
nsubfig = 16;
nfig = ceil(npar/nsubfig);

xim = 0.06;
yim = 0.00;
[loc,spw,sph] = get_loc(nx,ny,xim,yim);

for ifig = 1:nfig
    figure(ifig);
    isubfig = 1;

    while (isubfig <= nsubfig & ipar <= npar)

        if(ipar == npar-2)
            isubfig = isubfig + 1;
        end;

        subplot('Position',[loc(1,isubfig) loc(2,isubfig) spw sph]);
        box on; hold on;

        ipar = ipar +1;
        isubfig = isubfig + 1;

    end;
end;
return;
