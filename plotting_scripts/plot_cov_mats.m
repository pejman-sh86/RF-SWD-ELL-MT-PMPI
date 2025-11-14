
itrue_covs = 0;

fig=figure;
cmap = colormap( parula(1000) );
% cmap = colormap( bone(1000) );
% cmap = flip(cmap,1);
% % cmap(1:1,:) = ones(1,3);%[1 1 1];
colormap(cmap);
filebase = 'HON';
Npanel = 4;   
font_size = 30;
% covlabel = {'RF', 'SWD', 'H/V', 'MT'};
covlabel = {'MT', 'SWD', 'RF'};
% covlabel = {'SWD', 'RF','MT'};
   nx = Npanel;
   ny = 2;
   ipanely=2;
   xim = 0.01;
   yim = 0.08;
   xymarg = [0.1 0.04 0.04 0.1];
   [loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);
   loc(2,:) = loc(2,:) + 0.1;
%    spw = spw - 0.005;
   spw2 = spw;
   spw = spw * 7./8.;
   spw = spw - 0.005;
% for ipanel=2:nx
%         loc(1,ipanel+(ipanely-1)*nx) = loc(1,ipanel+(ipanely-1)*nx) - (ipanel-1)*(spw2-spw);   
% end
%    XTick = { [7.e-5, 1.e-4, 1.5e-4], [7.e-3, 1.e-2, 1.5e-2], [7.e-3, 1.e-2, 1.5e-2] };
%    exponent = [-4, -2, -2]; 
% XTick = { [1., 2., 3.], [1., 1.5], [1., 1.1] }; 
% newlims = [-3.7795e-09   2.0568e-08; ...
%            -4.4514e-09   4.4514e-09; ...
%            -1.9861e-06   8.7434e-06; ...
%            -5.9396e-03   1.4458e-02];  %real
% newlims = [-5.9396e-03   1.4458e-02; ...
%            -4.4514e-09   4.4514e-09; ...
%            -1.9861e-06   8.7434e-06; ...
%            -5.9396e-03   1.4458e-02];

%%
parfile  = strcat(filebase,'_parameter.dat');
[IMAP,ICOV,iar,i_varpar,irv,itr,iswd, iell, imt,izmt, ivref,ivpvs,ISMPPRIOR,...
          IEXCHANGE,idip,NDAT_SWD,NMODE,NDAT_ELL, NMODE_ELL, NDAT_MT,NTIME,NSRC,NVMN,NVMX,ICHAINTHIN,...
          NKEEP,NPTCHAINS1,dTlog,hmx,hmin,armxH,armxV,TCHCKPT,shift,...
          sampling_dt,dVs,dVpVs, sigmamin, sigmamax, sdmn, sdmx, ...
          ntr,baz]=rf_read_parfile(parfile);
NRF = ntr;
%%

NCOV = 0;
if (irv==-1); NCOV=NCOV+1; end
if (iswd==1); NCOV=NCOV+1; end
if (iell==1); NCOV=NCOV+1; end
if (imt==1); NCOV=NCOV+2; end
%%
covmats = cell(1, NCOV);
covmatlabels = cell(1, NCOV);
ipos = 1;
if (imt==1)
    if itrue_covs
        CdMT2 = load(strcat(filebase,'_CdZMT_true.dat'));
    else
        CdMT2 = load(strcat(filebase,'_CdZMT.dat'));
    end
    CdMT = complex(  CdMT2(:,1:NDAT_MT), CdMT2(:,NDAT_MT+1:2*NDAT_MT) );
    CdMTR = CdMT2(:,1:NDAT_MT);
%     covmats{ipos} = CdMTR;
    covmats{ipos} = flip( flip(CdMTR), 2);
    covmatlabels{ipos} = 'REAL(C_{MT})';
    ipos = ipos + 1;
    CdMTI = CdMT2(:,NDAT_MT+1:2*NDAT_MT);
%     covmats{ipos} = CdMTI;
    covmats{ipos} = flip( flip(CdMTI), 2);
    covmatlabels{ipos} = 'IMAG(C_{MT})';
    ipos = ipos + 1;
end
if (iswd==1) 
    if itrue_covs
        CdSWD = load(strcat(filebase,'_CdSWD_true.dat')); 
    else
        CdSWD = load(strcat(filebase,'_CdSWD.dat'));
    end
    covmats{ipos} = CdSWD;
    covmatlabels{ipos} = 'C_{SWD}';
    ipos = ipos + 1;
end
if (irv==-1)
    if itrue_covs
        CdRF = load(strcat(filebase,'_Cd_true.dat'));
    else
        CdRF = load(strcat(filebase,'_Cd.dat'));
    end
    covmats{ipos} = CdRF;
    covmatlabels{ipos} = 'C_{RF}';
end

%%
lims = zeros(NCOV,2);
ilabel = 1;

for i=1:nx
     h1 = subplot('Position',[loc(1,i+(ipanely-1)*nx) loc(2,i+(ipanely-1)*nx) spw sph]);
     hold on; box on;
     set(gca,'FontSize',font_size);
     c = colorbar('northoutside');
     caxis(h1,newlims(i,:))
     c.Limits = newlims(i,:);
%      if(i > 1);set(gca,'YTickLabel',[]);end;
     if(i > NCOV);set(gca,'YTickLabel',[]);end
     if(i>NCOV); set(gca,'XTickLabel',[]);end

     if i<=NCOV
         imagesc(covmats{i})
         set(gca,'XLim',[1 size(covmats{i},1)],'LineWidth',1);
         set(gca,'YLim',[1 size(covmats{i},1)],'YDir','reverse');
         lims(i,:) = c.Limits;
         
         if(i == 1);ylabel('Data points');end
         xlabel({'Data points',covmatlabels{i}});
%          label = strcat( '{',covlabel(ilabel) );
%          label = strcat( label, '}' );
%          xlabel( strcat('s_', label) )
%          ilabel = ilabel+1;


     end
end

