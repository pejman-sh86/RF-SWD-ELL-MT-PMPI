% File name:  hist95com.m
% plot histogram of pdfs
% prepared for display the result from FGS
%
% Yongmin Jiang, 2006, August. 08 

% Paper type of Figure
set(0,'DefaultFigurePaperType','a4');
set(0, 'DefaultFigurePaperPosition', [0 0 8.5 7]);

% input .m file's (root) name 
FileRoot=input('FileRoot ? ','s');

% execute .m file to get number of model parameters/true value/search bound 
eval(FileRoot);

%number of model parameters+1
npar=ndim+1;

%number of bins
nrbin = 50;

%number of dataset
nset = 3;

% Read from file
for iidata = 1:nset
    for ifile=1:2
        FileName=[FileRoot '_hist.mat' int2str(iidata)];  % result from sampling 1
        load(FileName);
        A1=A;
    end
    switch iidata
        case 1
            A = [A1];
            npx1=size(A,2);      % total sample numbers for each parameter
            [Emin,p]=min(A(1,:));  %Find the element of minimum energy
            xmap1=A(2:npar,p);     %Find the MAP model according to p 
        case 2
            B = [A1];
            npx2=size(B,2);      % total sample numbers for each parameter
            [Emin,p]=min(B(1,:));  %Find the element of minimum energy
            xmap2=B(2:npar,p);     %Find the MAP model according to p 
        case 3
            C = [A1];
            npx3=size(C,2);      % total sample numbers for each parameter
            [Emin,p]=min(C(1,:));  %Find the element of minimum energy
            xmap3=C(2:npar,p);     %Find the MAP model according to p 
    end
end;   

% kick out the energy
dat1=A(2:npar,:);
dat1=dat1';

% name of parametrs
%       1           2             3             4             5            6
%       7           8             9             10            11           12
%       13          14            15            16            17           18 
%       19          20            21            22            23           24 

Tp0=['cp1      '; 'cp2      '; 'cp2(g)   '; 'cs1      '; 'cs2      '; 'cs2(g)   '; ...
     '\alphap1 '; '\alphap2 '; '\alphas1 '; '\alphas2 '; '\rho1    '; '\rho2    '; ...
     '\rho2(g) '; 'WD       '; 'H        '; 'Range    '; 'SD       '; 'RD       '; ...
     'tilt     '; 'cpb      '; 'csb      '; '\alphapb '; '\alphasb '; '\rhob    '];
% unit of parameters
Tu0=['(m/s)     '; '(m/s)     '; '(m/s/m)   '; '(m/s)     '; '(m/s)     '; '(m/s/m)   '; ...
     '(dB/m/kHz)'; '(dB/m/kHz)'; '(dB/m/kHz)'; '(dB/m/kHz)'; '(g/cm^3)  '; '(g/cm^3)  '; ...
     '(g/cm^3/m)'; '(m)       '; '(m)       '; '(km)      '; '(m)       '; 'RD(m)     '; ...
     '(\circ)   '; '(m/s)     '; '(m/s)     '; '(dB/m/kHz)'; '(dB/m/kHz)'; '(g/cm^3)  '];

% number of subplots in a row
if ndim<=3
    div=1;
elseif (ndim==4) 
    div=2;
elseif ((ndim<=6) & (ndim>=5))
    div=3;
elseif ((ndim<=8) & (ndim>=7))
    div=4;
elseif (ndim==9) 
    div=3;
elseif (ndim==10)
    div=5;  
elseif ((ndim<=12) & (ndim>10))
    div=5;
else
    div=3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              Plot                               % 
% history plot
% plot the individual sampling history
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)

FontWeight='normal';
FontSize=10;

% number of subplots in a column
ysub=ceil(ndim/div);
clf;
ipl=0;
jplot=0;    % subplot index

xdiff = abs(maxlim - minlim);
delta = xdiff/nrbin;
N = length(A);
for i =0:nrbin
  edges(i+1,:) = minlim + xdiff/nrbin * i;
end

ipl = 0;
for iidx=[14 15 1 20 7 22 11 24 17 16]    % 1 layer, homo a,rho
    ipl=ipl+1;
    jplot=jplot+1;
    subplot(ysub,div,jplot);

    % plot the outline of the marginal distribution
    %--Smooth Style
    % mn--vector contains the frequency counts;  
    % xout-- vector contains the bin locations 
    [mn1,xout1]=hist(A(ipl+1,:),edges(:,ipl)); 
    if ipl==1
        ymax = 2.5*max(mn1);
    end
    mn1 = mn1/ymax;
    h = area(xout1,mn1 + 2,2);
    set(h,'FaceColor',[0.6 0.6 0.6]);
    set(h,'LineWidth',1.0);
    set(gca,'XLim',[minlim(ipl) maxlim(ipl)],'YLim',[0 3],...
        'YTick',[0 1 2],'YTickLabel',{'c';'b';'a'},'FontSize',11);
    hold on;
    plot([xmap1(ipl) xmap1(ipl)],[2 3],'-.r','LineWidth',1.5);
    exp_val = mean(A(ipl+1,:));
    plot([exp_val exp_val],[2 3],'--b','LineWidth',1.5);

    [mn2,xout2]=hist(B(ipl+1,:),edges(:,ipl));  
    if ipl==1
        ymax = 3.0*max(mn2);
    end
    mn2 = mn2/ymax;
    h = area(xout2,mn2 + 1,1);
    set(h,'FaceColor',[0.6 0.6 0.6]);
    set(h,'LineWidth',1.0);
    plot([xmap2(ipl) xmap2(ipl)],[1 2],'-.r','LineWidth',1.5);
    exp_val = mean(B(ipl+1,:));
    plot([exp_val exp_val],[1 2],'--b','LineWidth',1.5);

    [mn3,xout3]=hist(C(ipl+1,:),edges(:,ipl));  
    if ipl==1
        ymax = 3.0*max(mn3);
    end
    mn3 = mn3/ymax;
    h = area(xout3,mn3);
    set(h,'FaceColor',[0.6 0.6 0.6]);
    set(h,'LineWidth',1.0);
    plot([xmap3(ipl) xmap3(ipl)],[0 1],'-.r','LineWidth',1.5);
    exp_val = mean(C(ipl+1,:));
    plot([exp_val exp_val],[0 1],'--b','LineWidth',1.5);

    Text1=[char(cellstr(Tp0(iidx,:))) char(cellstr(Tu0(iidx,:)))];
    xlabel(Text1,'FontSize',11);
    hold off;
end;

% Text  
axes('position',[0,0,1,1]); axis('off')
FileDesc=strrep(FileRoot,'_',' ');
text(.05,.02, ['File: ',FileDesc,'; Date: ',date],'fontsize',8)

