% ---------------------------------------------------------------
%             BL processing for Holland data
% ---------------------------------------------------------------
%
%    Record:
%       Date            by            last  change
%     ========    ==============      ============
%      Apr 05      jand@uvic.ca        Aug 14 07
%
% ---------------------------------------------------------------
%    Description: Processing is based on Holland, JASA 114 Oct 2003 
%                 pp. 1861-1873
%                 the data is sampled at 0.042 ms which results in a 
%                 f_ny = 12 kHz. The Boomer has a band up to
%                 10 kHz.
%     input:
%    output:
% ---------------------------------------------------------------

function [] = process_r();

set(0, 'DefaultFigurePaperPosition', [0 0 6 7.5]);
format('long','g');
fig = 1;

%-----------------------------------------------------------------------
% Switches:
%-----------------------------------------------------------------------

i_fgauss  = 0;  % Gaussian freq ave (Harrison range ave paper)
i_resamp  = 2;  % Resample to even angle spacing:
                % 0 = keep orig. uneven spacing
                % 1 = resample onto even grid by linear interp.
                % 2 = resample by averaging bins

i_syn     = 0;	% Are we running a Simulation? (then add noise)
i_range   = 0;  % Rave R as function of frequency and RANGE
i_linint  = 0;	% Use linear interpolation instead
i_verbose = 0;
i_spher   = 0;	% Compare to spher forward model?

niter     = 4;	% Number of iterations for SL polynom fit
snr_cut   = 6;	% Cut off SNR level in dB

r_tol = 0.001;  % Raytracing accuracy (m)
sd    = 0.05;	% RMS value for trace 133 of Charles' data
		% for SNR 12

%-----------------------------------------------------------------------
% Read input file:
%-----------------------------------------------------------------------

[trial,fileroot,fileroot2,fileext,zs,F,nang_re,freq,frbw,freq_low,...
 freq_up,cnmo,t0,rmx,svpfile] = lfa1998108141219_3_18_2lay_geo();

%
% 1/3 oct bands:
%
%freq_low = freq.*2.^(-1/6);
%freq_up  = freq.*2.^(1/6);

%-----------------------------------------------------------------------
% Plot stuff
%-----------------------------------------------------------------------
nx = 4;
ny = 4;
xim = 0.03;
yim = 0.28/ny;
xymarg = [0.04 0.04 0.04 0.14];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

rawfile  = strcat(fileroot,'.mat');

%-----------------------------------------------------------------------
% Output files:
%-----------------------------------------------------------------------

% Main R output:
outfile  = strcat(fileroot2,fileext,'.mat');
outfilea = strcat(fileroot2,fileext,'.txt');

%  Raw data plots:
tr_file  = strcat(fileroot2,fileext,'_trace.eps');
red_file = strcat(fileroot2,fileext,'_red.eps');

%  Bottom loss plots:
r_file   = strcat(fileroot2,fileext,'_r.eps');
ovs_file = strcat(fileroot2,fileext,'_ovs.eps');

%  Log parameters:
log_file = strcat(fileroot2,fileext,'.log');
log_id = fopen(log_file,'w');

fprintf(log_id,'cnmo: %10.4f\n',cnmo);
fprintf(log_id,'t0: %10.6f\n',t0);
fprintf(log_id,'max range: %10.4f\n',rmx);
fclose(log_id);


%-----------------------------------------------------------------------
% Data to load:
%-----------------------------------------------------------------------

svp = load(svpfile);        % SVP file (cannot handle constant vel.)
load(rawfile);              % Raw data
load('s04_source.mat');     % Source level data from Site04
x.c = svp(2,end);

if i_spher == 1
    load(spherfile);		% Plane replica data fav
    B_fav = F;
end

nz = length(svp);
nr = length(x.r);
nt = length(x.t_s);
nf = length(freq);
c  = svp(:,2)';
z  = svp(:,1)';

%if i_verbose == 1
%    figure(fig)
%    fig=fig+1;
%    box on;
%    plot(c,z);
%    set(gca,'YDir','reverse');
%    title('SVP');xlabel('sound-speed (m/s)');ylabel('depth (m)')
%end

if (min(x.r) < 0)
    x.r  = -1 * x.r;	  %  shot ranges
end
x.fout = 1/x.t_ax(2); %  sampling frequency


%-----------------------------------------------------------------------
% limit the range in case windowing is not clean:
%-----------------------------------------------------------------------

x.r = x.r(find(x.r<rmx));
F(1).r = x.r;
nr = length(x.r);


%-----------------------------------------------------------------------
%  Rayfast computes traveltime (t), angles (ang); pathlength to receiver 
%  along different paths (len), pathlength to x_d for direct (len2) and
%  surface reflected paths (len3), distance x_d (xd). Status outputs 
%  success of the raytracer.
%  Set paths to get arrivals 1) direct
%                            2) surface reflected
%                            3) bottom reflected
%                            4) surface reflected bottom reflected
%-----------------------------------------------------------------------

figure(fig);
fig = fig+1;
hold on; box on;
plot(0,x.zr,'ok');
plot([0 rmx],[z(end) z(end)],'-k');
plot([0 rmx],[x.zr x.zr],':k');
set(gca,'YDir','reverse');
set(gca,'YLim',[0 z(end)+10]);
set(gca,'XLim',[0 rmx]);

paths = [1 1 1 1 0 0];

ieven = 1;
for ir = 1:nr	% loop over shots (ranges)

    [t(:,ir),ang(:,ir),len(:,ir),len2(ir),len3(ir),xd(ir),theta(ir),... 
    status(:,ir)] = rayfast(paths,z,c,zs,x.zr,x.r(ir),r_tol);

%    if ir == 100%ieven
%        plot([0 x.r(ir)],[x.zr zs],'--k');
%        plot([x.r(ir)-xd(ir) x.r(ir)],[x.zr zs],'--k');
%        ieven = ieven +4;
%    end
%    x.r(ir),theta(ir)

end		% end loop over shots
 
%return;
%----------------------------------------------------------------------
%  Compute attenuation after Jensen et al [dB/km]
%----------------------------------------------------------------------

alpha = atten_fit(freq);
%alpha = alpha/1000;		% in dB/m
alpha = alpha .* 0;

%----------------------------------------------------------------------
%  Compute transmission coefficients
%----------------------------------------------------------------------

for ir = 1:nr		% loop over shots
  for ifr = 1:nf	% loop over frequencies

    g01(ir,ifr) = 1/len(3,ir)*exp(-alpha(ifr)*len(3,ir)); % bottom refl. path
    g02(ir,ifr) = 1/len(4,ir)*exp(-alpha(ifr)*len(4,ir)); % surf. bottom refl.
   
%   These next ones need to be evaluated at x_d!

    gd1(ir,ifr) = 1/len2(ir)*exp(-alpha(ifr)*len2(ir)); %  dir. path
    gd2(ir,ifr) = 1/len3(ir)*exp(-alpha(ifr)*len3(ir)); % surf. refl. path

  end			% end loop over frequencies
end			% end loop over shots

%-----------------------------------------------------------------------
%  Try to account for surface reflection:
%  w = g1*w1 + g2*w2 approx. g * (w1 + w2)
%  w:  total field
%  w1: direct wave		g1: transmiss. coeff direct path
%  w2: surface reflected wave	g2: transmiss. coeff surf. refl. path
%-----------------------------------------------------------------------

%g0 = (g01+g02)/2; gd = (gd1+gd2)/2;
g0 = g01; gd = gd1;
g = gd./g0;

if i_verbose == 1
    figure(fig);
    fig = fig+1;
    box on; hold on;
    title ('Transmission Coefficients');
    xlabel('Range (m)');ylabel('\gamma');
    plot(x.r,gd(:,5),'--b');
    plot(x.r,g0(:,5),'--r');
%    plot(x.r,g(:,5),'--k');
    legend('gd1','g0');

    figure(fig);
    fig = fig+1;
    box on; hold on;
    title ('Path Lengths');
    xlabel('Range (m)');ylabel('path length');
    plot(x.r,len2,'-b',x.r,len3,'--b');
    plot(x.r,len(1,:),'-r',x.r,len(2,:),'--r');
    plot(x.r,len(3,:),'-k',x.r,len(4,:),'--k');
    legend('len2','len3','len(1,:)','len(2,:)','len(3,:)','len(4,:)');
end

%-----------------------------------------------------------------------
%  Compute windows to window out arrivals and plot windows to check:
%-----------------------------------------------------------------------

tsq_lim=[0.02  0.4];
t_lim=[0.02 0.6];

for iw = 1:length(cnmo)

    tnmo(iw,:) = t0(iw) *(1 + (abs(x.r)/(cnmo(iw)*t0(iw))).^2).^0.5;

end

tsq2 = sqrt(tnmo.^2-repmat(x.r.^2/x.c^2,4,1));	% Reduced time

%-----------------------------------------------------------------------
%  Plot data
%-----------------------------------------------------------------------

[figh_tr] = raw_plot(x, t_lim, 2e0,tnmo,cnmo,fig);
saveas(figh_tr,tr_file,'epsc2');
fig = fig+1;

if(trial == 1)
    [figh_red] = reduce_sqrt(x,tsq_lim,2e0,tnmo,tsq2,cnmo,fig);
    saveas(figh_red,red_file,'epsc2');
    fig = fig+1;
    return;
end

%-----------------------------------------------------------------------
%  Window out data:
%-----------------------------------------------------------------------

% this zero pads all windows to same length...
i1 = find(tnmo(1,1) <= x.t_ax & x.t_ax <= tnmo(2,1));
i2 = find(tnmo(3,1) <= x.t_ax & x.t_ax <= tnmo(4,1));
i3 = find(0 < x.t_ax & x.t_ax <= tnmo(1,end)-0.01);
ii = max(max(length(i1),length(i2)),length(i3));
if(ii<512)
    ii=512;
elseif(ii<1024)
    ii=1024;
else
    ii=8192;
end

win1.dat = zeros(nr,ii);win1.t = zeros(nr,ii);
win2.dat = zeros(nr,ii);win2.t = zeros(nr,ii);
win3.dat = zeros(nr,ii);

%if i_verbose == 1
%   figure(fig);box on;hold on;
%   fig=fig+1;
%   xlabel('Range (m)');ylabel('Time re trigger (s)');
%   set(gca,'YDir','reverse');title ('Traces');
%end

fact = 2e2;

for ir = 1:nr	% Loop over shots

%  Window 1:
    i1 = find(tnmo(1,ir) <= x.t_ax & x.t_ax <= tnmo(2,ir));
    win1.dat(ir,1:length(i1)) = x.t_s(ir,i1);
    %  Remove mean of traces (eliminate drift)
    win1.dat(ir,:) = win1.dat(ir,:) - mean(win1.dat(ir,:),2);
    
    if(i1(1)+ii-1 <= length(x.t_ax))
        win1.t(ir,1:ii) = x.t_ax(i1(1):i1(1)+ii-1);
    else
        win1.t(ir,:) = x.t_ax(i1(1))+mean(diff(x.t_ax))*[0:ii-1];
    end
    win1.dt(ir)              = win1.t(ir,2) - win1.t(ir,1);
    win1.dT(ir)              = win1.t(ir,end) - win1.t(ir,1);

%    if i_verbose == 1
%        win1.datp(ir,:) = fact*win1.dat(ir,:) + x.r(ir);
%        plot(win1.datp(ir,:), win1.t(ir,:), 'b');
%    end

%  Window 2:
    i2 = find(tnmo(3,ir) <= x.t_ax & x.t_ax <= tnmo(4,ir));
    win2.dat(ir,1:length(i2)) = x.t_s(ir,i2);
    %  Remove mean of traces (eliminate drift)
    win2.dat(ir,:) = win2.dat(ir,:) - mean(win2.dat(ir,:),2);
    if(i2(1)+ii-1 <= length(x.t_ax))
        win2.t(ir,1:ii) = x.t_ax(i2(1):i2(1)+ii-1);
    else
        win2.t(ir,:) = x.t_ax(i2(1))+mean(diff(x.t_ax))*[0:ii-1];
    end
    win2.dt(ir)              = win2.t(ir,2) - win2.t(ir,1);
    win2.dT(ir)              = win2.t(ir,end) - win2.t(ir,1);

%    if i_verbose == 1
%        win2.datp(ir,:) = fact*win2.dat(ir,:) + x.r(ir);
%        plot(win2.datp(ir,:), win2.t(ir,:), 'b');
%    end

%  Noise window (t=0 to first break for each offset):
    i3 = find(0 <= x.t_ax & x.t_ax <= tnmo(1,ir)-0.01);
    win3.dat(ir,1:length(i3)) = x.t_s(ir,i3);
%    win3.datp(ir,:) = fact*win3.dat(ir,:) + x.r(ir);
%    plot(win3.datp(ir,:), x.t_ax(1:length(win3.dat(ir,:))), 'b');

    win1.correct(ir) = 10*log10(length(i1)/length(i3));
    win2.correct(ir) = 10*log10(length(i2)/length(i3));
    win3.correct(ir) = 10*log10(length(i2)/length(i1));
    
end		% End loop over shots

%-----------------------------------------------------------------------
%  Compute FFTs of Windows as function of angle
%-----------------------------------------------------------------------

for ir = 1:nr

    win1.DAT(ir,:) = (fft(win1.dat(ir,:),ii));
    win2.DAT(ir,:) = (fft(win2.dat(ir,:),ii));
    win3.DAT(ir,:) = (fft(win3.dat(ir,:),ii));
    
    qd(ir,:) = abs(win1.DAT(ir,:));
    pr(ir,:) = abs(win2.DAT(ir,:));
    tmp3(ir,:) = abs(win3.DAT(ir,:));

end;


%-----------------------------------------------------------------------
%  Compute FFT frequencies (same for both windows since zero padded)
%-----------------------------------------------------------------------

f  = 1/mean(win1.dT) * (0:ii/2);
f_ny = ii/(2*mean(win1.dT));
win1.df = 1/mean(win1.dT);
win1.f = f;

%-----------------------------------------------------------------------
%  Calculate SNR:
%-----------------------------------------------------------------------

for ifr = 1:nf
    
    idx = find(f <= freq_up(ifr) & f >= freq_low(ifr));
    win1.oct(:,ifr) = 10*log10(2/x.fout^2*sum(qd(:,idx),2)/length(idx));
    win2.oct(:,ifr) = 10*log10(2/x.fout^2*sum(pr(:,idx),2)/length(idx));
    win3.oct(:,ifr) = 10*log10(2/x.fout^2*sum(tmp3(:,idx),2)/length(idx));

end

snr1 = win1.oct-win3.oct-repmat(win1.correct',1,nf);
snr2 = win2.oct-win3.oct-repmat(win2.correct',1,nf);

%
%  Average SNR with sliding window:
%
nsmooth = 7;
for ir = 4:nr-3
    
    snr1(ir,:) = sum(snr1(ir-3:ir+3,:),1)/nsmooth;
    snr2(ir,:) = sum(snr2(ir-3:ir+3,:),1)/nsmooth;
    
end

clear tmp3;

%-----------------------------------------------------------------------
%  RMS Frequency ave over bands (specified above)
%-----------------------------------------------------------------------

if(i_fgauss == 0)
    disp('Starting band average:')
    clear idx;
    for ifr = 1:nf
    
        idx = find(f <= freq_up(ifr) & f >= freq_low(ifr));
%        qd_ave(:,ifr) = sqrt(mean(qd(:,idx).^2,2));
%        pr_ave(:,ifr) = sqrt(mean(pr(:,idx).^2,2));

        freq(ifr)
        sen(ifr) = tsystsen(x.vga(1),x.pre,300,freq(ifr));
        qd_ave(:,ifr) = 1/x.fout^2*win1.df/(freq_up(ifr)-freq_low(ifr)) ...
                        *sum(qd(:,idx).*conj(qd(:,idx)),2) ...
                        *1./(10^(sen(ifr)/10));
        pr_ave(:,ifr) = 1/x.fout^2*win1.df/(freq_up(ifr)-freq_low(ifr)) ...
                        *sum(pr(:,idx).*conj(pr(:,idx)),2) ...
                        *1./(10^(sen(ifr)/10));

    end
    disp('Band average done.')
else
    disp('Starting Gaussian average (Harrison):')
    df = mean(diff(f));
    nf2 = find(freq(end)<f);
    nf2 = nf2(1) + floor(nf2(1)*.5);
    for ir = 1:nr
    for ifr = 2:nf2

        f0 = f(ifr);
        qd_ave_tmp(ir,ifr) = sum(qd(ir,2:nf2) .* ... 
                             exp(-(f(2:nf2)-f0).^2/...
                             (frbw*f0)^2)*df)/...
                             sum(exp(-(f-f0).^2/...
                             (frbw*f0)^2)*df); 
        pr_ave_tmp(ir,ifr) = sum(pr(ir,2:nf2) .* ... 
                             exp(-(f(2:nf2)-f0).^2/...
                             (frbw*f0)^2)*df)/...
                             sum(exp(-(f-f0).^2/...
                             (frbw*f0)^2)*df); 

    end
    end
    clear idx;
    for ifr = 1:nf
        idx = find((f>(freq(ifr)-df)) & (f<(freq(ifr)+df)));
        if(abs(freq(ifr)-f(idx(1))) < abs(freq(ifr)-f(idx(2))))
            idx = idx(1);
        else
            idx = idx(2);
        end
        qd_ave(:,ifr) = qd_ave_tmp(:,idx);
        pr_ave(:,ifr) = pr_ave_tmp(:,idx);
    end
    disp('Gaussian average done.')
end

%if i_verbose == 1
%    figure(fig);box on;hold on;
%    fig=fig+1;
%    plot(f,pr(100,1:ii/2+1)); title('Magnitude');
%
%    figure(fig);box on;hold on;
%    fig=fig+1;
%    plot(unwrap(Angle(qd(100,:)))); title('Phase');
%end

%-----------------------------------------------------------------------
%  Interpolate q_d (polynomial fit or linear interpolation)
%-----------------------------------------------------------------------

figure(fig);box on;hold on;
fig=fig+1;
%polygrade = [6 6 6 6 6 6 6 6 6];
%polygrade = [9 9 9 9 9 9 9 9 9];
%polygrade = [3 3 3 3 3 3 3 3 3];
polygrade = 6;
for iter = 1:niter
    for ifr = 1:nf
    if i_linint == 0

        [pol,S,mu] = polyfit(s04_ang',s04_qd_ave(:,ifr),polygrade);
        [qd_ave_new(:,ifr),delta(:,ifr)]  = polyval(pol,ang(1,:),S,mu);

        [pol,S,mu] = polyfit(x.r',qd_ave_new(:,ifr),polygrade);
        [qd_avei(:,ifr),delta(:,ifr)]  = polyval(pol,xd,S,mu);
        [qd_avei2(:,ifr),delta2(:,ifr)] = polyval(pol,x.r',S,mu);

    end

    qd_diff = qd_ave(:,ifr)-qd_avei2(:,ifr);
    stdev(ifr) = std(qd_diff);
    mxdev(ifr) = max(qd_diff);

    if(ifr<=16)
    subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);
    box on;hold on;
%    plot(x.r,qd_ave(:,ifr),'.',xd,qd_avei(:,ifr),'+k');
    plot(ang(1,:),log10(qd_ave(:,ifr)),'.k');
    plot(s04_ang,log10(s04_qd_ave(:,ifr)),'--');
    set(gca,'YLim',[4 8]);
    set(gca,'XLim',[5 40]);
    end

    clear idx;
    idx = find(abs(qd_diff) > 2*stdev(ifr));
    pr_ave(idx,ifr) = NaN;	% 2 sigma screening
    idx_snr1 = find(snr1(:,ifr) < snr_cut);
    idx_snr2 = find(snr2(:,ifr) < snr_cut);
    pr_ave(idx_snr1,ifr) = NaN;	% SNR screening
    pr_ave(idx_snr2,ifr) = NaN;	% SNR screening
    
    end
end

%F(1).xd = xd;
%F(1).qd_avei = qd_avei;
%-----------------------------------------------------------------------
%  Calculate R
%-----------------------------------------------------------------------

R(:,:) = (abs(pr_ave)./qd_avei) .* g;

%R2(:,:) = abs(pr_ave ./ qd_ave) .* (qd_avei2 ./ qd_avei) .* g;

theta = fliplr(theta);
R  = flipud(R);
%R2  = flipud(R2);
for ifr = 1:nf
    idx = find(isnan(R(:,ifr))~=1);
    idx2 = find(isnan(R(:,ifr))==1);
    Rnan(ifr).idx = idx;
    Rnan(ifr).idx2 = idx2;
    Rnan(ifr).theta = theta(idx);
    Rnan(ifr).dat = R(idx,ifr);
%    Rnan(ifr).dat2 = R2(idx,ifr);
end

%-----------------------------------------------------------------------
%  Resample onto even grid
%-----------------------------------------------------------------------
if(i_resamp > 0)

    disp('RESAMPLE ONTO EVEN GRID')
    theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];

    Rex = ones(length(theta_even),nf);
    theta_diff = diff(theta);
    for iang = 2:length(theta)-1
        if(theta_diff(iang)>2*theta_diff(iang-1))
            idx = find((theta_even>=theta(iang)) & ...
                       (theta_even<=theta(iang+1)));
            %
            % This Rex is for data gaps
            %
            Rex(idx(2:end-1),:) = 0;
        end
    end
    for ifr = 1:nf
    if(length(Rnan(ifr).idx2)>0)
    for iidx = 1:length(Rnan(ifr).idx2)
        if(Rnan(ifr).idx2(iidx)+1 <= length(theta) & Rnan(ifr).idx2(iidx)-1 >= 1)
            idx = find((theta_even>theta(Rnan(ifr).idx2(iidx)-1)) & ...
                      (theta_even<theta(Rnan(ifr).idx2(iidx)+1)));
        else
            if(Rnan(ifr).idx2(iidx)-1 < 1)
                idx = find((theta_even>theta(Rnan(ifr).idx2(iidx))) & ...
                          (theta_even<theta(Rnan(ifr).idx2(iidx)+1)));
            else
                idx = find((theta_even>theta(Rnan(ifr).idx2(iidx)-1)) & ...
                          (theta_even<theta(Rnan(ifr).idx2(iidx))));
            end
        end
        Rex(idx(1:end-1),ifr) = 0;
    end;end;end;

    clear delta;

    for ifr = 1:nf

        R_even(:,ifr)   = interp1(Rnan(ifr).theta,Rnan(ifr).dat,theta_even);
%        R_even2(:,ifr)   = interp1(Rnan(ifr).theta,Rnan(ifr).dat2,theta_even);

    end

  if(i_resamp == 2)
%
%  This will average the data where there is data and take the interpolation
%  from above where there was no data or poor data. interpolated regions
%  are excluded from the inverion through Rex
%

    disp('RESAMPLE ONTO EVEN GRID: BIN AVERAGE')
    theta_even = [theta(1):(theta(end)-theta(1))/nang_re:theta(end)];

    Rex = ones(length(theta_even),nf);
    theta_diff = diff(theta);
    for iang = 2:length(theta)-1
        if(theta_diff(iang)>2*theta_diff(iang-1))
            idx = find((theta_even>=theta(iang)) & ...
                       (theta_even<=theta(iang+1)));
            %
            % This Rex is for data gaps
            %
            Rex(idx(2:end-1),:) = 0;
        end
    end
    for ifr = 1:nf
    if(length(Rnan(ifr).idx2)>0)
    for iidx = 1:length(Rnan(ifr).idx2)
        if(Rnan(ifr).idx2(iidx)+1 <= length(theta) & Rnan(ifr).idx2(iidx)-1 >= 1)
            idx = find((theta_even>theta(Rnan(ifr).idx2(iidx)-1)) & ...
                      (theta_even<theta(Rnan(ifr).idx2(iidx)+1)));
        else
            if(Rnan(ifr).idx2(iidx)-1 < 1)
                idx = find((theta_even>theta(Rnan(ifr).idx2(iidx))) & ...
                          (theta_even<theta(Rnan(ifr).idx2(iidx)+1)));
            else
                idx = find((theta_even>theta(Rnan(ifr).idx2(iidx)-1)) & ...
                          (theta_even<theta(Rnan(ifr).idx2(iidx))));
            end
        end
        Rex(idx(1:end-1),ifr) = 0;
    end;end;end;

    clear delta;

    for ifr = 1:nf
       for iang = 2:length(theta_even)-1
           idx = find(Rnan(ifr).theta > theta_even(iang-1) & ...
                      Rnan(ifr).theta < theta_even(iang+1));
           if(idx)
              R_even(iang,ifr)   = mean(Rnan(ifr).dat(idx));
%              R_even2(iang,ifr)  = mean(Rnan(ifr).dat2(idx));
           else
              Rex(iang,ifr) = 0;
           end
       end
       idx = find(Rnan(ifr).theta > theta_even(1) & ...
                  Rnan(ifr).theta < theta_even(2));
       if(idx)
          R_even(1,ifr)   = mean(Rnan(ifr).dat(idx));
%          R_even2(1,ifr)  = mean(Rnan(ifr).dat2(idx));
       else
          Rex(1,ifr) = 0;
       end
       idx = find(Rnan(ifr).theta > theta_even(end-1) & ...
                  Rnan(ifr).theta < theta_even(end));
       if(idx)
          R_even(nang_re+1,ifr)   = mean(Rnan(ifr).dat(idx));
%          R_even2(nang_re+1,ifr)  = mean(Rnan(ifr).dat2(idx));
       else
          Rex(iang,ifr) = 0;
       end
    end

  end

  clear theta;
  clear R;
  R = R_even;
%  R2 = R_even2;
  theta = theta_even;
else
    theta_even = theta;
    R_even = R;
    Rex = ones(length(theta_even),nf);
    for ifr = 1:nf

        Rex(Rnan(ifr).idx2,ifr) = 0;
        R_even(Rnan(ifr).idx2,ifr) = 0;

    end
    R = R_even;
end

%-----------------------------------------------------------------------
%  Add noise to synthetic data
%-----------------------------------------------------------------------

if(i_syn == 1)
  for ifr = 1:nf	% Loop over shots
      randn('state',sum(100*clock));

      noise = sd * randn(length(R(:,ifr)),1);
      R(:,ifr) = R(:,ifr) + noise;

  end
end

%    correct corrects for window width (???)
BL = -20*log10(abs(R));
%BL = -20*log10(abs(R))-repmat(win3.correct',1,nf);
%BL2 = -20*log10(abs(R2));
%BL2 = -20*log10(abs(R2))-repmat(win3.correct',1,nf);

%-----------------------------------------------------------------------
%  Plot end result:
%-----------------------------------------------------------------------

ii = 1;
for ifr = 1:nf

    if(i_verbose == 1) 
        figh_r = figure(fig);
        if(ifr<=16)
            sub9(ifr) = subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);
            box on;hold on;
            plot(theta,R(:,ifr),'.b');
%            plot(Rnan(ifr).theta,0.5*ones(size(Rnan(ifr).theta)),'or');
            plot(theta_even,Rex(:,ifr),'or');
            plot(theta_even,R_even(:,ifr),'.-k');
            set(gca,'XLim',[10 40],'YLim',[0 1.1]);
            set(gca,'XTick',[10 20 30 40],'XTickLabel',{'10';'20';'30';'40'});
            set(gca,'YTick',[0 0.5 1],'YTickLabel',{'0';'0.5';'1'});
        end;
    end

    fig9 = figure(fig+1);
    if(ifr<=16)
        sub9(ifr) = subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);
    box on;hold on;
%    plot(theta,BL(:,ifr),'.');
    plot(theta,R(:,ifr),'.');
    plot(theta_even,R_even(:,ifr),'.-k');
    legend('Eq. 4','Even');
    if i_spher == 1
        plot(B_fav(ifr).ang,B_fav(ifr).dat,'--k');
%        plot(B_fav(ifr).ang,B_fav(ifr).plane,':r');
        legend('Eq. 4','spher');
    end
    xlabel('Angle');ylabel('BL')
    title('Bottom Loss');
%    set(gca,'XLim',[10 80],'YLim',[0 30]);
    set(gca,'XLim',[10 40],'YLim',[0 1]);
    set(gca,'XTick',[10 20 30 40],'XTickLabel',{'10';'20';'30';'40'});
    end;
    
    fig12 = figure(fig+3);
    if(ifr<=16)
        sub9(ifr) = subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);
    box on;hold on;
    plot(theta,BL(:,ifr),'-');
    if i_spher == 1
        plot(B_fav(ifr).ang,B_fav(ifr).dat,'--k');
%        plot(B_fav(ifr).ang,B_fav(ifr).plane,':r');
        legend('Eq. 4','OASR','spher');
    end
    xlabel('Angle');ylabel('BL')
    title('Bottom Loss');
    set(gca,'XLim',[10 40],'YLim',[0 30]);
    set(gca,'XTick',[10 20 30 40],'XTickLabel',{'10';'20';'30';'40'});
    end;

    fig13 = figure(fig+4);
    if(ifr<=16)
        sub9(ifr) = subplot('Position',[loc(1,ifr) loc(2,ifr) spw sph]);
    box on;hold on;
    plot(snr1(:,ifr),'--b');
    plot(snr2(:,ifr),'--r');
    plot([1,length(snr1(:,ifr))],[snr_cut,snr_cut],'--k');
    end;
    
end
if(i_verbose == 1)
    label_subplot(figh_r,sub9,'Angle (Deg.)','|R|',...
    'Reflection Coefficient',nf,nrp,nc,spw,sph,loc(1,k),loc(2,k))
    saveas(figh_r,r_file,'epsc2');
    %exportfig(figh_r,oasr_seis_file);
%    saveas(figh_r,ovs_file,'epsc2');
end
saveas(fig12,ovs_file,'epsc2');

%-----------------------------------------------------------------------
%   Prepare data for saving
%-----------------------------------------------------------------------


F(1).freq = freq;
clear tmp;
tmp = BL(1:end,:);
BL = tmp;
clear tmp;
tmp = theta(1:end);
theta = tmp;
for ifr = 1:nf;
%    F(ifr).dat = BL(:,ifr);
    F(ifr).dat = R(:,ifr);
    F(ifr).ang = theta;
end

F(1).nmod = length(F(1).maxlim);

minlim = F(1).minlim;
maxlim = F(1).maxlim;
mtrue  = F(1).mtrue;
F(1).Rex = Rex;

for i = 1:nf
    F(i).dat(find(isnan(R(:,i)))) = 0;
    F(1).Rex(find(isnan(R(:,i))),i) = 0;
end

save(outfile,'F');
save(outfilea,'minlim','-ascii');
save(outfilea,'maxlim','-ascii','-append');
save(outfilea,'mtrue','-ascii','-append');
for i = 1:nf

    tmp = F(i).dat';
    save(outfilea,'tmp','-ascii','-append');

end
ang = F(1).ang(1:size(Rex,1));
save(outfilea,'ang','-ascii','-append');

for i = 1:nf

    tmp = F(1).Rex(:,i)';
    save(outfilea,'tmp','-ascii','-append');

end


%range = fliplr(x.r(1:length(Rex)));
%save(outfilea,'range','-ascii','-append');

return;

%=============================================================================

function [figh] = reduce_sqrt(x, tsq_lim, fact, tnmo, tsq2, cnmo, fig);

figh = figure(fig);
for i = 1:length(x.r)
     tsq = sqrt(x.t_ax.^2-x.r(i)^2/x.c^2);
     i1 = find(tsq_lim(1) <= tsq & tsq <= tsq_lim(2));
     line('Color', 'k', ...
          'XData', fact*x.t_s(i, i1)+x.r(i), ...
          'YData', tsq(i1), 'LineWidth',.1)
end
xlabel('Range (m)','FontSize',14); ylabel('Reduced Time (s)','FontSize',14);
axis([min(x.r) max(x.r) tsq_lim]); box on;
set(gca,'YDir','reverse','FontSize',14);% title ('Traces');
hold on;

%for iw = 1:length(cnmo)
%
%  plot(x.r,tsq2(iw,:),'k');
%
%end
  plot(x.r,tsq2(1,:),'k');
  plot(x.r,tsq2(2,:),'k');
  plot(x.r,tsq2(3,:),'k');
  plot(x.r,tsq2(4,:),'--k');
%  plot(x.r,tsq2(5,:),'--k');

return;


%=============================================================================

function [figh] = raw_plot(x, t_lim, fact,tnmo,cnmo, fig)

figh = figure(fig);
hold on;box on;
i = find(t_lim(1) <= x.t_ax & x.t_ax <= t_lim(2));
i = i(1:4:end);
x1 = x.t_s(:, i);
for i1 = 1:2:length(x.r)
     x1(i1, :) = fact*x1(i1, :)+x.r(i1);
end
plot(x1, x.t_ax(i), 'k');
axis([min(x.r(1), x.r(end)),max(x.r(1), x.r(end)), x.t_ax(i(1)), x.t_ax(i(end))]);

xlabel('Range (m)');
ylabel('TWT (s)');
set(gca,'YDir','reverse');
set(gca,'XTick',[100 200 300 400 500 600 700]);
%title ('Traces');


for iw = 1:length(cnmo)-1

  plot(x.r,tnmo(iw,:),'k');

end
plot(x.r,tnmo(end,:),'r');
%load bla2.mat
%plot(x.r,tnmo(end,:),'r');
%load bla4.mat
%plot(x.r,tnmo(end,:),'r');
%save bla tnmo;
%---------------------------------------------------------
% Interactive change:
%---------------------------------------------------------
i_inter = 0;

if i_inter == 1
    n = 4;
    for ip = 1:2*n

        [xi(ip),yi(ip)] = ginput(1);
        plot(xi(ip),yi(ip),'ro');

    end
    for iw = 1:n
        
        dr = xi(iw*2) - xi((iw*2)-1)
        c(iw) = sqrt(dr^2/(yi(iw*2)^2-yi((iw*2)-1)^2))
        tnmo2(iw,:) = yi((iw-1)*2+1) *(1 + (abs(x.r-x.r(1))/(c(iw)*...
                      yi((iw-1)*2+1))).^2).^0.5;
        plot(x.r,tnmo2(iw,:),'--r');

    end

end


return;


%=============================================================================
%  This is the end my fiend...
%  EOF
%
