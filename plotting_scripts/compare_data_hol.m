function [] = compare_data_hol()

load rep_15_fav1.mat;
Fjan = F;
load rep_15_fav2.mat;
Frep = F;
load hol_id2b.mat;
nf = length(F(1).freq);
% margins: left, right, middle, top, bottom, Plots per row and per column
nr = 3; nc = 3;
lm = .05; rm = .03; mm = .05; tm = .05; bm = .05;
[width,height,llh,llw] = calc_subplot(lm,rm,mm,tm,bm,nr,nc);

ii = 1;
for ifr = 1:nf

    if ifr <= 3
      sp = [llw(ii) llh(1) width height];
      ii = ii + 1;
    elseif ifr <= 6
      sp = [llw(ii) llh(2) width height];
      ii = ii + 1;
    else
      sp = [llw(ii) llh(3) width height];
      ii = ii + 1;
    end
    
    fig1 = figure(1);
    sub1(ifr) = subplot('Position',sp);box on;hold on;
    plot(F(ifr).ang,F(ifr).dat,'xb');
    plot(Frep(ifr).ang,Frep(ifr).dat,'--r');
    plot(Fjan(ifr).ang,Fjan(ifr).dat,'--k');
    xlabel('angle');ylabel('BL')
    title('Bottom Loss');legend('Hol','No Fav','Fav');
    set(gca,'XLim',[10 80],'YLim',[0 30]);
    
    
end

label_subplot(fig1,sub1,'Angle (Deg.)','|R|','Reflection Coefficient',nf,nr,nc,width,height,llh,llw)
%exportfig(fig9,oasr_seis_file);
%saveas(fig9,oasr_seis_file,'epsc2');
    
%save data tsq2 snr BL BL2;
return;
