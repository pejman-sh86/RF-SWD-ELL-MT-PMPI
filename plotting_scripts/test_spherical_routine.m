function [] = test_spherical_routine();

load geo_sin;
%f = freq;
f = [1000];
cw = 1511;
z_h = 120;

%
% Compute OASES reference:
%

  freq_low = f.*2.^(-1/6);
  freq_up  = f.*2.^(1/6);
  freq_ave = [freq_low:freq_up];

  oases_file = 'input.trc'
  [R_oases,ang_oases,f_oases] = read_oases_r(oases_file);
  BL_oases = -20*log10(abs(R_oases));
  clear idx;
  for ifr = 1:length(f)
      idx = find(f_oases <= freq_up(ifr) & f_oases >= freq_low(ifr));
      BL_oases_ave(:,ifr) = sum(BL_oases(:,idx),2)/length(idx);
      R_oases_ave(:,ifr) = sum(abs(R_oases(:,idx)),2)/length(idx);
  end


%
%  Plot with original routine:
%

thd_orig = 8:.005:78;

ref_orig=ref_nlay3(thd_orig,geo_sin,f);
[thd_deg_orig,spher_r_num_orig]=spherical_ref_Brek_num_orig(thd_orig,...
                     cw,ref_orig,f,z_h);
ref_orig         = -20*log10(abs(ref_orig));
spher_r_num_orig = -20*log10(abs(spher_r_num_orig));

figure;hold on;box on;
plot(thd_orig,ref_orig(1,:),'--k');
plot(thd_deg_orig,spher_r_num_orig(1,:),'-b');
legend('orig plane','orig spher');
title('Charles orig routine');


%
%  Plot with Jan's simple routine:
%

thd_un = F(5).ang;
ref_un=ref_nlay3(thd_un,geo_sin,f);
[thd_deg_un,spher_r_num_un]=spherical_ref_Brek_num(thd_un,...
                     cw,ref_un,f,z_h);
ref_un         = -20*log10(abs(ref_un));
spher_r_num_un = -20*log10(abs(spher_r_num_un));

thd_ev = thd_orig;
ref_ev=ref_nlay3(thd_ev,geo_sin,f);
[thd_deg_ev,spher_r_num_ev]=spherical_ref_Brek_num(thd_ev,...
                     cw,ref_ev,f,z_h);
ref_ev         = -20*log10(abs(ref_ev));
spher_r_num_ev = -20*log10(abs(spher_r_num_ev));

figure;hold on;box on;
plot(thd_un,ref_un(1,:),'-.k');
plot(thd_deg_un,spher_r_num_un(1,:),'.k');
plot(thd_ev,ref_ev(1,:),'-.r');
plot(thd_deg_ev,spher_r_num_ev(1,:),'.r');
legend('uneven plane','uneven spher','even plane','even spher')
title('Jans routine; also for uneven spacing; sampling artifacts')

%
%  Plot with Jan's better derivs routine:
%

thd2 = 8:.1:78;
d2 = 0.01;
%d = (thd(10) - thd(9))/10;
for iang = 1:length(thd2)
  thd2m(iang) = thd2(iang)-d2;
  thd2p(iang) = thd2(iang)+d2;
end
ref2=ref_nlay3(thd2,geo_sin,f);
ref2m=ref_nlay3(thd2m,geo_sin,f);
ref2p=ref_nlay3(thd2p,geo_sin,f);
[thd_deg2,spher_r_num2]=spherical_ref_Brek_num_b(thd2,...
                     thd2m,thd2p,cw,ref2,ref2m,ref2p,f,z_h);
ref2         = -20*log10(abs(ref2));
ref2m        = -20*log10(abs(ref2m));
ref2p        = -20*log10(abs(ref2p));
spher_r_num2 = -20*log10(abs(spher_r_num2));

figure;hold on;box on;
plot(thd2,ref2(1,:),'--k');
plot(thd_deg2,spher_r_num2(1,:),'.k');
plot(thd_deg_orig,spher_r_num_orig(1,:),'-b');
legend('plane better','spher better','reference');
title('Better derivs routine fine angles even spacing');


%
%  Plot with Jan's better derivs routine
%  with coarse uneven angle spacing  :
%

thd3 = F(5).ang;
%thd3 = 8:1:78;
d3 = .01;
for iang = 1:length(thd3)
  thd3m(iang) = thd3(iang)-d3;
  thd3p(iang) = thd3(iang)+d3;
end
ref3=ref_nlay3(thd3,geo_sin,f);
ref3m=ref_nlay3(thd3m,geo_sin,f);
ref3p=ref_nlay3(thd3p,geo_sin,f);
[thd_deg3,spher_r_num3]=spherical_ref_Brek_num_b(thd3,...
                     thd3m,thd3p,cw,ref3,ref3m,ref3p,f,z_h);


ref3         = -20*log10(abs(ref3));
ref3m        = -20*log10(abs(ref3m));
ref3p        = -20*log10(abs(ref3p));
spher_r_num3 = -20*log10(abs(spher_r_num3));

figure;hold on;box on;
plot(thd3,ref3(1,:),'--k');
plot(thd_deg_orig,spher_r_num_orig(1,:),'-b');
plot(ang_oases,BL_oases_ave(:,ifr),'--r');
legend('plane','steepest descent','OASR');
title('Comparing plane wave refl coeff with steepest descent and OASR');

figure;hold on;box on;
plot(thd_deg_un,spher_r_num_un(1,:),'.r');
plot(thd_deg3,spher_r_num3(1,:),'.k');
plot(thd_deg_orig,spher_r_num_orig(1,:),'-b');
legend('uneven simple','uneven better','reference');
title('Comparing steepest descent refl coeff.');

return;
