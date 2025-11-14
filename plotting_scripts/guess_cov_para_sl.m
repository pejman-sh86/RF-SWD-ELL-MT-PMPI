load hol_ida;
A = load('syn_cov');
A = getfield(A,'F');

nfreq = length(F(1).freq);
nfreq = 1;

%F(1).ang = F(6).ang;
%sigma_start = [.1];
%lambda_start = [2000];
%lambda_start = [0.6];
sigma_start = [1 1 1 1 1 1 1 1];
lambda_start = [.9 .9 .9 .9 .9 .9 .9 .9];
beta = [.9 .9 .9 .9 .9 .9 .9 .9];
%load hol_syn_map.mat
%sigma_start2 = xmap(1:8);
%lambda_start2 = xmap(9:16);
%sigma_start = [2.2339      2.2402      2.2724      1.9648      2.2574 ... 
%               2.062      1.7534      1.8479];
%lambda_start = [0.86869     0.57247     0.36647     0.73758     0.62455 ...
%                0.51178     0.57116     0.69485];


%  for ifreq=1:nfreq
%    F(ifreq).csave = A(ifreq).csave;
%    F(ifreq).c_inv = inv(F(ifreq).csave);
%  end;

  figure(1);
  for ifreq=1:nfreq
%    B(ifreq).csave = cov_est_par_a(F(ifreq).ang,sigma_start(ifreq),...
%                     lambda_start(ifreq));
%    B(ifreq).csaveb = cov_est_par_a(F(ifreq).ang,sigma_start2(ifreq),...
%                     lambda_start2(ifreq));
%    B(ifreq).csavec = cov_est_par_c(F(ifreq).ang,sigma_start(ifreq),...
%                     lambda_start(ifreq));
    B(ifreq).csave = cov_est_par_d(F(ifreq).ang,sigma_start(ifreq),...
                     lambda_start(ifreq),beta(ifreq));

    subplot(3,3,ifreq);
    hold on;
%    plot(B(ifreq).csave(1,:),'k-o');
    plot(B(ifreq).csave(round(length(F(ifreq).ang)/2),:),'k-o');
%    plot(B(ifreq).csaveb(round(length(F(ifreq).ang)/2),:),'r-o');
%    plot(B(ifreq).csavec(round(length(F(ifreq).ang)/2),:),'g-o');
%    plot(A(ifreq).csave(round(length(F(ifreq).ang)/2),:),'b-o');

    legend('true','map')
  end;

clear F
F=B;
save syn_covbmf F;


