function [] = rearrange_oasr();

  oasrfile = 'input_dag4oasr.trc'
  repoasrfile = 'dag4_1oasr.mat'

%  freq     = [500 630 800 1000 1250 1600 2000 2511 3162 4000];
%  freq     = [315 400 500 630 800 1000 1250 1600 2000];
  freq     = [50 63 80 100 125 160 200 250 315];
  freq_low = freq.*2.^(-1/6);
  freq_up  = freq.*2.^(1/6);

  nf = length(freq);
  [R_oases,ang_oases,f_oases] = read_oases_r(oasrfile);
  BL_oases = -20*log10(abs(R_oases));
  clear idx;
  for ifr = 1:nf
      idx = find(f_oases <= freq_up(ifr) & f_oases >= freq_low(ifr));
      R_oases_ave(:,ifr)  = sqrt(mean(abs(R_oases(:,idx)).^2,2));
      BL_oases_ave(:,ifr) = -20*log10(R_oases_ave(:,ifr));
  end

  save(repoasrfile,'R_oases_ave','BL_oases_ave','BL_oases','ang_oases','f_oases');

return;
