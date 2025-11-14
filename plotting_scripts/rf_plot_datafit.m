function []=rf_plot_datafit_rf();
iar = 0;

%%
%%  Data fit and residual plots:
%%
dat_all=load('RIY_SELECT5_REG1.txt');
rep=load('RIY_SELECT5_REG1_rep.dat');
rep_all=repmat(rep,18,1);
if(iar == 1);dar=load('RIY_SELECT5_REG1_ar.dat');end;
figure();hold on;
for i=1:18;
  plot(dat_all(i,:));
end;
plot(rep,'r','LineWidth',2)

res_raw=dat_all-rep_all;
if(iar == 1);res_tot=dat_all-rep_all-dar;end;

figure();hold on;
for i=1:18;
  plot(res_raw(i,:));
end;
set(gca,'YLim',[-.1 .1]);

if(iar == 1)
  figure();hold on;
  for i=1:18;
    plot(res_tot(i,:));
  end;
  set(gca,'YLim',[-.1 .1]);

  figure();hold on;
  for i=1:18;
    plot(dar(i,:));
  end;
  set(gca,'YLim',[-.1 .1]);
end;

%%
%% Runstest of residuals
%%
disp('          Raw                Total     ');
disp('      H         p         H         p  ');
for i=1:18;
  [H1,p1]=runstest(res_raw(i,:));
  if(iar == 1)
    [H2,p2]=runstest(res_tot(i,:));
  else;
    H2=1; p2=0;
  end;
  disp([H1,p1,H2,p2]);
end;

return;
