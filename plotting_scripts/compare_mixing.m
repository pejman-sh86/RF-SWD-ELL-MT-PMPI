function []=compare_mixing();
set(0, 'DefaultFigurePaperPosition', [0 0 11 6]);

load ../test_PT2/x_s02_1_300_1600_mind_sample.mat
B=A;
pwd
load ../test_PT/x_s02_1_300_1600_mind_sample.mat
pwd
kmin = 3.75;
kmax = 8.75;


NC = 2;
fig1=figure();
nx = NC;
ny = 2;
xim = 0.01;
yim = 0.01;
xymarg = [0.08 0.01 0.01 0.1];
[loc,spw,sph] = get_loc(nx,ny,xim,yim,xymarg);

N = zeros(2*NC,1);
for i=1:NC;

  clear idx1 idx2;
  idx1=find(A(:,112)==i);
  idx2=find(B(:,112)==i);
  for j = 2:length(idx1);
    if(A(idx1(j),4)~=A(idx1(j-1),4));N(i) = N(i)+1;end;
  end;
  for j = 2:length(idx2);
    if(B(idx2(j),4)~=B(idx2(j-1),4));N(i+NC) = N(i+NC)+1;end;
  end;

  h1 = subplot('Position',[loc(1,i) loc(2,i) spw sph]);
  plot(A(idx1,4));
  set(gca,'XLim',[0 length(idx1)],'YLim',[kmin kmax],'FontSize',14);
  if(i == 1);ylabel('k');end;
  set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
  if(i ~= 1);set(gca,'YTickLabel',[]);end;
  set(gca,'XTickLabel',[]);
  text(length(idx1)/30,kmax-.5,['k jumps: ',num2str(N(i)),';  acceptance: ',num2str(N(i)/length(idx1),'%6.4f')],'FontSize',14);

  h1 = subplot('Position',[loc(1,i+NC) loc(2,i+NC) spw sph]);
  plot(B(idx2,4));
  set(gca,'XLim',[0 length(idx1)],'YLim',[kmin kmax],'FontSize',14);
  if(i == 1);ylabel('k');end;
  set(gca,'YTick',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]);
  if(i ~= 1);set(gca,'YTickLabel',[]);end;
  text(length(idx2)/30,kmax-.5,['k jumps: ',num2str(N(i+NC)),';  acceptance: ',num2str(N(i+NC)/length(idx2),'%6.4f')],'FontSize',14);
  xlabel('rjMCMC sweeps (no chain thinning)');


end;

saveas(fig1,'mixing.png','png');
return;
