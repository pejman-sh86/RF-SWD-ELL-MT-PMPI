for i = 1:169;
  name = ['p',num2str(i,'%03i'),'_0975_2700.txt'];
  emPlotAuvSummary(name,4,false,'track_environment_z.dat',false,false,255,4,i);
end;

