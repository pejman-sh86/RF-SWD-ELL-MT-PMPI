function [] = write_oasr_input();

file = 'input.dat';
load geo_sin;


fid = fopen(file,'w');
fprintf(fid,'Holland Test\n N T\n');
fprintf(fid,'102\n');
fprintf(fid,'0.0         1511.00    0   0.00000000  0  1.029000  0\n');
z = 150;
for iline = 2:length(geo_sin)-1;
  
    
    v = geo_sin(iline,2);
%    a = v/1000 * geo_sin(iline,3);
    a = 1.5 * geo_sin(iline,3);
    r = geo_sin(iline,4);

    fprintf(fid,'%10.4f %10.4f 0 %10.6f 0 %10.6f 0 \n',z,v,a,r);
%    fprintf(1,'%10.4f %10.4f 0 %10.6f 0 %10.6f 0 \n',z,v,a,r)

    z = z + geo_sin(iline,1);
end

%z = z + geo_sin(end-1,1);
v = geo_sin(end,2);
a = geo_sin(end,3);
r = geo_sin(end,4);
fprintf(fid,'%10.4f %10.4f 0 %10.6f 0 %10.6f 0 \n',z,v,a,r);
fprintf(fid,'100 4500 241 0\n');
fprintf(fid,'0 90 361 0\n');
fprintf(fid,'90 0 10 10\n');
fprintf(fid,'0 15 12 5\n');
fprintf(fid,'0 90 12 15\n');
fprintf(fid,'100 2500 1 1\n');

fclose(fid);
return;
