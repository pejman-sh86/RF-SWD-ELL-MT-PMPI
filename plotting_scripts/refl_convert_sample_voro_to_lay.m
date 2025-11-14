function [] = refl_convert_sample_voro_to_lay(sample,NPV,NVMX,NMISC,NFREQ)

%NPV = 6;
%NVMX = 20;

load(sample);
filebase = strrep(sample,'voro_sample.mat','');
sample2  = strcat(filebase,'sample.mat');
NSMP = size(A,1);
A2 = A;
A2(:,5:5+NVMX*NPV-1) = 0.;
k = A(:,4);
voro = zeros(NVMX,NPV,NSMP);
voroidx = zeros(size(voro));
for ismp = 1:NSMP;
for ivo = 1:k(ismp);
  ipar2 = (ivo-1)*NPV+1+4;
  voro(ivo,:,ismp) = A(ismp,ipar2:ipar2+NPV-1);
  voroidx(ivo,:,ismp) = 1;
end;
end;
voroidx(find(voro == -100)) = 0;

for ismp = 1:NSMP;
%for ismp = 1:1;
%  disp('voro:')
%  voro(1:k(ismp),:,ismp)
  [par,nunique] = refl_voro_to_lay(voro(:,:,ismp),voroidx(:,:,ismp),k(ismp),NPV,NVMX);
  A2(ismp,2) = k(ismp);
  A2(ismp,3) = sum(sum(voroidx(:,:,ismp)));
  A2(ismp,4) = nunique;
  A2(ismp,5:5+length(par)-1) = par;
  clear par;
end;
A2 = [A2,A(:,end-5:end)];
A=A2;
save(sample2,'A');

%save tmp.mat voro voroidx
return;
