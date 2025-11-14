function [] = refl_convert_sample_laynode_to_lay(sample,NPV,NVMX,NMISC,NFREQ)

%NPV = 6;
%NVMX = 20;
%NFREQ = 11;
%NMISC = 4;

load(sample);
filebase = strrep(sample,'voro_sample.mat','');
sample2  = strcat(filebase,'sample.mat');
NSMP = size(A,1);
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
  [par,nunique(ismp)] = refl_laynode_to_lay(voro(:,:,ismp),voroidx(:,:,ismp),k(ismp),NPV,NVMX);
end;
save tmp.mat nunique
nmax = max(nunique)
A2 = zeros(NSMP,nmax*NPV+4+NFREQ+NFREQ+NMISC+6);
A2(:,1) = A(:,1);

for ismp = 1:NSMP;
%for ismp = 1:1;
%  disp('voro:')
%  voro(1:k(ismp),:,ismp)
  [par,nunique(ismp)] = refl_laynode_to_lay(voro(:,:,ismp),voroidx(:,:,ismp),k(ismp),NPV,NVMX);
  A2(ismp,2) = k(ismp);
  A2(ismp,3) = sum(sum(voroidx(:,:,ismp)));
  A2(ismp,4) = nunique(ismp);
  A2(ismp,5:5+length(par)-1) = par;
  clear par;
end;

nmax*NPV+4+1
size(A2,2)
length([nmax*NPV+4+1:size(A2,2)])
NVMX*NPV+4+1
size(A,2)
length([NVMX*NPV+4+1:size(A,2)])

A2(:,nmax*NPV+4+1:end) = A(:,NVMX*NPV+4+1:end);
A=A2;
save(sample2,'A');

%save tmp.mat voro voroidx
return;
