function [] = rf_convert_sample_laynode_to_lay(sample,NPL,NVMX,NMISC,NRF,NMODE,NMODE_ELL,imt,iseis,ivref,ivpvs,hmax);
%NPL = 6;
%NVMX = 20;
%NMISC = 4;
load(sample);
%if(size(A,2)==136);
%  A=[A(:,1:127),A(:,128:129),A(:,128:end)];
%  A(:,128)=A(:,127);
%end;
filebase = strrep(sample,'voro_sample.mat','');
sample2  = strcat(filebase,'sample.mat');
velreffile  = strcat(filebase,'vel_ref.txt');

%if(ivref == 1);
%  tmp = dlmread(velreffile);
%  vel_ref = tmp(2:end,:);
%  vel_ref
%end;
if(ivref == 1);
  tmp = dlmread(velreffile);
  NREF = tmp(1,1);
  ntmp = tmp(1,2);
  NPREM = ntmp-NREF;
  vel_ref = tmp(2:NREF+1,:);
%  vel_ref(NREF+1,:) = [hmax,tmp(end,2:end)];
  vel_prem = tmp(NREF+1:NREF+1+NPREM,:);

NREF
NPREM
vel_ref
end
NSMP = size(A,1);
k = A(:,4);
voro = zeros(NVMX,NPL,NSMP);
voroidx = zeros(size(voro));
for ismp = 1:NSMP;
for ivo = 1:k(ismp);
  ipar2 = (ivo-1)*NPL+1+4;
  voro(ivo,:,ismp) = A(ismp,ipar2:ipar2+NPL-1);
  voroidx(ivo,:,ismp) = 1;
end;
end;
voroidx(find(voro < -99)) = 0;

% Test ref model
%ztmp = [0.:150./1000.:150.];
%for i =1:1000;
%  [vref(i)]=rf_getref(ztmp(i),vel_ref);
%end;
%save tmp.mat ztmp vref;
if iseis == 1 % my change
if(ivref == 1);
  for ismp = 1:NSMP;
  for ivo = 1:k(ismp);
    if(voroidx(ivo,2,ismp) == 1);
      %voro(ivo,1,ismp),vel_ref
      [vref,vpvsref]=rf_getref(voro(ivo,1,ismp),vel_ref);
      voro(ivo,2,ismp) = vref + voro(ivo,2,ismp);
    end
    if (ivpvs == 1) %my change
    if(voroidx(ivo,3,ismp) == 1);
      [vref,vpvsref]=rf_getref(voro(ivo,1,ismp),vel_ref);
      voro(ivo,3,ismp) = vpvsref + voro(ivo,3,ismp);
      %disp([voro(ivo,1,ismp),vref,voro(ivo,2,ismp)]);
    end
    end
  end;
  end;
end;
end

% if (imt==1); voro(:,end,:) = 10.0 .^ (-voro(:,end,:)); end;
if (imt==1); voro(:,end,:) = -voro(:,end,:); end;

for ismp = 1:NSMP;
    %if (ismp==1)%%%%%%%%%%%%%%%%%%%%%%
  [par,nunique(ismp)] = rf_laynode_to_lay(voro(:,:,ismp),voroidx(:,:,ismp),k(ismp),NPL,NVMX);
    %end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;
%par%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nmax = max(nunique);
nmax = max(nunique)+1; %%%%%%%%%%%%%%%%%%%Pejman
A2 = zeros(NSMP,nmax*NPL+4+6*NRF+2*NMODE+2*NMODE_ELL+NMISC+6+1); %I have added 1 for MT
A2(:,1) = A(:,1);

for ismp = 1:NSMP;
%for ismp = 1:1;
%  disp('voro:')
%  voro(1:k(ismp),:,ismp)
  [par,nunique(ismp)] = rf_laynode_to_lay(voro(:,:,ismp),voroidx(:,:,ismp),k(ismp),NPL,NVMX);
  A2(ismp,2) = k(ismp);
  A2(ismp,3) = sum(sum(voroidx(:,:,ismp)));
  A2(ismp,4) = nunique(ismp);
  A2(ismp,5:5+length(par)-1) = par;
  clear par;
end;

nmax
NVMX
NPL
save tmp.mat
A2(:,nmax*NPL+4+1:end) = A(:,NVMX*NPL+4+1:end);
A=A2;
save(sample2,'A');

%save tmp.mat voro voroidx
return;
