function []=cmt_noise_floor();

%% Read moment rate function:
dimensions = load('input_dimensions.dat');

NMRF  = dimensions(1);
NMETA = dimensions(2);
NSTN  = dimensions(3);
NDAT  = dimensions(4);

%% Read meat data:
obsdat = load('input_obsdat.dat');
dobs = obsdat;
%% Read NTSMP:
ntrace = load('input_ntrace.dat');
NTSMP = ntrace;
csum = [0;cumsum(ntrace)];
pred_py=load('py_pred.dat');
res=dobs(3,:)-pred_py';

for istn=1:NSTN;

  sd(istn) = std(res(csum(istn)+1:csum(istn+1)));

end;
save('input_sdmin.dat','sd','-ascii')
%% Read meat data:
%OPEN(UNIT=20,FILE='metafile.dat',FORM='formatted',STATUS='OLD',ACTION='READ')
%Do ii = 1,NSTN
%  READ(20,*) stn_name(ii),metadata_arr(ii,:)
%  tmpname = stn_name(ii)
%  ipos = SCAN(tmpname, 'ZNE12', back=.true.)
%  comp(ii) = tmpname(ipos:ipos)
%  %PRINT*,comp(ii),'   ',stn_name(ii)
%ENDDO

return;
