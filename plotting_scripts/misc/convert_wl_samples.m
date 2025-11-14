function [] = convert_wl_samples(filebase,N,M,NTW,IVRUP);

sample = strcat(filebase,'_sample.txt');
samplemat = strrep(sample,'.txt','.mat');

sample_sl1 = strcat(filebase,'_wl_sl1_sample.txt');
samplemat_sl1 = strrep(sample_sl1,'.txt','.mat');
sample_sl2 = strcat(filebase,'_wl_sl2_sample.txt');
samplemat_sl2 = strrep(sample_sl2,'.txt','.mat');
sample_sl3 = strcat(filebase,'_wl_sl3_sample.txt');
samplemat_sl3 = strrep(sample_sl3,'.txt','.mat');

sample_sr = strcat(filebase,'_wl_sr_sample.txt');
samplemat_sr = strrep(sample_sr,'.txt','.mat');

convert_sample(sample);
convert_wl_sample(sample_sl1);
if(NTW > 1);convert_wl_sample(sample_sl2);end;
if(NTW > 2);convert_wl_sample(sample_sl3);end;
if(IVRUP == 1);convert_wl_sample(sample_sr);end;

load(samplemat);
A(1:N,:)=[];
figure();plot(A(:,1));
A=A(1:M:end,:);
save(samplemat,'A');
save(sample,'A','-ascii');

load(samplemat_sl1);
A(1:N,:)=[];
A=A(1:M:end,:);
save(samplemat_sl1,'A');
save_wl_txt_sample(samplemat_sl1,A);
if(NTW > 1);
  load(samplemat_sl2);
  A(1:N,:)=[];
  A=A(1:M:end,:);
  save(samplemat_sl2,'A');
  save_wl_txt_sample(samplemat_sl2,A);
end;
if(NTW > 2);
  load(samplemat_sl3);
  A(1:N,:)=[];
  A=A(1:M:end,:);
  save(samplemat_sl3,'A');
  save_wl_txt_sample(samplemat_sl3,A);
end;

if(IVRUP == 1);
  load(samplemat_sr);
  A(1:N,:)=[];
  A=A(1:M:end,:);
  save(samplemat_sr,'A');
  save_wl_txt_sample(samplemat_sr,A);
end;
return;