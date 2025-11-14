function [] = save_wl_txt_sample(sample);

load sample;

txtfile = strrep(sample,'.mat','.txt');

fid = fopen(txtfile,'w');
fmt1 = '%18.8e\n';
N = size(A,1);
for i=1:N;
  for j=2:2:2*A(i,1)+1;
    fprintf(fid,fmt1,A(i,j),A(i,j+1));
  end;
end;

return;
