function [] = save_wl_txt_sample(sample,A);

%load(sample);

txtfile = strrep(sample,'.mat','.txt');

fid = fopen(txtfile,'w');
fmt1 = '%10i';
fmt2 = '%10i%15.9f';
N = size(A,1);
for i=1:N;
  fprintf(fid,fmt1,A(i,1));
  for j=2:2:2*A(i,1)+1;
    fprintf(fid,fmt2,A(i,j),A(i,j+1));
  end;
  fprintf(fid,'\n');
end;

return;
