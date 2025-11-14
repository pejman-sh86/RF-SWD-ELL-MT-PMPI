load tohoku_data_wl_sl1_sample.mat

NSMP = length(A);
k= A(:,1);
NK = 512;

A(:,1)=[];
B=A(:,1:2:end);
C=A(:,2:2:end);
NB = size(B,2);

for i=1:NB;
    for j=1:NK;
        idx=find(B(:,i)==j);
        if(idx);
            cf(j).m = [cf(j).m,C(idx,i)];
        end;
end;end;
