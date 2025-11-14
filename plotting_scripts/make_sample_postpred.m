
load HON_voro_sample.mat
AS = A;
filebase = 'HON';
irandom = 0;
step = 10;
NSMP = 5000;

%% 


sample_postpred_file = strcat(filebase,'_sample_postpred.dat');

if irandom
    r = randi([1, size(AS,1)], 1, NSMP);
    AS = AS(r,:);
else
    AS = AS(1:step:end,:);
end 

save(sample_postpred_file, 'AS', '-ascii')