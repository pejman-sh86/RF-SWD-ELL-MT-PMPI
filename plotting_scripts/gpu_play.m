function []=gpu_play();
s1 = gpuArray(rand(4000));
s2 = gpuArray(rand(4000));
s3 = gpuArray(rand(4000));

t1 = tic;
[o1, o2] = arrayfun(@aGpuFunction, s1, s2, s3);

d = gather(o2);
t2 = toc(t1);
disp(['time: ',num2str(t2)]);

return;

function [o1, o2] = aGpuFunction(a, b, c)
o1 = a + b;
o2 = o1 .* c + 2;
return;
