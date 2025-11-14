function [] = save_exmpl()

file = 'bla.bin';
fid = fopen(file,'w+');

a = zeros(100,11);
j = 1;
for i = 1:1000
    
    a(j,:) = [i 1 2 3 4 5 6 7 8 9 10];
    
    j = j + 1;
    if(j > 100)
%        fprintf(fid,'%10.4f\n',a);
        fwrite(fid,a(:,1),a(:,2),'float');
        j = 1;
    end

end;

fclose(fid);
return;

