% function jmbench

disp('Matlab benchmark')
disp('Jan Mandel, December 2000')
disp('Last updated August 2004')

disp(['Matlab version ',version])

for iii=1:3
disp(' ')

f='%7.3f\n';
fprintf('1. large LU: ')
n=1000;m=3;
A=ones(n)+eye(n);
tic
for i=1:m,R=chol(A); end
t=toc; fprintf(f,t)

fprintf('2. small LU: ')
n=100;m=3000;
A=ones(n)+eye(n);
tic
for i=1:m,R=chol(A);end
t=toc; fprintf(f,t)

fprintf('3. sparse  : ')
n=2000;m=50;k=10;
B=ones(n,2*m+1);
B(:,m+1)=m+1;
d=[-m:m];
tic
A=spdiags(B,d,n,n);
R=chol(A);
R=A-R'*R;
t=toc; fprintf(f,t)

pause(1)

end
disp(' ')
disp('end of jmbench')

