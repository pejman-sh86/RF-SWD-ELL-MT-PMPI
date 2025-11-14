track=dlmread('track_environment.dat');
for i=1:170;k(i)=track(i,1);h(i,1:k(i))=track(i,[2:4:2+(k(i)-1)*4]);end;
for i=1:170;z(i,1:k(i))=cumsum(h(i,1:k(i)));end;
for i=1:170;track(i,[2:4:2+(k(i)-1)*4])=z(i,1:k(i));end;
save('track_environment.dat','track','-ascii');
