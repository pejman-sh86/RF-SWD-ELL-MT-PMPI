
load V0Th0.mat;

data(1,:) = sum(V0(7:9,:))/3;
data(2,:) = sum(V0(9:11,:))/3;
data(3,:) = sum(V0(11:13,:))/3;
data(4,:) = sum(V0(13:16,:))/4;
data(5,:) = sum(V0(16:20,:))/5;
data(6,:) = sum(V0(20:25,:))/6;
data(7,:) = sum(V0(25:31,:))/7;
data(8,:) = sum(V0(32:39,:))/8;

fr = freq([8 10 12 14 18 23 28 35]);
ang = Theta0;

%save bla data fr ang;
