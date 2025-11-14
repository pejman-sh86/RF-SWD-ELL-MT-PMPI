
 d(1).g(:,:) = load('CWG202.dat');
 d(2).g(:,:) = load('CWG203.dat');
 d(3).g(:,:) = load('CWG205.dat');
 d(4).g(:,:) = load('CWG214.dat');
 d(5).g(:,:) = load('CWG215.dat');
 d(6).g(:,:) = load('CWG219.dat');
 d(7).g(:,:) = load('CWG504.dat');
 d(8).g(:,:) = load('CWG602.dat');
 d(9).g(:,:) = load('CWG607.dat');
 d(10).g(:,:) = load('CWG613.dat');
 d(11).g(:,:) = load('CWG618.dat');
 d(12).g(:,:) = load('CWG625.dat');
 d(13).g(:,:) = load('CWG626.dat');
 d(14).g(:,:) = load('CWG901.dat');
 d(15).g(:,:) = load('DART21401_notide.dat');
 d(16).g(:,:) = load('DART21413_notide.dat');
 d(17).g(:,:) = load('DART21418_notide.dat');
 d(18).g(:,:) = load('DART21419_notide.dat');
 d(19).g(:,:) = load('GPS801.dat');
 d(20).g(:,:) = load('GPS802.dat');
 d(21).g(:,:) = load('GPS803.dat');
 d(22).g(:,:) = load('GPS804.dat');
 d(23).g(:,:) = load('GPS806.dat');
 d(24).g(:,:) = load('GPS807.dat');
 d(25).g(:,:) = load('GPS811.dat');
 d(26).g(:,:) = load('GPS812.dat');
 d(27).g(:,:) = load('GPS813.dat');
 d(28).g(:,:) = load('GPS815.dat');
 d(29).g(:,:) = load('GPSB801.dat');
 d(30).g(:,:) = load('GPSB802.dat');
 d(31).g(:,:) = load('GPSB803.dat');
 d(32).g(:,:) = load('GPSB804.dat');
 d(33).g(:,:) = load('GPSB806.dat');
 d(34).g(:,:) = load('GPSB807.dat');
 d(35).g(:,:) = load('GPSB811.dat');
 d(36).g(:,:) = load('GPSB813.dat');
 d(37).g(:,:) = load('OBPG1002.dat');
 d(38).g(:,:) = load('OBPG1006.dat');
 d(39).g(:,:) = load('OBPG2672.dat');
 d(40).g(:,:) = load('OBPG2673.dat');
 d(41).g(:,:) = load('OBPG5741.dat');
 d(42).g(:,:) = load('OBPG5742.dat');
 d(43).g(:,:) = load('OBPG5861.dat');
 d(44).g(:,:) = load('OBPG5862.dat');
 d(45).g(:,:) = load('OBPGN2672.dat');
 d(46).g(:,:) = load('OBPGN5742.dat');
 d(47).g(:,:) = load('OBPGN5861.dat');
 d(48).g(:,:) = load('OBPGN5862.dat');
 d(49).g(:,:) = load('TDG618.dat');
 d(50).g(:,:) = load('TDG625.dat');
 d(51).g(:,:) = load('TDG626.dat');

figure();
for isub=1:49;
  h(isub)=subaxis(7,7,isub,'Spacing',0.005,'Padding',0,...
                 'ML', 0.04,'MR',0.005,'MB',.06,'MT',.005);
  a = 50;
  b = 100;

  plot(d(isub).g(:,1),d(isub).g(:,2));
  %plot(d(isub).g(a:b,1),d(isub).g(a:b,2));
  
  %sd(isub) = std(d(isub).g(a:b,2));
end;

