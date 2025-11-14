function [] = plot_model_hol_3();
m1 = [1.8488498      1519.4803      1.2316569 0.081320648      1.3387829 1568.1905      1.9674834  3.0537816e-17 1540.2047      2.0743406 0.00024013637];
m2 = [1.8654916      1527.2777      1.3633663    0.099923128      1.9661114 1559.2479      2.0966138    0.054608573     0.88882545       1517.635 1.942428     0.14655505      1544.4705 2.0722928     0.12387149];

nl(1) = 3;
nl(2) = 4;
symb = ['- '; '--'];
c = cellstr(symb);

for i = 2:2
if i == 1
  m = m1;
else 
  m = m2;
end;
figure(1);
hold on, box on;

line = strcat(char(c(i)),'k');
% Interfaces:
if(nl(i) == 2)
  plot([1 2.24],[m(1) m(1)],line);
elseif(nl(i) == 3)
  plot([1 2.24],[m(1) m(1)],line);
  plot([1 2.24],[m(1)+m(5) m(1)+m(5)],line);
else
  plot([1 2.24],[m(1) m(1)],line);
  plot([1 2.24],[m(1)+m(5) m(1)+m(5)],line);
  plot([1 2.24],[m(1)+m(5)+m(9) m(1)+m(5)+m(9)],line);
end

% Velocity
line = strcat(char(c(i)),'b');
if(nl(i) == 2)
plot([m(2)/1000 m(2)/1000],[0 m(1)],line);
plot([m(5)/1000 m(5)/1000],[m(1) m(1)+0.4],line);
elseif(nl(i) == 3)
plot([m(2)/1000 m(2)/1000],[0 m(1)],line);
plot([m(6)/1000 m(6)/1000],[m(1) m(1)+m(5)],line);
plot([m(9)/1000 m(9)/1000],[m(1)+m(5) m(1)+m(5)+0.4],line);
else
plot([m(2)/1000 m(2)/1000],[0 m(1)],line);
plot([m(6)/1000 m(6)/1000],[m(1) m(1)+m(5)],line);
plot([m(10)/1000 m(10)/1000],[m(1)+m(5) m(1)+m(5)+m(9)],line);
plot([m(13)/1000 m(13)/1000],[m(1)+m(5)+m(9) m(1)+m(5)+m(9)+0.4],line);
end

% Density
line = strcat(char(c(i)),'r');
if(nl(i) == 2)
plot([m(3) m(3)],[0 m(1)],line);
plot([m(6) m(6)],[m(1) m(1)+0.4],line);
elseif(nl(i) == 3)
plot([m(3) m(3)],[0 m(1)],line);
plot([m(7) m(7)],[m(1) m(1)+m(5)],line);
plot([m(10) m(10)],[m(1)+m(5) m(1)+m(5)+0.4],line);
else
plot([m(3) m(3)],[0 m(1)],line);
plot([m(7) m(7)],[m(1) m(1)+m(5)],line);
plot([m(11) m(11)],[m(1)+m(5) m(1)+m(5)+m(9)],line);
plot([m(14) m(14)],[m(1)+m(5)+m(9) m(1)+m(5)+m(9)+0.4],line);
end

set(gca,'YDir','reverse');
%set(gca,'XTick',[1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2],...
%         'XTickLabel',{'1';'';'';'';'';'1.5';'';'';'';'';'2.0';'';''});
%set(gca,'YTick',[0 0.5 1 1.5 2],'YTickLabel',{'0';'';'1';'';'2'});
if(nl(i) == 2)
  set(gca,'XLim',[1 2.24],'YLim',[0 m(1)+0.4]);
elseif(nl(i) == 3)
  set(gca,'XLim',[1 2.24],'YLim',[0 m(1)+m(5)+0.4]);
else
  set(gca,'XLim',[1 2.24],'YLim',[0 m(1)+m(5)+m(9)+0.4]);
end;

end;

return;

figure(2);
for j = 1:11
  subplot(3,4,j)
  plot(mm(:,j))
end;
figure(3);
hold on;
for j = 11:23
  plot([mm(j,2)/1000 mm(j,2)/1000 mm(j,6)/1000 mm(j,6)/1000 ... 
        mm(j,9)/1000 mm(j,9)/1000],[0 mm(j,1) mm(j,1) ... 
        mm(j,1)+mm(j,5) mm(j,1)+mm(j,5) mm(j,1)+mm(j,5)+0.4]);
%  plot([mm(j,6)/1000 mm(j,6)/1000],[mm(j,1) mm(j,1)+mm(j,5)]);
%  plot([mm(j,9)/1000 mm(j,9)/1000],[mm(j,1)+mm(j,5) mm(j,1)+mm(j,5)+0.4]);
end;
set(gca,'YDir','reverse');

figure(4);
hold on;
for j = 11:23
  plot([mm(j,3) mm(j,3) mm(j,7) mm(j,7) mm(j,10) mm(j,10)],...
       [0 mm(j,1) mm(j,1) mm(j,1)+mm(j,5) mm(j,1)+mm(j,5) ... 
        mm(j,1)+mm(j,5)+0.4]);
end;
set(gca,'YDir','reverse');

return;

