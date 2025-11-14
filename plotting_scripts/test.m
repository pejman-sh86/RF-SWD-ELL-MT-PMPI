%some data
[x,y] = meshgrid(0:0.2:2,0:0.2:2);
u = cos(x).*y;
v = sin(x).*y;

%quiver plots
figure('Position',[10 10 1000 600],'Color','w');
hax_1 = subplot(1,2,1);

%left version (regular)
hq1 = quiver(x,y,u,v);

%get the line position (first handle)
hkid = get(hq1,'children');
X = get(hkid(1),'XData');
Y = get(hkid(1),'YData');
box on;

%right version (with annotation)
hax_2 = subplot(1,2,2);

jj = 1
for ii = 1:3:length(X)-1
  len(jj) = sqrt((Y(ii+1)-Y(ii))^2+(X(ii+1)-X(ii))^2); %get the angle
  jj = jj + 1;
end;
len = round(len/max(len)*127)+1;
cmap = jet(128); %colormap, 116 because angles goes up to 115 degrees

jj = 1;
for ii = 1:3:length(X)-1
  headwidth = 200 * sqrt((X(ii+1)-X(ii)).^2 + (Y(ii+1)-Y(ii)).^2); % set the headWidth, function of length of arrow
  ah = annotation('arrow',...
      'Color', cmap(len(jj),:),...
      'headStyle','cback1','HeadLength',2*headwidth,'HeadWidth',headwidth);
  set(ah,'parent',gca);
  set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
  jj = jj + 1;
end
box on;

%linkaxes([hax_1 hax_2],'xy');
