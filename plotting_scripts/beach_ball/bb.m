function [xx,yy,X,Y]=bb(fm, centerX, centerY, diam, ta, color,ifill,alph)
%function bb(fm, centerX, centerY, diam, ta, color)
%
%	Draws beachball diagram of earthquake double-couple focal mechanism(s). S1, D1, and
%		R1, the strike, dip and rake of one of the focal planes, can be vectors of
%		multiple focal mechanisms.
%	fm - focal mechanism that is either number of mechnisms (NM) by 3 (strike, dip, and rake)
%		or NM x 6 (mxx, myy, mzz, mxy, mxz, myz - the six independent components of
%		the moment tensor). The strike is of the first plane, clockwise relative to north.
%		The dip is of the first plane, defined clockwise and perpedicular to strike, 
%		relative to horizontal such that 0 is horizontal and 90 is vertical. The rake is 
%		of the first focal plane solution. 90 moves the hanging wall up-dip (thrust),
%		0 moves it in the strike direction (left-lateral), -90 moves it down-dip
%		(normal), and 180 moves it opposite to strike (right-lateral).
%	centerX - place beachball(s) at position centerX
%	centerY - place beachball(s) at position centerY
%	diam - draw with this diameter.  If diam is zero, beachball
%		is drawn on stereonet.
%	ta - type of axis. If 0, this is a normal axis, if 1 it is a map axis. In case
%		of the latter, centerX and centerY are Lons and Lats, respectively.
%	color - color to use for quadrants of tension; can be a string, e.g. 'r'
%		'b' or three component color vector, [R G B].
%


[ne,n] = size(fm);
if n == 6
	for j = 1:ne
		[s1(j),d1(j),r1(j)] = mij2sdr(fm(j,1),fm(j,2),fm(j,3),fm(j,4),fm(j,5),fm(j,6));
	end
else
	s1 = fm(:,1);
	d1 = fm(:,2);
	r1 = fm(:,3);
end

r2d = 180/pi;
d2r = pi/180;
ampy = cos(mean(centerY)*d2r);

if ne > 1
	[ds,i] = sort(diam,1,'descend');
	diam = diam(i);
	s1 = s1(i);
	d1 = d1(i);
    r1 = r1(i);
 	centerX = centerX(i);
	centerY = centerY(i);
end

mech = zeros(ne,1);
j = find(r1 > 180);
r1(j) = r1(j) - 180;
mech(j) = 1;
j = find(r1 < 0);
r1(j) = r1(j) + 180;
mech(j) = 1;

% Get azimuth and dip of second plane
[s2,d2,r2] = AuxPlane(s1,d1,r1);

if diam(1) > 0
	hold on
end
for ev = 1:ne
	S1 = s1(ev);
	D1 = d1(ev);
	S2 = s2(ev);
	D2 = d2(ev);
	P = r1(ev);
	CX = centerX(ev);
	CY = centerY(ev);
	D = diam(ev);
	M = mech(ev);

if M > 0
   P = 2;
else
   P = 1;
end

if D1 >= 90
   D1 = 89.9999;
end
if D2 >= 90
   D2 = 89.9999;
end

phi = 0:.01:pi;
d = 90 - D1;
m = 90;
l1 = sqrt(d^2./(sin(phi).^2 + cos(phi).^2 * d^2/m^2));

d = 90 - D2;
m = 90;
l2 = sqrt(d^2./(sin(phi).^2 + cos(phi).^2 * d^2/m^2));

if D == 0
   stereo(phi+S1*d2r,l1,'k')
   hold on
   stereo(phi+S2*d2r,l2,'k')
end

inc = 1;
[X1,Y1] = pol2cart(phi+S1*d2r,l1);
if P == 1
   lo = S1 - 180;
   hi = S2;
   if lo > hi
      inc = -inc;
   end
   th1 = S1-180:inc:S2;
   [Xs1,Ys1] = pol2cart(th1*d2r,90*ones(1,length(th1)));
   [X2,Y2] = pol2cart(phi+S2*d2r,l2);
   th2 = S2+180:-inc:S1;
else
   hi = S1 - 180;
   lo = S2 - 180;
   if lo > hi
      inc = -inc;
   end
   th1 = hi:-inc:lo;
   [Xs1,Ys1] = pol2cart(th1*d2r,90*ones(1,length(th1)));
   [X2,Y2] = pol2cart(phi+S2*d2r,l2);
   X2 = fliplr(X2);
   Y2 = fliplr(Y2);
   th2 = S2:inc:S1;
end
[Xs2,Ys2] = pol2cart(th2*d2r,90*ones(1,length(th2)));

X = cat(2,X1,Xs1,X2,Xs2);
Y = cat(2,Y1,Ys1,Y2,Ys2);

if D > 0
   X = ampy*X * D/90 + CY;
   Y = Y * D/90 + CX;
   phid = 0:.01:2*pi;
   [x,y] = pol2cart(phid,90);
   xx = x*D/90 + CX;
   yy = ampy*y*D/90 + CY;
   if ta == 0
      if(ifill == 1);
        %fill(xx,yy,'w')
        fill(Y,X,color,'FaceAlpha',alph,'EdgeColor',color,'EdgeAlpha',alph)
      end;
      %line(xx,yy,'color','k','linewidth',0.5,);
      %line(X,Y,'color','k','linewidth',0.5,);
   else
      if(ifill == 1);
        fillm(yy,xx,'w')
        fillm(X,Y,color)
      end;
      linem(yy,xx,'color','k','linewidth',0.5);
   end
else
   if ta == 0
      if(ifill == 1);
        fill(X,Y,color)
      end
      line(xx,yy,'color','k','linewidth',1);
      line(X,Y,'color','k','linewidth',1);
   else
      if(ifill == 1);
        fillm(Y,X,color)
      end
      linem(yy,xx,'color','.k','linewidth',1);
      linem(Y,X,'color','.k','linewidth',1);
   end
   view(90,-90)
end

end
