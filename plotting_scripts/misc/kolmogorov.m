%
% Script which is based on Kolmogorov - Smirnov method
% 
% From the given dataset the mean and standard deviation
% is calculated. For this mean and standard deviation a
% 'theoretical normal distribution' is calculated. After that
% the calculated en the given dataset are compared using the
% Kolmogorov - Smirnov method.
% 
% 
% usage: kolmogorov(data)
% 
% example usage: kolmogorov(data)
%
% Remark: 
%          - Data is supposed to be a columnvector thus:
%                    data = | value 1 |
%                           | value 2 |
%                           |    ..   |
%                           | value n |
%
%
% The writer of this M-file doesn't take any responsibility which 
% occur due to usage of this script in any way.
%
% A. Saglam
% asaglam@worldonline.nl
%
% 21/09/00
% 
%
function[pass,value,bereik2,opp,dplot]=kolmogorov(data,k);
%
%
%
[a b]=size(data);
nwk=round(3000/a);
%
% Override the mean to actually test for G(0,sd)!
%
%gem=mean(data);
gem = 0;
s=std(data);
temp=zeros(a,1);
for i=1:1:a
    temp(i,1)= (data(i,1)-gem)/s;
end
temp2=sort(temp);
data=temp2;
%gem=mean(data);
%s=std(data);
%
%
% The theoretical gaussian curve is
% divided in three parts.
%
% Part 1: area between -10 and min(data)
% Part 2: area between min(data) and max(data)
% Part 3: area between max(data) and +10
%
%
% calculate area of part 1
%
temp2=zeros(1000,1);
minimaal=-10.0;
maximaal=min(data);
bereik=abs(maximaal-minimaal);
xstap=abs(bereik/(1000));
p=minimaal;
for i=1:1:1000
   temp2(i,1)= (1/(s*sqrt(2*pi))) * exp(-0.5*((p-gem)/s)^2);
   p=p+xstap;
end

x=[minimaal+xstap:xstap:maximaal]';

nv1=trapz(x(:,1),temp2(:,1));
%
%
% calculate area of part 2
% calculate stepwise area of part 2
%
temp2=zeros(nwk*a,1);
minimaal=min(data);
maximaal=max(data);
bereik=abs(maximaal-minimaal);
xstap=bereik/(nwk*a);
p=minimaal;
for i=1:1:nwk*a
      temp2(i,1)= (1/(s*sqrt(2*pi))) * exp(-0.5*((p-gem)/s)^2);
      p=p+xstap;
end
x=[minimaal:xstap:maximaal]';
n=2;
for i=2:1:a*nwk
         ropp(i,1)=trapz(x((1:n),1),temp2((1:n),1));
         ropp(i,1)=ropp(i,1)+nv1;
         n=n+1;
end;
nv2=ropp(nwk*a,1);
oppverloop=ropp;
%
%
% calculate area of part 3
% 
temp2=zeros(1000,1);
minimaal=max(data);
maximaal=10;
bereik=abs(maximaal-minimaal);
xstap=abs(bereik/(1000));
p=minimaal;
for i=1:1:1000
      temp2(i,1)= (1/(s*sqrt(2*pi))) * exp(-0.5*((p-gem)/s)^2);
      p=p+xstap;
end
x=[minimaal+xstap:xstap:maximaal]';
nv3=trapz(x,temp2)+nv2;
%
%
% Correct area of part 2 by use of total
% area of part 1, 2 and 3 
%
maxopp=nv3;
for i=1:1:a*nwk
   opp(i,1)=oppverloop(i,1)/maxopp;
end
%
%
% find values which are the same for
% the normal distribution and the 
% distribution to be tested
%
minimaal=min(data);
maximaal=max(data);
bereik=abs(maximaal-minimaal);
xstap=bereik/(nwk*a);
bereik2=minimaal:xstap:maximaal;
for i=1:1:a
    uverschil=abs(bereik2-data(i,1));
    minu=min(uverschil);
    [c d]=find(uverschil==minu);
    tpos(i,1)=d;
end
%
%
% calulate differences between 
% normal distribution and distribution
% to be tested
%
verschil=zeros(a,1);
maxverschil=0;
for i=1:1:a
   ttpos=tpos(i,1);
   if ttpos >= (nwk*a)
      ttpos=nwk*a; 
   end
      verschil(i,1)= abs( opp(ttpos,1)-(i/a));
      if verschil(i,1)>=maxverschil
         rood=(i/a);
         blauw=opp(ttpos,1);
         xas=i;
         maxverschil=verschil(i,1);
      end
end
%
%
% prepare data for being plot
%
for i=1:1:a
   d2=i-1;
   if i==1
      d2=1;
   end
   d=tpos(d2,1);
   e=tpos(i,1);
   dplot(d:e,1)=((i-1)/a);
   if i==xas
      xas1=e;
      xas2=e+1;
   end
   d=tpos(d2,1);
   e=tpos(i,1);
end
[aa bb]=size(dplot);
dplot=dplot(1:aa-1,1);
minimaal=min(data);
maximaal=max(data);
bereik=abs(maximaal-minimaal);
xstap=bereik/(nwk*a);
bereik2=[minimaal+xstap:xstap:maximaal]';

%
%
% Calculated Kolmogorov-Smirnov value
%
ksverschil=max(verschil);
value = ksverschil * sqrt(a);
%
%
% Theoretical value of Kolmogorov-Smirnov
% for a=0.05
%
tksverschil=0.90/sqrt(a);
%tksverschil=1.36/sqrt(a);
%
%
% Output
%
d1=sprintf('Calculated value : %0.5g',ksverschil);
%disp(d1);
d2=sprintf('Theoretical value : %0.5g',tksverschil);
%disp(d2);
%disp('warning: Theoretical value is only valid if the number of measurementsis greater than 30!');
%disp('So, is the experimental data normally distributed ? '); 
if ksverschil<= tksverschil 
  pass = 1;
%  disp('yes');
end; 
if ksverschil > tksverschil 
  pass = 0;
%  disp('no');
end; 
%
% 
% end of script
