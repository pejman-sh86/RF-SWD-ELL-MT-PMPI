function [width,height,llh,llw] = ...
         calc_subplot(lm,rm,mm,tm,bm,nr_r,nr_c);

width = (1 - lm - rm - (nr_r-1)*mm)/nr_r;
width = floor(width*1000)/1000;
height = (1 - tm - bm - (nr_c-1)*mm)/nr_c;
height = floor(height*1000)/1000;

posl(1) = lm; 
posh(1) = bm; 
for i = 2:nr_r

    posl(i) = posl(i-1) + width + mm; 
    posl(i) = floor(posl(i)*1000)/1000;

end

for i = 2:nr_c

    posh(i) = posh(i-1) + height + mm; 
    posh(i) = floor(posh(i)*1000)/1000;

end

llw = [repmat(posl,1,nr_c)];
llh = [fliplr(posh)];

return;
