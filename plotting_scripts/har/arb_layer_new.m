%===========================================================================
   function [t1,V] = arb_layer_new(c1,h1,r1,a1,f,ang);
%===========================================================================

%Set up array sizes
nang=length(ang);
nfreq=length(f);
%Generate angle variable
t1 = fliplr(ang); %In degrees for output
thetar=t1*pi/180;       %In radians for calculation
%Change the orientation of the layer variables and get the number of layers
c=c1';
r=r1';
h=h1';
alfa=a1';
last=length(c);
%Include the attenuation in a complex sound speed
const=2*pi*20*log10(exp(1));
c=c.*(1-i*alfa./sqrt(const*const-alfa));
%Set up variables to be used in explicit matrix multiplication
ma = ones(nang,nfreq);
mb = zeros(nang,nfreq);
mc = zeros(nang,nfreq);
md = ones(nang,nfreq);
%Set up a and gamma variables for each layer and angle
a = c*cos(thetar)/c(1);	
gamma = sqrt(1-a.*a)./(c*ones(1,nang));
%Set up gamma2 variable for each layer, angle and frequency
gamma2 = zeros(last,nang,nfreq);
for ii = 1:last
    gamma2(ii,:,:) = 2.*pi.* (squeeze(gamma(ii,:)).') * f;
end
%Loop over layers from the one before the halfspace to the water, going
%upwards towards the water
for m = last-1:-1:1
    %Set up "array elements"
    b0 = 2*r(m+1)     .*( gamma(m+1,:).') *ones(size(f));
    b1 =   r(m)       .*( gamma(m+1,:).') *ones(size(f));
    b2 =   r(m+1)     .*( gamma(m,:)  .') *ones(size(f));
    e0 = exp(i*h(m)   .*(squeeze(gamma2(m,:,:)))       );
    if nfreq < 2
      e0 = e0';
    end 
    me = (b1+b2)./e0./b0;
    mf = (b1-b2).*e0./b0;
    mg = (b1-b2)./e0./b0;
    mh = (b1+b2).*e0./b0;
    ma1 = ma.*me + mb.*mg;, mb1 = ma.*mf + mb.*mh;
    mc1 = mc.*me + md.*mg;, md1 = mc.*mf + md.*mh;
    ma = ma1;
    mb = mb1;
    mc = mc1;
    md = md1;
end
V = (-mb./ma).';
return


