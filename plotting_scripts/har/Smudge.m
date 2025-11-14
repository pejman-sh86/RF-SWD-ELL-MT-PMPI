%============================================================================
%                    SMUDGE.M
%============================================================================

function [RLcalc,angle,RLnew,newang]=Smudge(freqs,angs,RL,cs,cr,cb,ship,...
                                            Bft,interp_back,newang);


% Take RL (dB) from layer model and perform beam smearing on it
% assuming sheet dipole noise sources. Writing
% everything as a function of sin theta (at array) the beam smear is a
% true convolution which can be done by FTing.
% Beams are calculated separately and FTed in Beams.m
%
% C H Harrison Saclantcen May 2003
%
%	freqs	= array of freqs
%	angs	= array of angles (at the bottom)
%	RL		= RL(freqs,angs) is RL (dB) from a model
%	cs		= sound speed at the surface
%	cr		= sound speed at the array
%	cb		= sound speed at the bottom
%  ship  = ship spike (a fudge) ...  1 spoils RL contrast; 0 = no ship
%	Bft	= FT of beam shape from Beams.m (essentially a LUT)
%  interp_back = interpolate back to uniform angle spacing 'yes'/'no'
%  newang = array of desired angles for result

nf = length(freqs);
nb = length(angs);

% work in power
R2=10.^(-RL/10);
M=length(Bft(1,:))/2;	% Predetermined!!!
M2=M/2;	% vert to horiz
MM=2*M;	% expand to avoid wrap around
%
% work in linear sin(th)
% na is number of angles at array
%
na = M2+1;      
sinth = linspace(0,1,na);	% sin(theta) at array is linear 
th = asin(sinth);
thdeg = th*180/pi;

% Set up angles at surface and bottom using Snell. Note high cb results
% in refraction turning point with no loss. Therefore thb=0 works.
thb = acos(min(1,cb/cr*cos(th)));	%takes care of high or low cb

% Note high cs results in refraction turning point with no noise
% contribution. Therefore ths=0 works.

ths = acos(min(1,cs/cr*cos(th)));	%takes care of high or low cs
% interpolate R2 to angles at the array for all freqs
Rb = interp1(angs,R2',thb*180/pi)';
%
% For single frequency need to do this!
%
if nf < 2
  Rb = Rb';
end
Rs = .99; 	% fake surface loss
% absorption per ray cycle

a = freqs'*1e-7;	% atten in water (0.1dB/km@1kHz= 10^-7 dB/m/Hz)
A = 10.^(-200*a*(1./(sinth+eps))/10); % approx absorption


% N goes from 1 (upwards) to M2 (horiz) to M (down) followed by blank to MM

N = zeros(nf,MM);
ship = ship./(sin(ths)+.01);        % afudge that works
N(:,M2+1:-1:1) = (ones(nf,1)*(sin(ths)+ship))./(1-(Rb.*A)*Rs);% UP
N(:,M2+1:M+1) = N(:,M2+1:-1:1).*Rb;  % DOWN
Nft = fft(N,[],2);
AR = real(ifft(Bft.*Nft,[],2));

%*************************************************
% divide up by down

Rcalc=AR(:,M2+1:M)./AR(:,M2+1:-1:2);
%*************************************************

RLcalc = -10*log10(real(Rcalc));
% really 90 deg is at M2+1 so cheat

RLcalc(:,M2+1) = RLcalc(:,M2);
angle = asin((0:M2)/M2)*180/pi;
if strcmp(interp_back,'yes')   % interpolate back to uniform angle
   RLnew = interp1(angle,RLcalc',newang)';
   if nf < 2
     RLnew = RLnew';
   end
else
   RLnew  = 0;
   newang = 0;
end

return

