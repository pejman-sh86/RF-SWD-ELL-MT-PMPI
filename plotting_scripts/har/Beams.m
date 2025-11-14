%===========================================================================
%                   BEAMS.M
%===========================================================================

function [Bft] = Beams(z,freqs,cr,M);

%----------------------------------------------------------------------------
% Calculate Beams and FT
%
% C H Harrison Saclantcen May 2003
%
%  z	  = array of hp separations in m eg [.5:.5:16]
%  Bft	  = FT of beam pattern at each freq
%  freqs  = array of freqs
%  cr	  = sound speed at the array
%  M      = number of points in final range -pi/2 to pi/2
%----------------------------------------------------------------------------


M2 = M/2;	% vert to horiz
MM = 2*M;
nf = length(freqs);
% centre beams on M2 but place in MM array
z = z';
sinth = linspace(-2,2,MM+1); % crazy but necessary for convolved grating lobes
sinth = sinth(1:MM);
k = 2*pi*freqs/cr;
% shading
h = hamming(length(z))*ones(1,MM);

for f = 1:nf
  B(f,:)=ifftshift(abs(sum(h.*exp(-i*k(f)*(z*sinth)),1)).^2);
%  disp(num2str(f))
end

Bft = fft(B,[],2);
Bft = Bft./(Bft(:,1)*ones(1,MM));

return


