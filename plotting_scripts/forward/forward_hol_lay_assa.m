%================================================================
%   FORWARD.M
%================================================================
function [E,Fm] = forward_hol_lay_assa(m,F);
global icorr;
%
% Setting up the environment in the (nlay+2)x4 Array geo_sin:
%

nfreq = length(F(1).freq);
if(F(1).nmod == 7)
  % 2-Layer model
  h  = [ NaN    m(1) NaN ]';
  c  = [F(1).cw m(2) m(5)]';
  r  = [F(1).rw m(3) m(6)]';
  a  = [   0    m(4) m(7)]';
elseif(F(1).nmod == 11)
  %  3-Layer model
  h  = [ NaN    m(1) m(5) NaN  ]';
  c  = [F(1).cw m(2) m(6) m(9) ]';
  r  = [F(1).rw m(3) m(7) m(10)]';
  a  = [   0    m(4) m(8) m(11)]';
elseif(F(1).nmod == 15)
  %  4-Layer model
  h  = [ NaN    m(1) m(5) m(9)  NaN  ]';
  c  = [F(1).cw m(2) m(6) m(10) m(13)]';
  r  = [F(1).rw m(3) m(7) m(11) m(14)]';
  a  = [   0    m(4) m(8) m(12) m(15)]';
elseif(F(1).nmod == 19)
  %  5-Layer model
  h  = [ NaN    m(1) m(5) m(9)  m(13) NaN  ]';
  c  = [F(1).cw m(2) m(6) m(10) m(14) m(17)]';
  r  = [F(1).rw m(3) m(7) m(11) m(15) m(18)]';
  a  = [   0    m(4) m(8) m(12) m(16) m(19)]';
elseif(F(1).nmod == 23)
  %  6-Layer model
  h  = [ NaN    m(1) m(5) m(9)  m(13) m(17) NaN  ]';
  c  = [F(1).cw m(2) m(6) m(10) m(14) m(18) m(21)]';
  r  = [F(1).rw m(3) m(7) m(11) m(15) m(19) m(22)]';
  a  = [   0    m(4) m(8) m(12) m(16) m(20) m(23)]';
elseif(F(1).nmod == 27)
  %  7-Layer model
  h  = [ NaN    m(1) m(5) m(9)  m(13) m(17) m(21) NaN  ]';
  c  = [F(1).cw m(2) m(6) m(10) m(14) m(18) m(22) m(25)]';
  r  = [F(1).rw m(3) m(7) m(11) m(15) m(19) m(23) m(26)]';
  a  = [   0    m(4) m(8) m(12) m(16) m(20) m(24) m(27)]';
elseif(F(1).nmod == 31)
  %  7-Layer model
  h  = [ NaN    m(1) m(5) m(9)  m(13) m(17) m(21) m(25) NaN  ]';
  c  = [F(1).cw m(2) m(6) m(10) m(14) m(18) m(22) m(26) m(29)]';
  r  = [F(1).rw m(3) m(7) m(11) m(15) m(19) m(23) m(27) m(30)]';
  a  = [   0    m(4) m(8) m(12) m(16) m(20) m(24) m(28) m(31)]';
end

geo_sin=[h      c         a./(c/1000)      r];
%
% Compute reflectivity:
%
if((m > F(1).minlim) & (m < F(1).maxlim))
  E = 0;
  for i=1:nfreq
    [ref] = ref_nlay_fav3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
    ref         = -20*log10(abs(ref));
    Fm(i).dat = reshape(ref,length(ref),1);
    F(i).dat = reshape(F(i).dat,length(F(i).dat),1);
    if icorr == 1
        E = E + (Fm(i).dat-F(i).dat)'*F(i).c_inv*(Fm(i).dat-F(i).dat)/2;
    else
        E = E + length(F(i).dat)/2*log(sum((Fm(i).dat-F(i).dat).^2));
    end
  end
%  E = E/10^28;
else
  disp('Out of bounds')
  for i=1:nfreq
    [ref] = ref_nlay3(F(i).ang,geo_sin,F(1).freq(i));% compute Reflection
    ref         = -20*log10(abs(ref));
    Fm(i).dat = ref';
  end
  E = 10^50.;
end

return
% ---------------------------------------------------------------
% ...this is the end my fiend.
% EOF
