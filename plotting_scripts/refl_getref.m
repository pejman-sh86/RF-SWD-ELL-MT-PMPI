function [vref,rref]=refl_getref(z,vel_ref);
%%
%% Get linear interp
%%

NVELREF = size(vel_ref,1);

if(z >= vel_ref(NVELREF,1));
  vref = vel_ref(NVELREF,2);
  rref = vel_ref(NVELREF,3);
elseif(z == 0.);
  vref = vel_ref(1,2);
  rref = vel_ref(1,3);
else
  iint = 0;
  for ipar=1:NVELREF
    if((z-vel_ref(ipar,1)) <= 0.);
      break;
    end;
    iint = iint + 1;
  end;
  grad = (vel_ref(iint+1,2)-vel_ref(iint,2))/(vel_ref(iint+1,1)-vel_ref(iint,1));
  dz = (z-vel_ref(iint,1));
  vref = vel_ref(iint,2) + dz*grad;
  grad = (vel_ref(iint+1,3)-vel_ref(iint,3))/(vel_ref(iint+1,1)-vel_ref(iint,1));
  rref = vel_ref(iint,3) + dz*grad;
end;
return;
