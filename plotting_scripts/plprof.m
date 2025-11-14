function [h] = plprof(m,zmax,npl,col,idx,offset);
col = char(col);
nl = (length(m)-(npl-1))/npl;
z = 0;z_old = 0;
v = 0;v_old = 0;
il2 = 1;
for il = 1:nl
  z = m(((il - 1) * npl) + 1);
  v = m(((il - 1) * npl) + 1 + idx);
  vplt(il2) = v+offset;
  zplt(il2) = z_old;
  vplt(il2+1) = v+offset;
  zplt(il2+1) = z;
  z_old = z;
  v_old = v;
  il2 = il2 + 2;
end
v = m(end - (npl-2) + (idx-1));
vplt(il2) = v+offset;
zplt(il2) = z_old;
vplt(il2+1) = v+offset;
zplt(il2+1) = zmax;
h = plot(vplt,zplt,col,'Linewidth',1);
return;
