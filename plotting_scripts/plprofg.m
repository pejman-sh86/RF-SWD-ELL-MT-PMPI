function [] = plprofg(m,hmax,col,idx);
npl = 4;
col = char(col);
nl = (length(m)-4)/npl;
plot([m(3) m(4)],[0 m(1)],col,'Linewidth',2);
plot([m(4) m(8)],[m(1) m(1)],col,'Linewidth',2);
h = m(1);h_old = m(1);
v = 0;v_old = m(4);
for il = 2:nl
    h = h + m(((il - 1) * npl) + 2);
    v = m(((il - 1) * npl) + 2 + idx);
    if(il > 1)
        plot([v_old v],[h_old h_old],col,'Linewidth',2);
    end
    plot([v v],[h_old h],col,'Linewidth',2);
    h_old = h;
    v_old = v;
  end
v = m(end - 2 + (idx-1));
plot([v v],[h_old hmax],col,'Linewidth',2);
plot([v_old v],[h_old h_old],col,'Linewidth',2);
return;
