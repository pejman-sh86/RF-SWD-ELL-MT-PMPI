function [scsm,dscsm,sdang,dC] = acov_non_b(res,ang,widths,widthe,ifreq)

  res = res - mean(res);
  nang = length(ang);
  Q = (nang*(nang+1)/2)-nang+1;

% Variance:
  c(1) = sum(res.^2)/nang;
  dang(1) = 0;
% Off diagonal terms
  l = 2;
  for i = 1:nang-1
    for j = 1:nang-i
      c(l) = res(i) * res(i+j);
      dang(l) = ang(i+j) - ang(i);
      l = l + 1;
    end
  end
%
% Sorting the lags:
%
  [sdang,IDX] = sort(dang);
  sc = c(IDX);

  scsm(1) = sc(1);
  sd(1) = std(res);
  width = widths;
  for l = 2:Q
    th0 = sdang(l);
    thmn = th0-width/2;
    thmx = th0+width/2;
% Loop comps the weights in area around the midpoint
    if(thmn < 0)
      kmn = 2;
      k = 1;
      while(sdang(l+k)<thmx);kmx = l+k; k = k+1;end;
    elseif(thmx > sdang(end))
      k = 1;
      while(sdang(l-k)>thmn);kmn = l-k; k = k+1;end;
      kmx = Q;
    else
      k = 1;
      while(sdang(l-k)>thmn);kmn = l-k; k = k+1;end;
      k = 1;
      while(sdang(l+k)<thmx);kmx = l+k; k = k+1;end;
    end
    newwidth = sdang(kmx)-sdang(kmn);
    th = sdang(kmn:kmx);
    w = weight(th,th0,newwidth);
    w = w./sum(w);
    n(l) = length(w);
    scsm(l) = sum(sc(kmn:kmx).*w);
    sd(l) = std(sc(kmn:kmx))/sqrt(kmx - kmn);
    width = width+((widthe-widths)/Q);
    clear w;
  end % end l loop
%
% Set more damping for problem childs:
%
  if(ifreq == 1);    fact = 1;elseif(ifreq == 2);fact = 1;
  elseif(ifreq == 3);fact = 1;elseif(ifreq == 4);fact = 1;
  elseif(ifreq == 5);fact = 1;elseif(ifreq == 6);fact = 1;
  elseif(ifreq == 7);fact = 1.45 ;elseif(ifreq == 8);fact = 1.65 ;
  else;fact = 1;end;
  x = exp(-fact*[1:Q]/(Q));
  dscsm = scsm .* x;
%
% Resort and fill in matrix
%
  [sIDX,IDX2] = sort(IDX);
  rdscsm = dscsm(IDX2);
  rscsm = scsm(IDX2);
  l = 2;
  for i = 1:nang
    dC(i,i) = rdscsm(1);
    C(i,i) = rscsm(1);
    for j = i+1:nang
      dC(i,j) = rdscsm(l);
      dC(j,i) = dC(i,j);
      C(i,j) = rscsm(l);
      C(j,i) = C(i,j);
      l = l+1;
    end
  end
return;

function [w] = weight(th,th0,newwidth)

w = cos((th-th0)*pi/180*90/newwidth);

return;
