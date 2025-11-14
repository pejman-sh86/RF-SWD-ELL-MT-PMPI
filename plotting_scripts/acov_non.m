function [scnew,dscnew,sdang,dC,sc,rscnew,rdscnew,sd,x,C,fact] = ...
         acov_non(res,ang,widths,widthe,ifreq);

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

  scnew(1) = sc(1);
  scnew_orig = sc(1);
  sd(1) = std(res);
  width = widths;
  
  for l = 2:Q
    th0 = sdang(l);
    thmn = th0-width/2;
    thmx = th0+width/2;
% Loop comps the weights in area around the midpoint
    if(thmn < 0)
      kmn = 2;
      idx = find(sdang<thmx);
      kmx = idx(end);
    elseif(thmx > sdang(end))
      idx = find(sdang>thmn);
      kmn = idx(1);
      kmx = Q;
    else
      idx = find(sdang>thmn);
      kmn = idx(1);
      idx = find(sdang<thmx);
      kmx = idx(end);
    end
    newwidth = sdang(kmx)-sdang(kmn);
    th = sdang(kmn:kmx);
    w = weight(th,th0,newwidth);
    w = w./sum(w);
    n(l) = length(w);
    scnew(l) = sum(sc(kmn:kmx).*w);
    if(scnew(l) >= scnew(1))
        scnew(l) = scnew(1)-l/Q*scnew(1);
    end
    sd(l) = std(sc(kmn:kmx))/sqrt(kmx - kmn);
    width = width+((widthe-widths)/Q);
    clear w;
  end % end l loop

%------------------------------------------------
%  Damping
%------------------------------------------------

fact = 10;
while(fact >= 1.)
  fact = fact-.1;
  x = real(cos(pi*[1:Q]./(2*Q)).^fact);
  dscnew = scnew .* x;
%
% Resort and fill in matrix
%
  [sIDX,IDX2] = sort(IDX);
  rdscnew = dscnew(IDX2);
  rscnew = scnew(IDX2);
  l = 2;
  for i = 1:nang
    dC_tmp(i,i) = rdscnew(1);
    C(i,i) = rscnew(1);
    for j = i+1:nang
      dC_tmp(i,j) = rdscnew(l);
      dC_tmp(j,i) = dC_tmp(i,j);
      C(i,j) = rscnew(l);
      C(j,i) = C(i,j);
      l = l+1;
    end
  end
  if(min(eig(dC_tmp))>0)
    dC = dC_tmp;
  else
    if(fact>3)
      load_fact = 0.001;
      dC_tmp2 = dC_tmp;
      while(min(eig(dC_tmp))<0)
        EE = load_fact*diag(ones(1,length(dC_tmp)));
        dC_tmp = (EE+ones(size(dC_tmp))).*dC_tmp2;
        load_fact = load_fact+0.001;
      end
      scnew(1) = (1+load_fact) * scnew_orig;
      if(load_fact >= .05)
        break;
      else
        dC = dC_tmp;
        fprintf(1,'\n WARNING:diag loading factor %10.6f\t%10.6f \n',...
                load_fact,fact);
      end
    else
      break;
    end
  end
end

return;

function [w] = weight(th,th0,newwidth)

w = cos((th-th0)*pi / 180. * 90. / newwidth);


return;
