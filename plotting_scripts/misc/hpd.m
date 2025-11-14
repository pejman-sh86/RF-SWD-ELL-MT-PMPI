function [nf_int] = hpd(m,nrbin,percent)

percent = 100-percent;
xmin = min(m);
xmax = max(m);
xdiff = xmax - xmin;
for i =0:nrbin
  edges(i+1) = xmin + xdiff/nrbin * i;
end

n1 = histc(m,edges);
nintyfive = sum(n1)-percent*sum(n1)/100;
lb = 0;
rb = 0;

p = 1;
for i = 1:length(n1)
  s = 0;
  stopj = 0;
  for j = i:length(n1)
    if(stopj == 0)
      s = s + n1(j);
      if(s >= nintyfive)
        lb(p) = edges(i);
        rb(p) = edges(j);
        p = p+1;
        stopj = 1;
      end
    end
  end
end
bdiff = lb-rb;
[foo,nf_index] = min(abs(bdiff));
nf_int(1) = lb(nf_index);
nf_int(2) = rb(nf_index);

return;
