function [nf_int] = hpd2(n1,edges,percent)

n_ci = sum(n1)*percent/100.;
if(max(n1) > 0);
%% Find first element:
stopi = 0;
for i=1:length(n1);
  if(stopi == 0);
    if(n1(i) >0);
      lb(1) = i;
      lbmin = i;
      stopi = 1;
    end;
  end;
end;
%% Find last element:
stopi = 0;
for i=0:length(n1)-1;
  if(stopi == 0);
    if(n1(end-i) >0);
      rbmax = length(n1)-i;
      stopi = 1;
    end;
  end;
end;

%% Find first X % interval:
stopi = 0;
for i=lb(1):length(n1);
  if(stopi == 0);
    if(sum(n1(lb(1):i)) >= n_ci);
      rb(1) = i;
      stopi = 1;
    end;
  end;
end;
stopj = 0;
for j = 2:rbmax;
  if(stopj == 0);
    lb(j) = lb(j-1) + 1;
    stopi = 0;
    for i=lb(j):rbmax;
      if(stopi == 0);
        if(sum(n1(lb(j):i)) >= n_ci);
          rb(j) = i;
          stopi = 1;
        end;
      end;
    end;
    if(stopi == 0);rb(j) = rbmax;end;
    if(rb(j) >= rbmax);
      stopj = 1;
    end;
  end;
end;

bdiff = rb-lb;

stopj = 0;
[foo(1),nf_index2(1)] = min(abs(bdiff));
bdiff(nf_index2(1)) = 1e10;
for j = 2:length(rb);
  if(stopj == 0);
    [foo(j),nf_index2(j)] = min(abs(bdiff));
    bdiff(nf_index2(j)) = 1e10;
    if(foo(j)>foo(j-1));
      stopj = 1;
      foo(j) = [];
      nf_index2(j) = [];
    end;
  end;
end;

nf_int(1) = mean(edges(lb(nf_index2)));
nf_int(2) = mean(edges(rb(nf_index2)+1));
else;
nf_int(1) = edges(1);
nf_int(2) = edges(1);
end;
%figure();
%stairs(edges,n1);hold on;
%plot([nf_int(1) nf_int(1)],[0  300],'k');
%plot([nf_int(2) nf_int(2)],[0  300],'k');
%plot([edges(lb(1)) edges(lb(1))],[0  300],'--r');
%plot([edges(rb(1)) edges(rb(1))],[0  300],'--r');
%plot([edges(rbmax) edges(rbmax)],[0  300],'--y');
%set(gca,'XLim',[edges(lbmin) edges(rbmax)]);
return;
