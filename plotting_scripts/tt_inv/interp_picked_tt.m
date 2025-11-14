function [] = interp_picked_tt()

tt_file = 'picked_tt.mat';
tt_file2 = 'picked_tt_int.txt';
tt_file3 = 'picked_tt_int.mat';
NLAY = 6;

load(tt_file);
r = tt(end,1:end)
tt = tt(1:end-1,1:end);
r_old = r;

for iter = 1:5
dr = diff(r)
for ir = 1:length(dr);

  if(dr(ir) > 8.0)
 
    nr = floor(dr(ir)/(dr(ir-1)+1))
    clear r_ins;
    for iins = 1:nr;
      r_ins(iins) = r(ir)+iins*dr(ir-1);
    end
    r = [r(1:ir) r_ins r(ir+1:end)];
    break;

  end

end
end
%plot(diff(r))

for ilay = 1:NLAY

    tt_int(ilay,:)   = interp1(r_old,tt(ilay,:),r);
    figure(1);hold on;box on;
    plot(r,tt_int(ilay,:),'+');
    plot(r_old,tt(ilay,:),'+r');

end

tt = tt_int;
save(tt_file2,'tt','r','-ascii');
save(tt_file3,'tt','r');

return;
