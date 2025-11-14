function [] = label_subplot(fig,sub9,xlab,ylab,figtitle,np,nr,nc,...
              width,height,llh,llw);

figure(fig);

ii = 1;
for ip = 1:np

    if ip <= nc
      sp = [llw(ii) llh(1) width height];
      ii = ii + 1;ir = 1;
    elseif ip <= 2*nc
      sp = [llw(ii) llh(2) width height];
      ii = ii + 1;ir = 2;
    elseif ip <= 3*nc
      sp = [llw(ii) llh(3) width height];
      ii = ii + 1;ir = 3;
    else
      sp = [llw(ii) llh(4) width height];
      ii = ii + 1;ir = 4;
    end


    subplot(sub9(ip));
%    subplot('Position',sp);
    if (ip <= nc)
        title(figtitle);
    end
    if (ir == nr)
        xlabel(xlab);
    else
        set(gca,'XTickLabel',{});
    end
    if (ip == 1 | ip == nc+1 | ip == 2*nc+1 | ip == 3*nc+1)
        ylabel(ylab);
    else
        set(gca,'YTickLabel',{});
    end

end

return;
