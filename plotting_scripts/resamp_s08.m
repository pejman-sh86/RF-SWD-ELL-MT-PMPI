function [] = resamp(filename);

filebase = strrep(filename,'.mat','');

i_inter = 0;  %% 0 -> leave raw
              %% 1 -> linear interpolation
              %% 2 -> bin averaging
nang_re = 90;
load(filename);

%% Frequencies:
%%1000 1200 2000 2200 2400 2600 2800 3000 3200
%idx = [9,11,13,15];
%fidx = [26:5:81]
%fidx = [26:5:56]
%fidx = [[26:5:56],[26:5:56]+2];
%fidx=sort(fidx);
%fidx = [26:5:61]
%fidx = [46:5:56]
%fidx = [[26:5:81],[26:5:76]+2];
%fidx=sort(fidx);
%fidx = [36:5:76]
%% Elba old data (1/3 octave from 1997)
fidx = [9:1:15]

freq = round(x.pref(fidx,1))

nf = length(fidx);
theta = x.ang(1:1:end);
R = x.R(fidx,1:1:end);
%R = R';
%R = (x.R(fidx,:)+x.R(fidx2,:))/2.;
z_t = x.z_t;
cw = x.cw;
rw = 1.029;
hmx = x.hmax;

outfile  = strcat(filebase,'_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.mat')
outfilea  = strcat(filebase,'_',num2str(freq(1),'%04i'),'_',num2str(freq(end),'%04i'),'.txt')

thdiff = diff(theta);
if(i_inter > 0);
   dtheta = (theta(end)-theta(1))/(nang_re-1);
   theta_even = [theta(1):dtheta:theta(end)];
else
   theta_even = theta;
end;
if(i_inter == 2);
  theta_lo = theta_even-dtheta/2;
  theta_hi = theta_even+dtheta/2;
end;

if(i_inter == 1);
  for ifr = 1:nf
    R_even(ifr,:) = interp1(theta,R(ifr,:),theta_even);
  end
elseif(i_inter == 2);
  for ifr = 1:nf
    for iang = 1:nang_re
      idx = find(theta>theta_lo(iang) & theta<theta_hi(iang) );
      if(idx);
        R_even(ifr,iang) = nanmean(R(ifr,idx));
      else;
        R_even(ifr,iang) = NaN;
      end;
    end
  end
else;
  for ifr = 1:nf
    R_even(ifr,:) = R(ifr,:);
  end
end;
nanidx = isnan(R_even);
Rex = (-1*isnan(R_even))+1;
R_even(nanidx) = 0.;

save tmp.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Save data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F(1).cw = cw
F(1).rw = rw;

F(1).freq = x.pref(fidx,1)
for ifr = 1:nf;
    F(ifr).dat = R_even(ifr,:);
    F(ifr).ang = theta_even;
end
F(1).Rex = Rex;


save(outfile,'F');

filebase = strrep(outfile,'.mat','');
fid = fopen('filebase.txt','w');
filebase
fprintf(fid,'%i\n',length(filebase));
fprintf(fid,'%s',filebase);
fclose(fid);

save(outfilea,'z_t','-ascii');
save(outfilea,'cw','-ascii','-append');
save(outfilea,'rw','-ascii','-append');
save(outfilea,'hmx','-ascii','-append');
f1=figure();
f2=figure();
for i = 1:nf

   tmp = F(i).dat;
   save(outfilea,'tmp','-ascii','-append');
   tmp = F(i).dat;
   tmp(find(Rex(i,:) == 0)) = NaN;
   figure(f1);
   subplot(4,4,i);hold on;
   plot(F(i).ang,tmp,'.-k');
   figure(f2);
   subplot(4,4,i);hold on;
   plot(tmp,'.-k');

end
save(outfilea,'theta_even','-ascii','-append');

for i = 1:nf

    tmp = F(1).Rex(i,:);
    save(outfilea,'tmp','-ascii','-append');

end

return;
