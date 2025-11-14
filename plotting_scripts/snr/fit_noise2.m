function [ins2,ysnrall] = fit_noise2(tldsnr,frq_column,snr,minrange,maxrange,rangemax)
%function [ins2,ysnrall] = fit_noise2(tldsnr,frq_column,snr,minrange,maxrange,rangemax)
%
% this fitting uses a sliding window to smooth the data before screening for SNR
%
num_smooth_points=7;

 inan=find(isnan(tldsnr(:,frq_column ))==0 ); %screen out NaNs and ranges where SNR is unreliable
 [ysnr] = sliding_win(tldsnr(inan,frq_column),0,num_smooth_points);
 ysnrall=ones(size(tldsnr,1),1)*nan; ysnrall(inan)=ysnr;
 [ins2]=find(ysnrall>snr & tldsnr(:,frq_column)>snr &...
    tldsnr(:,1)>minrange & tldsnr(:,1)<maxrange & tldsnr(:,1)<rangemax);
