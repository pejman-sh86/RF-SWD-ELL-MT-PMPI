%%
%% Performs ping averaging and computed sigma for AUV
%%
function [] = auv_put_noise();

isave = 1;
bands = [975. 1100. 1250. 2100. 2400. 2700.];
sig = 0.05;
nf = length(bands);
pstgl  = 1;
pendgl = 170;
fileext  = '.txt';
fileext2  = '.png';

for iping=1:pendgl;
   disp(iping);

   filebase1 = 'p';
   filebase2 = strcat('_',num2str(bands(1),'%04i'),...
                      '_',num2str(bands(end),'%04i'));
   outfile1 = strcat(filebase1,num2str(iping,'%03i'),filebase2,fileext);
   outfile2 = strcat(filebase1,num2str(iping,'%03i'),filebase2,fileext2);

   dat(:,:,iping) = dlmread(outfile1);
   NANG = size(dat,2);
   noise = sig*randn(nf,NANG);
   sim(:,:,iping) = dat(:,:,iping);
   sim(5:4+nf,:,iping) = dat(5:4+nf,:,iping)+noise;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%
   %% Save data
   %%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   z_t = dat(1,1,iping);
   cw  = dat(2,1,iping);
   rw  = dat(3,1,iping);
   hmx = dat(4,1,iping);
   save(outfile1,'z_t','-ascii');
   save(outfile1,'cw','-ascii','-append');
   save(outfile1,'rw','-ascii','-append');
   save(outfile1,'hmx','-ascii','-append');
   for j = 1:nf
      tmp = squeeze(sim(j+4,:,iping));
      save(outfile1,'tmp','-ascii','-append');
      clear tmp;
   end
   tmp = dat(4+nf+1,:,iping);
   save(outfile1,'tmp','-ascii','-append');
   clear tmp;
  if(isave == 1);
    fig = figure(1);hold off;
    for j = 1:nf
      subplot(2,3,j);hold off;
      plot(dat(4+nf+1,:,iping),dat(4+j,:,iping),'-b');hold on;
      plot(sim(4+nf+1,:,iping),sim(4+j,:,iping),'xk');
      set(gca,'XLim',[25,70],'YLim',[0 1]);
      hold off;
    end
    print(fig,outfile2,'-dpng');
  end;

end;
save tmp.mat dat sim;

return;
