function [] = mfi_test_residual_sample(represfile);

load(represfile);

NMOD  = size(zres,1);
NDAT2 = 4000;
NBAND = 9;
bands = [400., 450., 500., 550.,...
         600., 650., 700., 750., 800.];
%NBAND = 11;
%bands = [300., 350., 400., 450., 500., 550.,...
%         600., 650., 700., 750., 800.];

fid = fopen('residual_statistics.dat','w');
for ifr=1:NBAND
  zresp = squeeze(real(zres(:,ifr,:)));
  zresip = squeeze(imag(zres(:,ifr,:)));
  zresp_raw = squeeze(real(zresraw(:,ifr,:)));
  zresip_raw = squeeze(imag(zresraw(:,ifr,:)));

  meanzres=mean(zresp,2);
  meanzres=repmat(meanzres,1,size(zresp,2));
  zresp=zresp-meanzres;

  meanzresi=mean(zresip,2);
  meanzresi=repmat(meanzresi,1,size(zresip,2));
  zresip=zresip-meanzresi;

  meanzres=mean(zresp,2);
  meanzres=repmat(meanzres,1,size(zresp,2));
  zresp=zresp-meanzres;

  meanzresi_raw=mean(zresip_raw,2);
  meanzresi_raw=repmat(meanzresi_raw,1,size(zresip_raw,2));
  zresip_raw=zresip_raw-meanzresi_raw;

  ipass = 0;ifail = 0;iadpass = 0;iadfail = 0;
  ipassi = 0;ifaili = 0;iadpassi = 0;iadfaili = 0;
  ipass_raw = 0;ifail_raw = 0;iadpass_raw = 0;iadfail_raw = 0;
  ipassi_raw = 0;ifaili_raw = 0;iadpassi_raw = 0;iadfaili_raw = 0;

  xlim = 5.;
  dx = 0.2;
  edges = [-xlim-dx/2.:dx:xlim+dx/2.];
  ctrs = edges(1:end-1) + dx./2.;
  for i=1:NMOD;
 
    if(rem(i,1000) == 0);disp(i);end;

    [h(i),p(i)]=runstest(zresp(i,:),0.);
    if(h(i) == 1);ifail = ifail+1;end;
    if(h(i) == 0);ipass = ipass+1;end;

    [had(i),pad(i)]=AnDartest(zresp(i,:));
    if(had(i) == 1);iadfail = iadfail+1;end;
    if(had(i) == 0);iadpass = iadpass+1;end;

    [hi(i),pi(i)]=runstest(zresip(i,:),0.);
    if(hi(i) == 1);ifaili = ifaili+1;end;
    if(hi(i) == 0);ipassi = ipassi+1;end;

    [hadi(i),padi(i)]=AnDartest(zresip(i,:));
    if(hadi(i) == 1);iadfaili = iadfaili+1;end;
    if(hadi(i) == 0);iadpassi = iadpassi+1;end;

    [h_raw(i),p_raw(i)]=runstest(zresp_raw(i,:),0.);
    if(h_raw(i) == 1);ifail_raw = ifail_raw+1;end;
    if(h_raw(i) == 0);ipass_raw = ipass_raw+1;end;

    [had_raw(i),pad_raw(i)]=AnDartest(zresp_raw(i,:));
    if(had_raw(i) == 1);iadfail_raw = iadfail_raw+1;end;
    if(had_raw(i) == 0);iadpass_raw = iadpass_raw+1;end;

    [hi_raw(i),pi_raw(i)]=runstest(zresip_raw(i,:),0.);
    if(hi_raw(i) == 1);ifaili_raw = ifaili_raw+1;end;
    if(hi_raw(i) == 0);ipassi_raw = ipassi_raw+1;end;

    [hadi_raw(i),padi_raw(i)]=AnDartest(zresip_raw(i,:));
    if(hadi_raw(i) == 1);iadfaili_raw = iadfaili_raw+1;end;
    if(hadi_raw(i) == 0);iadpassi_raw = iadpassi_raw+1;end;

  end;


  runs_im(ifr)       = ipassi/(ipassi+ifaili);
  runs_real(ifr)     = ipass/(ipass+ifail);
  runs_im_raw(ifr)   = ipassi_raw/(ipassi_raw+ifaili_raw);
  runs_real_raw(ifr) = ipass_raw/(ipass_raw+ifail_raw);

  fprintf(fid,'\n');
  fprintf(fid,'FREQUENCY:         %12.4f %%\n',bands(ifr));
  fprintf(fid,'Runstest real:     %12.4f %% passed\n',ipass/(ipass+ifail));
  fprintf(fid,'Runstest imag:     %12.4f %% passed\n',ipassi/(ipassi+ifaili));
  fprintf(fid,'Runstest real raw: %12.4f %% passed\n',ipass_raw/(ipass_raw+ifail_raw));
  fprintf(fid,'Runstest imag raw: %12.4f %% passed\n',ipassi_raw/(ipassi_raw+ifaili_raw));

  fprintf(fid,'AD test real:      %12.4f %% passed\n',iadpass/(iadpass+iadfail));
  fprintf(fid,'AD test imag:      %12.4f %% passed\n',iadpassi/(iadpassi+iadfaili));
  fprintf(fid,'AD test real raw:  %12.4f %% passed\n',iadpass_raw/(iadpass_raw+iadfail_raw));
  fprintf(fid,'AD test imag raw:  %12.4f %% passed\n',iadpassi_raw/(iadpassi_raw+iadfaili_raw));
end;
fclose(fid);
save('runstest_results.mat','runs_im','runs_real','runs_im_raw','runs_real_raw');
return;
