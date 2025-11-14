%Processing for Reflection coefficients as a function of angle
% Does both 1/3 octave and Narrow band processing
%
% must set freq, ID (add NB for narrow band processing), x2, tlim
% if the swell filter has been run, it replaces the time, range and angle fields with the filtered results at the end of the processing

%freq=[100:50:2000 2050:100:10000];
freq=-[20:40]; %for third ctave processing the frequencies should be negative

%INPUT ID#
         ID='id4_neg'

%FIRST TIME get rid of any traces that are clipped or have non-monotonic ranges
%x2 = xremove_post(x, union(x.hi,kki) ); x2.hi=[];kki=[];
%x2=x;

% DETERMINE DUSS Sensitivity; started using DUSS in Oct 1999
 hyd_sens=-201; clear year; 
if exist('fnames','var'); ii=max(findstr(fnames(1,:),'duss')); 
    if length(ii)>0; year=fnames(1,ii+4:ii+7);
      else ii=max(findstr(fnames(1,:),'vla')); year=fnames(1,ii+3:ii+6);end
  else year=input('Enter year \n','s'); end

duss_sens=0;  
% if year=='1997' | year=='1998';duss_sens=0; end
 if year=='1999';  vla=input('enter 1 if Jan 1999 or 2 if Sept-Oct 1999');
     if vla==2; duss_sens=-192; end; end
 if year=='2000' | year=='2001'; duss_sens=hyd_sens+17+23.5; end; %NB: depending on phone there may be another gain
 if year=='2002' | year=='2004'; duss_sens=hyd_sens+21.5;  end

block_len=x2.fout/24e3*8192; %set block length

dr1='d:\experiment_data\';
%dirnm=[dr1 'Capraia SCARAB97\bottom loss\Site 2\'];
%dirnm=[dr1 'Capraia SCARAB97\bottom loss\Site 3\'];
%dirnm=[dr1 'Capraia SCARAB97\bottom loss\Site 5\'];

dirnm=[dr1 'Sicily SCARAB98\BL_VLA\Site2\'];tlimd=[3e-3, 17e-3];ampd=1000;tlimb=[14e-3, .035]; ampb=2500;
%dirnm=[dr1 'Sicily SCARAB98\BL_VLA\Site4\'];tlimd=[4e-4, 6e-3];ampd=1000;tlimb=[8e-3, .011]; ampb=3000;
%dirnm=[dr1 'Sicily SCARAB98\BL_HLA\UBDW3\'];

 %dirnm=[dr1 'boomer99\bottom loss\BLGW2 Site 2 Jan 29,99\'];
 %dirnm=[dr1 'boomer99\bottom loss\Site 2 Jan 26,99\'];
 %dirnm=[dr1 'boomer99\bottom loss\Site5\'];
 %dirnm=[dr1 'boomer99\bottom loss\Site1\'];
  %dirnm=[dr1 'Malta Geo99\bottom loss\Site6\'];
 %dirnm=[dr1 'Geoscat99\bottom loss\Site FVLA\'];
 %dirnm=[dr1 'Geoscat99\bottom loss\Site 4\'];

 %dirnm=[dr1 'Boundary 2000\bottom loss\Capraia_Site6\'];
 %dirnm=[dr1 'Boundary 2000\bottom loss\Site7\'];
 %dirnm=[dr1 'Boundary 2000\bottom loss\Site8\'];

 %dirnm=[dr1 'geoclutter2001\bottom loss\Site2b\'];
 %dirnm=[dr1 'geoclutter2001\bottom loss\Site2b_N\'];
 %dirnm=[dr1 'geoclutter2001\bottom loss\Site 1A\'];
 %dirnm=[dr1 'Boundary01\bottom loss\N1A\'];
 %dirnm=[dr1 'Boundary01\bottom loss\N2C_C\'];
 %dirnm=[dr1 'Boundary01\bottom loss\N3\'];
 %dirnm=[dr1 'Boundary01\bottom loss\N4\'];
 %dirnm=[dr1 'Boundary01\bottom loss\N6\'];
 %dirnm=[dr1 'Boundary01\bottom loss\N8\'];
 %dirnm=[dr1 'Boundary01\bottom loss\S1\'];
 %dirnm=[dr1 'Boundary01\bottom loss\S2\'];

%dirnm=[dr1 'Boundary02\bottom loss\Site 10\'];
%dirnm=[dr1 'Boundary04\bottom loss\Site 20\'];
 
prfx=[dirnm 'bl data\hp_results\bl' blnum '_' dirc ph];
rng2=(x2.r/x2.c).^2; 

for iw=1:-1:0
 clear fn
 if iw==1 ; 
'Pick bottom reflection limits of integration'
    fn='b'; 
    figure;box on
     xreduce_sqrt(x2, sqrt(tlimb), ampb);hold on; plot(x2.r, x2.t.^2-rng2,'c:')
     tsq_lim1 = xhyp(x2);
     tsq_lim1 = tsq_lim1.^2;

  else if iw==0 ;      
  'Pick direct path limits of integration (actually these are automatically set so they are consistent with direct)'
    fn='d'; 

 % first set direct path limits based on bottom path int time
    winstart_off=5e-5;
    rng2_r=repmat(rng2,2,1);
    tb=sqrt(tsq_lim1+rng2_r); int_time_bot=diff(tb); 
    tsq_limd(1,:)=x2.t.^2-rng2-winstart_off;
    tsq_limd(2,:)=(x2.t+int_time_bot).^2-rng2-winstart_off;
    figure; box on
    xreduce_sqrt(x2, sqrt(tlimd), ampd);hold on; plot(x2.r,sqrt( x2.t.^2-rng2),'c:')
    plot(x2.r,sqrt(tsq_limd)','m:');
    ans1=input('is the fit acceptable? (y or n)')
    if ans1=='y' ;
        tsq_lim1=tsq_limd;
      else
        tsq_lim1 = xhyp(x2); 
         tsq_lim1 = tsq_lim1.^2; 
     end
                 
 else;end
end    
    
% COMPUTE RECIEVED LEVEL (xoct2... uses a long time window for noise, xoct1 uses the same length as the signal)
     x1 = xoct2_all(x2, tsq_lim1, freq, 1e-2, block_len,duss_sens); %DUSS
     title(['bl' blnum '\_' dirc ph '\_rl' fn '\_' ID]) 

% WRITE OUT DIRECT PATH DATA or reflection data
  if exist('swell','var'); x1.r = x2.r_swell;x1.ang = x2.ang_swell;x1.t = x2.t_swell;end

 if freq<0; 
   octpr([prfx '_rl' fn '_' ID '.asc'], x1)
   snrpr([prfx '_rl' fn 'snr_' ID '.asc'], x1)
   save([prfx '_t' fn 'sqr_' ID], 'tsq_lim1')
  else
   save([prfx '_' fn '_' ID], 'x1')
 end
   
end