function [] = tornado();

%FileRoot='sim_A_2b_PDE';
%FileRoot='sim_A_2_4_16_3lay_PDE';
%FileRoot='sim_B6_1_7_40_3layrg';

%FileRoot='x_s16_1_17_50_3layrg';
FileRoot='x_s16_1_15_40_3layrg';
%FileRoot='x_s16_1_7_40_3layrg';
%FileRoot='x_s16_2_7_32_4layrg';

% Read variables from the .m file
FileLog=[FileRoot '_log'];
eval(FileLog);

FileName=[FileRoot '_hist.dat'];
B = load(FileName);
B(:,1)= B(:,1)+abs(min(B(:,1)))+1;
ebest = ebest + abs(ebest)+1;

npx  = size(B,1);   % number of iterations
nmod = size(B,2)
 
ymin = min(B(:,1));
ymax = max(B(:,1));

figure(1);
for imod=2:nmod
   subplot(5,4,imod-1);
   hold on;

   plot(B(:,imod),B(:,1),'.k','MarkerSize',0.5);      
   plot( mbest(imod-1), ebest,'+r','MarkerSize',10);
%   plot([mtrue(imod-1),mtrue(imod-1)],[ymin ymax],'--b');

   set(gca,'YScale','log');
   set(gca,'XLim',[minlim(imod-1) maxlim(imod-1)] );
   set(gca,'YLim',[ymin ymax]);
end;

return;
