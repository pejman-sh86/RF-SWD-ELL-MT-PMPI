%----------------------------------------------------------
% Compare Prior from PPD
%----------------------------------------------------------
function [] = compare_prior();

%prior = {'sim_A_1dB_1_prior.mat' 'sim_A_1dB_2_prior.mat' ...
%          'sim_A_1dB_3_prior.mat'};
prior = {'sim_A_1dB_4_nols_prior.mat' 'sim_A_1dB_4_prior.mat'};
%prior = {'sim_A_1dB_1_prior.mat'};
prior = char(prior);
nsubfig = 12;
ncomp = size(prior,1)
%P(1).idx = [1 2 3 4 6 7 8];
%P(2).idx = [1 2 3 4 5 6 7 8 9 10 11 12 2 3 4];
%P(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 10 11 12];
%P(2).idx = [1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 10 11 12];
P(1).idx = [1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 9 10 11 12 2 3 4];
P(2).idx = [1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 9 10 11 12 2 3 4];
%P(3).idx = [1 2 3 4 5 6 7 8 9 10 11 12 1 2 3 4 5 6 7 8 9 10 11 12 2 3 4];


mtrue = [1.5     1540    1.40    0.05 3.0     1520    1.55    0.10 1.0     1560    1.75    0.15 7.0     1600    1.90    0.10 1.0     1700    2.05    0.10 4.0     1600    2.10    0.10 1650    2.20    0.01];
axisl = {'h_1' 'c_1' '\rho_1' '\alpha_1' 'h_2' 'c_2' '\rho_2' '\alpha_2' 'h_3' 'c_3' '\rho_3' '\alpha_3' 'h_4' 'c_4' '\rho_4' '\alpha_4' 'h_5' 'c_5' '\rho_5' '\alpha_5' 'h_6' 'c_6' '\rho_6' '\alpha_6' 'c_b' '\rho_b' '\alpha_b' };

minlim = [0.5     1500    1.20    0.00 2.0     1500    1.2    0.00 0.0     1500    1.2    0.0 6.0     1550    1.20    0.00 0.0     1650    1.2    0.00 3.0     1550    1.20    0.00 1600    1.20    0.00];
maxlim = [2.5     1600    2.30    0.5  4.0     1600    2.3    0.50 2.0     1600    2.3    0.5 8.0     1650    2.30    0.50 2.0     1750    2.3    0.50 5.0     1650    2.30    0.50 1700    2.30    0.5];

for i = 1:ncomp

    load(prior(i,:));
    P(i).lim(:,:) = lim;
    P(i).width = lim(end,:)-lim(1,:);
    P(i).delta = lim(3,:)-lim(2,:);
%    P(i).n1(:,:) = n1.*repmat(P(i).width,size(P(i).lim,1),1)./repmat(P(i).delta,size(P(i).lim,1),1);
    P(i).n1(:,:) = n1./repmat(P(i).delta,size(P(i).lim,1),1);
    P(i).n1max = max(P(i).n1,[],1);

end

save bla P;

%load(prior2);
%P(2).lim(:,:) = lim;
%P(2).n1(:,:) = n1;
clear lim; clear n1;
col = {'-b' '-k' '-r'};

col = char(col);

for i = 1:ncomp;
    
    npar(i) = size(P(i).lim,2);
    nbin(i) = size(P(i).n1,1);
    P(i).n1  = cat(1,zeros(1,npar(i)),P(i).n1,zeros(1,npar(i)));
    P(i).lim = cat(1,P(i).lim(1,:)-(P(i).lim(2,:)-P(i).lim(1,:)).* ...
               ones(1,npar(i)),P(i).lim,P(i).lim(nbin(i),:)+ ... 
               (P(i).lim(nbin(i),:)-P(i).lim(nbin(i)-1,:)).* ... 
               ones(1,npar(i)));
end

%----------------------------------------------------------
%  Plot histograms:
%----------------------------------------------------------


for i = 1:ncomp;

    nfig(i) = ceil(npar(i)/nsubfig);
    ipar = 1;

    for ifig = 1:nfig(i)

        figure(ifig); hold on; box on;
        isubfig = 1; 

        while (isubfig <= nsubfig & ipar <= npar(i))

            subplot(3,4,P(i).idx(ipar));hold on; box on;

            stairs(P(i).lim(:,ipar),P(i).n1(:,ipar),col(i,:));
            if(i == ncomp)
                itrue = ipar;
                if(ifig > 1)
                    itrue = ipar;
                elseif(ifig > 2)
                    itrue = ipar + 2;
                end
                plot([mtrue(itrue) mtrue(itrue)],[0 P(i).n1max(ipar)],'--k');
                set(gca,'XLim',[minlim(ipar) maxlim(ipar)]);
                xlabel(axisl(ipar));
%                set(gca,'YLim',[0 P(i).n1max(ipar)]);
            end

            ipar = ipar +1;
            isubfig = isubfig + 1;
 
        end;

    end; 
end;

return;
