%% Assets:
acc.marg = [150000., 127000.];
acc.rrsp = [170000., 0.];
acc.tfsa = [54000., 50000.];
acc.cash = [32000., 10000.];

%% House:
h.price = 410000.;  %% Purchase price
h.pay_yr = 40000;   %% Payement per year
h.term = 5;

%% Taxes
tax.marg = 0.36;

%% Interest & returns
int.prime = 0.027;
int.mort = 0.0184;

int.mort_m = int.mort/12;
int.heloc = int.prime + 0.005;

h.down = sum(acc.marg)+sum(acc.cash);
mort.tot = h.price - h.down;
h.term_m = h.term*12;
heloc.max = 0.65 * h.price;
heloc.ltv = 0.80 * h.price;
heloc.start = heloc.ltv - mort.tot;

months = [1:12*h.term];
years = months/12.;

mort.pay    = int.mort_m*mort.tot/(1-(1+int.mort_m)^(-h.term_m));
h.lump = h.pay_yr - mort.pay*12;
h.int_totla = h.term_m*mort.pay-mort.tot;

mort.bal(1)    = mort.tot;
mort.int(1)    = mort.bal(1)*int.mort_m;
mort.prin(1)   = mort.pay - mort.int(1);
heloc.bal(1)   = heloc.start;
heloc.int(1)   = heloc.bal(1)*int.heloc/12.;
mort.bal_acc(1)= mort.tot;

nw(1,:) = [years(1),months(1),mort.int(1),mort.prin(1),mort.bal(1)];


for im = 2:h.term_m;

    mort.bal(im) = mort.bal(im-1) - mort.prin(im-1);
    mort.int(im)  = mort.bal(im)*int.mort_m;
    mort.prin(im) = mort.pay - mort.int(im)

    if(heloc.ltv-mort.bal(im-1)<heloc.max);
        heloc.bal(im) = heloc.ltv-mort.bal(im-1);
    else;
        heloc.bal(im) = heloc.max;
    end;
    heloc.int(im) = heloc.bal(im)*int.heloc/12.;
    heloc.break(im) = 0.;
    if(mod(im,8)==0);
        lag = 11;
        if(im == 8);lag = 7;end;
        heloc.break(im) = sum(heloc.int(im-lag:im))*tax.marg;
    end;
    mort.lump(im) = 0.;
    if(mod(im,12)==0);
        mort.lump(im) = h.lump;
    end;
    mort.bal_acc(im) = mort.bal_acc(im-1) - mort.prin(im-1) ...
                     - heloc.break(im) - mort.lump(im);
    
    nw(im,:) = [years(im),months(im),mort.int(im),mort.prin(im),mort.bal(im)];
    
end;


