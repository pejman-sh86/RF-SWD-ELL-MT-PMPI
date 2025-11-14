function z=runtest(x)
% PURPOSE: check if the positive and negative runs in the vector 
% x is random or not.
% -----------------------------
% USAGE: z=runtest(x)
% where: z is the z stat for a run
%        x=a vector variable (nobs x 1)
% -----------------------------
%
% written by:
% Wei Li
% MBA 2003, University of Chicago Graduate School of Business
% wayneli@fastmail.fm
% Jan 29, 2004
%
% The function could be used on any data with a binomial distribution.
% But the user need to assign positivity and negativity to the two 
% states in the data.
%
% Reference: Siegal (1956), Nonparametric Statistics.

x = x - median(x);

if size(x,1)<=20; error('Too few observations'); end;
if size(x,2)>1;error('the input data must be a vector');end;

logic=(x>0); %extra only positive and negative sign, 1 positive, 0 negative
run=diff(logic); %assign 1 or -1 to the beginning of the run.
run=abs(run);  %turn the negative run count into positive
run=[1;run]; %add a starting run, complete the run count.
n=size(logic, 1); n1=sum(logic); n2=n-n1; r=sum(run);
u_r=(2*n1*n2)/(n1+n2)+1;
std_r=sqrt((2*n1*n2*(2*n1*n2-n1-n2))/(((n1+n2)^2)*(n1+n2-1)));
z=(r-u_r)/std_r;

