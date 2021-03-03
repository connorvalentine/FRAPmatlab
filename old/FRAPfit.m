%% fit function for FRAP
% need to add second half of this equation.
function y = FRAPfit(x,f,k,tau)
global Ifitparam
nsum = 4; % number of terms to evaluate sum to They use 8 in the paper

     syms n 
     ytemp = f*symsum(((-k).^n)./(1+n.*(1+ (2.*x./tau))*factorial(n)), n, 0, nsum)...
         +(1-f).*(Ifitparam);
     ytemp = double(ytemp);
     y = ytemp;
end