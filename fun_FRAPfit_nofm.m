%% fit function for FRAP
% need to add second half of this equation.
function y = fun_FRAPfit_nofm(x,k,tau_n)
global Ifitparam fm0
nsum = 12; % number of terms to evaluate sum to They use 8 in the paper
     tau = tau_n;
     syms n 
     ytemp = fm0*symsum(((-k).^n)./(1+n.*(1+ (2.*x./tau))*factorial(n)), n, 0, nsum)...
         +(1-fm0).*(Ifitparam);
     ytemp = double(ytemp);
     y = ytemp;
end