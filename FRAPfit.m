%% fit function for FRAP
% need to add second half of this equation.
function y = FRAPfit(x,f,k,tau)
global Ifit
nsum = 3; % number of terms to evaluate sum to They use 8 in the paper
% y = zeros(length(x),1);

% for i = 1:length(x)
%      xi = x(i);
     syms n 
     ytemp = f*symsum(((-k).^n)./(1+n.*(1+ (2.*x./tau))*factorial(n)), n, 0, nsum)...
         +(1-f).*(Ifit);
     ytemp = double(ytemp);
     y = ytemp;
end