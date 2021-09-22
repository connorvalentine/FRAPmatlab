function y = interstitial_hopping_fit_ln(x,B,G)
    y = log(B) - 2.*log(x) - G.*(x.^1.5);
     
end