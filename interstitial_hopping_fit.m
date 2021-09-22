function y = interstitial_hopping_fit(x,B,G)
    y = B.*(x.^(-2)) .* exp(-G .* (x.^(1.5)));
    disp(y)
     
end