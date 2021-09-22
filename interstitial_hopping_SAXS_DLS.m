
function y = interstitial_hopping_SAXS_DLS(x,B,G)
    global d a 
    
    y = B + log(a.^2) - G.*(d.^2).*(x.^(2/3));
     
end