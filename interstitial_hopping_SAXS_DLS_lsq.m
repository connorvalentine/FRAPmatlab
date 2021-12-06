function y = interstitial_hopping_SAXS_DLS_lsq(x,xdata)
    % d and csi are in nm.
    % a should be in microns for the prefactor so D is um^2/s
    global d pluronic_name T
    y_temp = zeros(length(xdata),1);
    B = x(1);
    G = x(2);
    
    rho = 1e-21; %g/nm^3 in the system
    rs = 0.1057; % nm^3/monomer of PPO
    Navg = 6.022e23;
    
    if strcmp(pluronic_name,'F87')
        MW = 7700; %g/mol
        PEOm = 61;
        PPOn = 40;
        
    elseif strcmp(pluronic_name,'F127')
        MW = 12600; %g/mol
        PEOm = 99;
        PPOn = 69;
    end
    for i = 1:length(xdata)
        a_calc = unit_cell_extrapolation(pluronic_name,T,xdata(i)*100); % nm
        nagg = rho*Navg*(a_calc^3)*xdata(i)/(2*MW);
        Vmicelle = (4/3)*pi()*((sqrt(3)/4)*a_calc)^3;
        Vcore = nagg*PPOn*rs;
        csi = ((Vmicelle-Vcore)./(2*nagg*PEOm)).^(1/3);
        
        y_temp(i) = B + log((a_calc/1000).^2) - (G.*(d./csi).^2)./(T+273);
    end
    y = y_temp;
end