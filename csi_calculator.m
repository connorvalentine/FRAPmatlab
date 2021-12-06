function [a,naggout,csi,vm, vcore] = csi_calculator(pluronic_name,xdata)

    y_temp1 = zeros(length(xdata),1);
    y_temp2 = y_temp1;
    y_temp3 = y_temp1;
    y_temp4 = y_temp1;
    y_temp5 = y_temp1;

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
        a_calc = unit_cell_extrapolation(pluronic_name,25,xdata(i)*100); % need this nm
        nagg = rho*Navg*(a_calc^3)*xdata(i)/(2*MW);
        Vmicelle = (4/3)*pi()*((sqrt(3)/4)*a_calc)^3;
        Vcore = nagg*PPOn*rs;
        csi = ((Vmicelle-Vcore)./(2*nagg*PEOm)).^(1/3);
        
        y_temp1(i) = a_calc;
        y_temp2(i) = nagg;
        y_temp3(i) = csi;
        y_temp4(i) = Vmicelle;
        y_temp5(i) = Vcore;
    end
    a = y_temp1;
    naggout = y_temp2;
    csi = y_temp3;
    vm = y_temp4;
    vcore = y_temp5;
end