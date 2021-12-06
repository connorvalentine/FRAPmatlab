function unit_cell_size = unit_cell_extrapolation(pluronic,T,c)

if strcmp(pluronic,'F87')
    if T == 25
        a = @(c) 0.0015675945.*(c.^2) - 0.1224856646.*c + 17.2532184510;
    elseif T == 35
        a = @(c) 0.0040238558.*(c.^2) - 0.3061125063.*c + 20.7590417232;
    elseif T == 45
        a = @(c) -0.0008388480.*(c.^2)+ 0.0593267649.*c + 14.3591039139;
    elseif T == 55
        a = @(c) -0.0058525633.*(c.^2) + 0.4328334281.*c + 7.8087166053;
    end
end

if strcmp(pluronic,'F127')
    if T == 25
        a = @(c) -0.2507420576.*c + 28.3759718341;
    elseif T == 35
        a = @(c) -0.2490817260.*c + 28.3801361940;
    elseif T == 45
        a = @(c) -0.2470966047.*c + 28.3698931315;
    elseif T == 55
        a = @(c) -0.2532182270.*c + 28.5782844159;

    end
end

unit_cell_size = a(c);
end