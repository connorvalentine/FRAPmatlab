% fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam
function [struct_out] = fun_intensity_fits(fits,alldata)

for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    disp(position)
    % adding fit parameters 
    global Ifitparam
        % Ifitparam is basically the initial intensity (normalized) of the
        % bleached region (first frame after bleaching).
        
        % find min indices - for some reason there is an initial drop in
        % intensity sometimes. This is due to poor laser focusing sometimes
        % where the profile is not guassian
        a = [alldata.(position).IN_ti];
        ind = find(a == min(a));
        I_t1 = alldata.(position)(ind).I_ti ;
        R_t1 = alldata.(position)(ind).ref_i;
        I_t0 = fits.(position).I_t0;
        R_t0 = fits.(position).ref_0;
        Ifitparam = (I_t1./I_t0).*(R_t0./R_t1);
    
    % use x for the variable instead of t, works for fit algorithm better.
    x = [alldata.(position).dt]';
    x = x(ind:end); % could have to adapt this later
    norm_i = [alldata.(position).IN_ti]';
    norm_i = norm_i(ind:end);
    
    % remove any bad frames where normalized intensity is not a number
    NANind = find(isnan(norm_i));
    norm_i(NANind) = [];
    x(NANind) = [];
    
    t_fit = linspace(0,x(end),250)'; % fake time data to put into the fit equation
    
    % fm is the mobile fraction of proteins. 
    % in this case, it is taking as the ratio of the amount of intensity
    % recovered (from first frame after photobleacing aka t1 to tend)
    % borrom of ratio is the amount of intensity lost from the
    % photobleaching.

    fm0 = (norm_i(end)-norm_i(1))/(1-norm_i(1));
    fm0_lb = 0.9*fm0;
    fm0_ub = 1.1*fm0;  
    if fm0_ub > 1
        fm0_ub = 1;
    else
    end

    % using this as a first guess point for the coefficients
    ft = fittype('fun_FRAPfit(x,f,k,tau_n)');
        options = fitoptions(ft);
        options.StartPoint = [fm0,1.3, 3];
        options.DiffMinChange = 0.001;
        options.TolFun = 1e-8;
        options.Algorithm = 'Levenberg-Marquardt';
        options.Robust = 'off';  
    % perform the fit to the data
    [f2,gof2,output2] = fit(x,norm_i, ft,options);
    
    % double fit process to place bounds on f.
    dlower = 0.5;
    dupper = 2;
    ft2 = fittype('fun_FRAPfit(x,f,k,tau_n)');
        options2 = fitoptions(ft2);
        options2.StartPoint = [fm0   ,f2.k     ,f2.tau_n    ];
        options2.Lower =      [fm0_lb,0.9*f2.k ,dlower*f2.tau_n];
        options2.Upper =      [1     ,1.1*f2.k ,dupper*f2.tau_n];
        options2.DiffMinChange = 0.001;
        options2.TolFun = 1e-8;
        options2.Robust = 'off';  
    % perform the fit to the data
    [f,gof,output] = fit(x,norm_i, ft,options2);    
    I_fit = f(t_fit);   
%     
    % confidence intervales from the fit, f
    ci = confint(f,0.95);
    x_ci = t_fit;
    y_ci = predint(f,x_ci,0.95);
    %shapec and xc can be used to plot the shaded confidence region
    xc = [x_ci;flip(x_ci)];
    curve1 = y_ci(:,1);
    curve2 = y_ci(:,2);
    shapec = [curve1; flip(curve2)];
    
   % calculated diffusivity value
    ri = fits.(position).FWHM .*0.5 .*fits.(position).pixel_size;
    dr = (fits.(position).dFWHM(2) - fits.(position).FWHM).*0.5 .*fits.(position).pixel_size;
    tau = 1000*f.tau_n;
    dtau = 1000*ci(:,end);
    D = (ri.^2)./(4*tau);
    Dlow = ((ri -dr).^2)./(4*max(dtau)); 
    Dhigh = ((ri+dr).^2)./(4*min(dtau)); 
    Dneg = D-Dlow;
    Dpos = Dhigh-D;
    errD = [Dneg Dpos];
    
    % adding the fit information to fits structure
    fits.(position).('fit_info') = f;
    fits.(position).('ci') = ci;
    fits.(position).('D') = D;
    fits.(position).('errD') = errD;
    fits.(position).('gof') = gof;
    fits.(position).('I_fit') = I_fit;
    fits.(position).('t_fit') = t_fit;
    fits.(position).('xc') = xc;
    fits.(position).('shapec') = shapec; 
    fits.(position).('Ifitparam') = Ifitparam;
    fits.(position).('fm0') = fm0;
    
     % output values from fit 
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(1000*f.tau_n)]);
    disp(ci)
end
struct_out = fits;
end
