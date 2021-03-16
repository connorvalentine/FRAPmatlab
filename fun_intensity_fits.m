%% Function description: 
% readins in the prebleach images and adds information to data structures
%% inputs and outputs
% inputs: id structure and fits structure. 

function [struct_out] = fun_intensity_fits(fits,alldata)
t1 = 2; %ignores frame 1 from frap data which is the laser bleaching pic

for field = fieldnames(alldata)'
    tic % timing start
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
    x = x(t1:end); % could have to adapt this later
    t_fit = linspace(0,x(end),250)';
    norm_i = [alldata.(position).IN_ti]';
    norm_i = norm_i(t1:end);
    
    % Set up the fit parameter search bounds and fit options
    % bounds on f are set by + or - 5% of the last intensity value now
%     fm0 = norm_i(end);
        % ind is from the Ifitparam calc - for some reason there is an initial drop in
        % intensity sometimes. This is due to poor circle radius choice I
        % believe
    
    % fm is the mobile fraction of proteins. 
    % in this case, it is taking as the ratio of the amount of intensity
    % recovered (from first frame after photobleacing aka t1 to tend)
    % borrom of ratio is the amount of intensity lost from the
    % photobleaching.
    
    fm0 = (norm_i(end)-norm_i(ind))/(1-norm_i(ind));
    fm0_lb = 0.95*fm0;
    fm0_ub = 1.05*fm0;
    
    if fm0_ub > 1
        fm0_ub = 1;
    else
    end
    disp(fm0)
    disp(fm0_lb)
    disp(fm0_ub)
    
    ft = fittype('fun_FRAPfit(x,f,k,tau)');
        options = fitoptions(ft);
        options.StartPoint = [fm0 1 1000];
        options.Lower = [fm0_lb,10,10];
        options.Upper = [fm0_ub,100,20000];
        options.DiffMinChange = 0.01;
        options.TolX = 0.001;
%         options.Robust = 'Bisquare';  
        options.Robust = 'off';  
    % perform the fit to the data
    [f,gof] = fit(x,norm_i, ft,options);
    I_fit = f(t_fit);
    
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
    ri = fits.(position).radius.*fits.(position).pixel_size;
    tau = f.tau;
    dtau = ci(:,3);
    D = (ri.^2)./(4*tau);
    Dlow = ((0.95*ri).^2)./(4*max(dtau)); % adds 5% variation to ri
    Dhigh = ((1.05*ri).^2)./(4*min(dtau)); % adds 5% variation to ri
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
    
     % output values from fit 
    b = round(toc); % timing
    disp([position,' data fit in ',num2str(b),' s.'])
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(round(f.tau))]);
end

disp('Diffusivity fits added')
struct_out = fits;

end
