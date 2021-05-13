% fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam
function [struct_out] = fun_intensity_fits_nofm(fits,alldata)

for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    disp(position)
    try
        % adding fit parameters 
        global Ifitparam fm0
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
        fm_ind = length(norm_i)-4; % 
        
        % fm is the mobile fraction of proteins. 
        % in this case, it is taking as the ratio of the amount of intensity
        % recovered (from first frame after photobleacing aka t1 to tend)
        % borrom of ratio is the amount of intensity lost from the
        % photobleaching.
        
        fm0 = (mean(norm_i(fm_ind:end))-norm_i(1))/(1-norm_i(1)); % 
        fm0err = std(norm_i(fm_ind:end)); % 

        if fm0 > 1
            fm0 = 1;
        else
        end
        
        % remove any bad frames where normalized intensity is not a number
        NANind = find(isnan(norm_i));
        norm_i(NANind) = [];
        x(NANind) = [];

        % fake time data to put into the fit equation
        t_fit = linspace(0,x(end),500)'; 

        % fit weights
        weights = ones(length(norm_i),1);
        weights(1: 20) = 100;

        k0 = 1.2;
        tau0 = 10;
        ft  = fittype('fun_FRAPfit_nofm(x,k,tau_n)');
            options = fitoptions(ft);
            options.StartPoint = [k0, tau0];
            options.Lower =      [0.8, 1]; % 55C P123
            options.Upper =      [10,30]; % 55C
%             options.Lower =      [0.8,0.06]; % 55C F127 and F87
%             options.Upper =      [10,100]; % 55C
%            options.Lower =      [0.5,1]; % 45C
%             options.Upper =      [10,100];%  45C
            options.DiffMinChange = 0.001;
            options.TolFun = 1e-8;
            options.Robust = 'off';  
            options.Weights = weights;


        % perform the fit to the data
        [f,gof,output] = fit(x,norm_i, ft,options);

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
        ri = fits.(position).radius .*fits.(position).pixel_size;
        dr = (fits.(position).err_radius(2) - fits.(position).radius) .*fits.(position).pixel_size;
        tau = f.tau_n;
        dtau = abs(ci(:,end));
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
        fits.(position).('fm') = fm0;
        fits.(position).('fm0err') = fm0err;
        fits.(position).('GoodFit') = 'y';
         % output values from fit 
        disp(["f= " + fm0 + ",k= "+string(f.k)+",tau= "+string(f.tau_n)]);
        disp(ci)
    catch  e
        warning('fit failed for ')
        fprintf(1,'error message \n%s',e.message);
        fits.(position).('fit_info') = [];
        fits.(position).('ci') = [];
        fits.(position).('D') = [];
        fits.(position).('errD') = [];
        fits.(position).('gof') = [];
        fits.(position).('I_fit') = [];
        fits.(position).('t_fit') = t_fit;
        fits.(position).('xc') = [];
        fits.(position).('shapec') = []; 
        fits.(position).('Ifitparam') = Ifitparam;
        fits.(position).('fm0') = fm0;
        fits.(position).('fm0err') = fm0err;
        fits.(position).('GoodFit') = 'n';
    end
end
struct_out = fits;
end