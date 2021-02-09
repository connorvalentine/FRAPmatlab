%% add the intro

%% loading the data struture from load

%% fitting the data
% 2. fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam

    % fitting the data to the equation in Cheng paper for FRAP
    % use x for the variable instead of t
    x = [alldata.(position).dt]';
    x = x(t1:end); % could have to adapt this later
    t_fit = linspace(0,x(end),250)';
    norm_i = [alldata.(position).IN_ti]';
    norm_i = norm_i(t1:end);
    
    % adding fit parameters 
    global Ifit
        I_t1 = alldata.(position)(2).I_ti; % skip laser frame
        R_t1 = alldata.(position)(2).refI_ti; %% check
        I_t0 = fits.(position).I_t0;
        R_t0 = fits.(position).refI_t0;
        Ifit = (I_t1./I_t0).*(R_t0./R_t1);
        disp(["Ifit= " + string(Ifit)])
    
    % now setting fit options and fitting
    ft = fittype('FRAPfit(x,f,k,tau)');
        options = fitoptions(ft); 
        options.StartPoint = [0.75 1 10];
        options.Lower = [0.15,0,0];
        options.Upper = [1,10,5000];
        options.Robust = 'Bisquare';  
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
    Dl = (ri.^2)./(4*tau + max(dtau));
    Dr = (ri.^2)./(4*tau - min(dtau));
    Dneg = D-Dl;
    Dpos = Dr-D;
    errD = [Dneg Dpos];
    % creating new data structure for fit data
    fits.(position).('fit_info') = f;
    fits.(position).('D') = D;
    fits.(position).('errD') = errD;
    fits.(position).('gof') = gof;
    fits.(position).('I_fit') = I_fit;
    fits.(position).('t_fit') = t_fit;
    fits.(position).('xc') = xc;
    fits.(position).('shapec') = shapec;
    
    % check fm value
    I_inf = alldata.(position)(end).IN_ti;
    I_1 = alldata.(position)(1).IN_ti;
    fm_temp = (I_inf - I_1)./(I_t0 - I_1);
    fits.(position).('fm_calc') = fm_temp;
    % timing 
    b = round(toc);
    disp([position,' data fit in ',num2str(b),' s.'])
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(round(f.tau))]);
    
    %% export the fits structure for later use
    a = fieldnames(id);
pp = a{1};

struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
alldata_name = [struct_name,'_','data','.mat'];
fits_name = [struct_name,'_','fits','.mat'];
    save(fits_name,'fits');
    
%% plotting the fits
%% plotting the analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array

    fig = figure('name',position,'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    
%     subplot(2,1,1)
    title(position)
    % choosing data to plot
    t = [alldata.(position).dt]';
    t = t(t1:end);
    I = [alldata.(position).IN_ti]';
    I = I(t1:end);
    % fit info
    t_fit = [fits.(position).t_fit]';
    I_fit = [fits.(position).I_fit]';
    xc = fits.(position).xc;
    shapec = fits.(position).shapec;
    
    % adding photobleaching effect
    refI_ti = [alldata.(position).refI_ti]'./[fits.(position).refI_t0]';
    refI_ti = refI_ti(t1:end);
    
    % plotting and marker choice, linestyle, etc
    fill(xc',shapec,'r','FaceAlpha',0.1,'LineStyle','--','linewidth',1);
    hold on
    plot(t_fit,I_fit,'--')
    plot(t,I,'x','color',C(40,:))
    

%     plot(t,refI_ti,'d')
    axis([0,max(t),0,1])
    xlabel('time');
%     subplot(2,1,2);
%         center = fits.(position).center;
%         radius = fits.(position).radius;
%         hold on 
%         image = alldata.(position)(end).image;
%         imshow(image,[min(image(:)),mean2(image)],'Border','tight','InitialMagnification', 'fit');
%         viscircles(center,radius,'linewidth',0.2,'color','g');
%         hold off 
end
    
%% diffusivity data vs conc

fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 

% need for loop to unpack
hold on
i = 0;
for field = fieldnames(alldata)'
    position = field{1};
    i = i+1;
    D = fits.(position).D/55.1;
    c = id.(position).plwt;
    Dall(i) = D; % from cheng to normalize by free soln BSA D
    Dneg = fits.(position).errD(1)/55.1;
    Dpos = fits.(position).errD(2)/55.1;
    plot(c,D,'db')
    errorbar(c,D,Dneg,Dpos,'ob')
end
    axis([0.8*min(conc) 1.2*max(conc) 0.8*min(Dall) 1.2*max(Dall)])
    set(gca,'YScale','log');
    xlabel('F127 wt%')
    ylabel('D/D_0 [\mum^2s^{-1}]')
%% saving the plots