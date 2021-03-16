% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - look up what the laser focus knob is actually doing
% - work on the master step3 file so can compare data faster

%%%%%%%%%%% primary
% - is photobleaching a permanent phenomena? irreversible?
% - why is radius so different than cheng??
% - cap fm off at one for the fits?

%%%%%%%%%%% thoughts
% - uniformity metrics for the bleached spot and the reference region
% - reference sample in each array? HDFL gel?

%%%%%%%%%%% Experiments to run %%%%%%%%%%%%%%%%%%%%%
% - energy balance on laser to make sure it isnt heating up

%%%%%%%%%%% Diffusivity fitting %%%%%%%%%%%%%%%%%%%%%5
% - other models for FRAP and normalization?
% - confirm that the diffusivity is calculated correctly
% - apply Ionic hopping mechanism model to frap data to see if that is a
% good metric?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%% Part 1: Initialize some basic parameters and defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Arial')
    set(0,'DefaultAxesFontSize',14)
% Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial','DefaultTextFontSize',14)
% default axis settings for better plots
    set(groot,'defaultaxeslinewidth',1)
    set(groot,'DefaultLineLineWidth',1)
    set(groot, 'DefaultAxesBox', 'on')
    set(groot, 'DefaultAxesYGrid', 'off')
    set(groot, 'DefaultAxesXGrid', 'off') 
% turn the error message beeps off
    beep off
% Reset the environment
    clear all;
    clear global all;
    clc
    close all;
% Establish the main directory, 
    global mainfolder 
    mainfolder = cd; 
    
%% %%%%%%%%%%%%%%%% Inputs Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything that must be selected when code is running smoothly 
% Part 1: choose the folders
% NOTE:         prebleach folder must be named 'prebleach'
%               frap folder must be named 'frap'

global folder1 folder2
    folder1 = 'P123_BSA_45C';
    folder2 = 'trial_1';

% Part 2: Would you like to save the plots generated? 
% Note:         Select 'y' or 'n'
save_im = 'y'; 

% Part 3: Are you troubleshooting the fits? 
% Note:         Select 'y' or 'n'
trouble = 'n'; 

%% %%%%%%%%%%%%%%%% Parameters Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Data Folder, structures, and sample ID data
global boxfolder outputfolder plotfolder datafolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    datafolder = fullfile(boxfolder,folder1,folder2); % full path where outputs structures will be saved
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,folder1,folder2,'z_plots'); % plot folder within the data folder for the defined experiment
    
%initialize the structures to store all data in 
alldata = struct();
fits = struct();
images = struct();

% load in the info structure. This strucutre is made by FRAP_data_cleaner
id_name = [folder1,'_',folder2,'_','info','.mat'];
temp_struct = load(fullfile(datafolder,id_name));
id = temp_struct.id;

%% Reading in pre-bleach image 
[fits] = fun_prebleach_frame(id,fits);
%% Reading in First Frap image to find bleached circle
close all
tic; % timing
% [id,fits,fig] = fun_first_bleached_frame(id,fits,save_im);
% define the global variables
global boxfolder folder1 folder2 plotfolder

% define other parameters 
folder3 = 'frap';
lim1 = 50; % minimum pixel radius to look for   
t1 = 2;% frame to analyze as first bleached frame. t1 = 2 because frame 1 is pic of laser bleaching

% parameters that may change from experiment to experiment during troubleshooting    
dy = 350;  % pixel dy shift to make reference ROI 
dx = 350;  % pixel dx shift to make reference ROI 
for field = fieldnames(id)' % iterate through the position list in id structure
    position = field{1};
    plotflag = 'y';

    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    
    % Load the pre-bleach frame
    im0 = fits.(position).pbim;
    im0sm = imgaussfilt(im0,2); % smoothed prebleach image
    
    % Load first bleached frame
    f1 = fullfile(data(t1).folder,data(t1).name);
    im1 = imread(f1);
    im1sm = imgaussfilt(im1,2); % now smooth image
    fits.(position).('imbds') = [min(im1(:)),mean2(im1)]; % parameters for plotting images
    
     % First we do some pre-analysis on the images 
    T = graythresh(im1sm); 
    % bw1 = imbinarize(im,T); % make image black and white
    bw1 = im2bw(im1sm,T); % old way to do it seems to work better
    bwf = imcomplement(bw1);% flip black and white    
    bw2 = bwareaopen(bwf,300); % remove specks and dots smaller than lim2 pixels in area
    bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area
    % filling in the boundaries in black and white image to make solid shapes
    [B,L] = bwboundaries(bw3,'noholes');
    
    % These objects are now analyzed by regionprops(). 
    % The output is a table, called stats, that has the information for each object
    stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','PixelList');
    % now we sort through the objects to make sure the bleached circle is found
    if height(stats) == 0 % if no objects meet our criteria
        disp(['no droplets found at ', position]) 
        % remove field name from id list
        % id = rmfield(id,position);
        temp_pos = fieldnames(id);
        temp_pos = temp_pos{1};
        disp(['using center of ',temp_pos, 'instead']);
        fits.(position).('center') = fits.(temp_pos).center; % center is the middle of the image
        fits.(position).('radius') = 10;% radius set to 10 pixels (too small to be a real object)
    else 
        % the main object is the most solid one usually
        main_idx = find(stats.Solidity == max(stats.Solidity));
    end
    
    % Add data to fits structure if it doesn't meet the other errors
    if stats.EquivDiameter(main_idx)/2 < lim1
        % object is too small to be real 
        disp(['droplet found, but not good at ', position]) 
        % id = rmfield(id,position);
        temp_pos = fieldnames(id);
        temp_pos = temp_pos{1};
        disp(['using center of ',temp_pos, 'instead']);
        fits.(position).('center') = fits.(temp_pos).center; % center is the middle of the image
        fits.(position).('radius') = 10; % radius set to 10 pixels (too small to be a real object)
    else % the main object chosen is fine. save the circle and reference regions into fits 
        fits.(position).('radius') = (stats.EquivDiameter(main_idx)/2);
        fits.(position).('center') = stats.Centroid(main_idx,:); 
    end
    fits.(position).('pixel_size') = id.(position).pixel_size_manual;
    
    % now we take some image profiles based on the center of the laser hole
    % First we make a vertical line profile across the radius
    center = fits.(position).center;
    radius = fits.(position).radius;
    profile_line_x = [center(1)-10 center(1)-10];
    profile_line_y = [0 size(im1,1)];
    
    % Now we record the profile on the first bleached image
    prof = improfile(im1sm,profile_line_x,profile_line_y); % line profiles 
    profile1 = prof(2:end);
    
    % Now we record the profile of the same line, but on the pre-bleached
    % image
    im0 = fits.(position).pbim;
    im0sm = imgaussfilt(im0,2); % now smooth image
    prof = improfile(im0sm,profile_line_x,profile_line_y); % line profiles 
    profile0 = prof(2:end);  
    % temporary normalization of profile 1 by profile 0 for the fit
    norm_ratio = mean(profile0(250:500))/mean(profile1(250:500));
    
    % now we normalize profile1 by profile0 and flip 
    y_norm = 1- (profile1./profile0)*norm_ratio;
    x_norm = linspace(1,2048,2048)'; 
    
    f = fit(x_norm,y_norm,'gauss1'); %number after gauss changes number of peaks in guass
    % a = height of peak
    % b = position of the center of the peak
    % c = sigma*sqrt(2)
    FWHM = 2*sqrt(log(2))*f.c1;
    peak_height = f.a1;
    peak_center = f.b1;
    
    ci = confint(f,0.95);
    FWHM_err = 2*sqrt(log(2))*abs(ci(1,1));
    
    % make a circle mask defined by the FWHM of the bleached area
    circleCenterX = center(1);
    circleCenterY = center(2); % might have these flipped
    circle_mask = false(2048,2048);
    [x_im, y_im] = meshgrid(1:2048,1:2048);
    circle_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (FWHM/2).^2) = true; 
    
    % now make a concentric circle mask to define the reference region
    r_out = 4;
    r_in = 2.5;
    reference_mask = false(2048,2048);
    [x_im, y_im] = meshgrid(1:2048,1:2048);
    reference_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (r_out*FWHM/2).^2 ...
        & (x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 >= (r_in*FWHM/2).^2) = true; 
    
    % now we mask the original image to only see pixels we want in the
    % prebleach image to find the intensity in the circle and the reference
    % region
    masked_image0 = im0; % Initialize with the entire image.
	masked_image0(~circle_mask) = 0; % Zero image outside the circle mask.
    masked_circle_indices = find(masked_image0); % find all pixels that are not == 0
    I_t0 = mean(masked_image0(masked_circle_indices));
    
    masked_reference_image0 = im0;
    masked_reference_image0(~ reference_mask) = 0;
    masked_reference_indices = find(masked_reference_image0); % find all pixels that are not == 0
    ref_0 = mean(masked_reference_image0(masked_reference_indices));
    
    % masking image 1 just for the plot
    masked_image1 = im1; % Initialize with the entire image.
	masked_image1(~circle_mask) = 0; % Zero image outside the circle mask.
    masked_reference_image1 = im1;
    masked_reference_image1(~ reference_mask) = 0;    
    
    % find the mean value of the pixels in this circle
    % add peak info to the fits 
    fits.(position).('FWHM') = FWHM;
    fits.(position).('dFWHM') = FWHM_err;
    fits.(position).('peak_height') = peak_height;
    fits.(position).('peak_center') = peak_center;
    fits.(position).('profile0') = profile0;
    fits.(position).('profile1') = profile1;
    fits.(position).('profile_line_x') = profile_line_x;
    fits.(position).('profile_line_y') = profile_line_y;
    fits.(position).('circle_mask') = circle_mask;
    fits.(position).('reference_mask') = reference_mask;
    fits.(position).('I_t0') = I_t0;
    fits.(position).('ref_0') = ref_0;
    
    fig = figure('name',position,'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
        imlb = fits.(position).imbds(1);
        imub = fits.(position).imbds(2);

        subplot(2,2,1);
        %show the first image after photobleaching
            hold on 
            imshow(im1,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            plot(profile_line_x,profile_line_y,'m')
            title('image1 with profile line')
            hold off 
        if plotflag == 'y'
        subplot(2,2,2)
        %show the laser region k circles
            hold on
            % first plot the black and white image
%             imshow(masked_image1,[imlb,imub],'Border','tight','InitialMagnification', 'fit')
            imshow(masked_reference_image1+masked_image1,[imlb,imub],'Border','tight','InitialMagnification', 'fit')
            h= viscircles(center,radius,'linewidth',0.2,'color','g');
            title('green imfind vs white improfile')
            hold off 
        subplot(2,2,3)
            hold on 
            plot(profile1,'b')
            plot(profile0,'g')
            title('prebleach profile vs im1 profile')
            axis([0 2024 2000  8000])
        subplot(2,2,4)
            hold on 
            plot(f,x_norm,y_norm)
            errorbar(peak_center,peak_height/2,FWHM/2,'horizontal')
            title('normalized im1 profile with guassian peak')
            axis([0 2024 0  1])
        else
        end
    % save images if selected above
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 6];             % define location to save the images 
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_frames'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png  
    else    
    end
end  
    
% designate the outputs properly
% struct_out1 = id;
% struct_out2 = fits;
% fig_out = fig;

disp("bleached circles found in " + string(round(toc)) +" s");

%% Reading in the rest of the frap data 
tic % timing
% [alldata] = fun_frap_frames(id,fits,alldata);
global boxfolder folder1 folder2
folder3 = 'frap';

for field = fieldnames(id)'
    position = field{1};
    
    % pull out the prebleach image info from fits structure
    
    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    profile_line_x = fits.(position).profile_line_x;
    profile_line_y = fits.(position).profile_line_y;
    ref_0 = fits.(position).ref_0;
    I_t0 = fits.(position).I_t0;
    circle_mask = fits.(position).circle_mask;
    reference_mask = fits.(position).reference_mask;

    % iterate through the images in data structure
    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time(data,t); 
   
        % calculate pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        imi = double(imread(f1));
        
        % find mean intensity inside of the circle
        masked_image_i = imi; % Initialize with the entire image.
        masked_image_i(~circle_mask) = 0; % Zero image outside the circle mask.
        masked_circle_indices = find(masked_image_i);
        I_ti = mean(masked_image_i(masked_circle_indices));
    
        % calculate mean intensity inside the reference region
        masked_reference_imagei = imi;
        masked_reference_imagei(~ reference_mask) = 0;
        masked_reference_indices = find(masked_reference_imagei);
        ref_i = mean(masked_reference_imagei(masked_reference_indices));
        norm_ratio = ref_0/ref_i;
        
        if ref_i ==0
            ref_i = 0.001;
            disp(position)
            disp('reference region mistake at frame #' +string(t))
        else
        end
        % normalized by prebleach Intensity(I_t0)
        % double normalized by the intensity ratio in the reference region
        % (refI_t0) divdided by the ref region ROI at time ti refI_ti
        normalized_i = (I_ti./I_t0).*(norm_ratio); 

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('ref_i') = ref_i;
        data(t).('ref_ratio') = norm_ratio;
        data(t).('IN_ti') = normalized_i;
    end
    
    % add elapsed time from received time of first frap frame to data struct.
    time1 = data(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(data)
        t = data(i).r_time;
        t = datevec(t);
        dt = etime(t,time1); % calculate elapsed time from the first image in seconds
        data(i).('dt') = dt; %save to data structure
    end
    
    % add the completed data structure to alldata
    alldata.(position) = data;
    
    disp([position,' loaded'])
end

disp(["Data loaded in " + string(round(toc)) + " s."])

%% Performing fits of the normalized data Cheng style
% close all
% fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam
tic;
t1 = 2; %ignores frame 1 from frap data which is the laser bleaching pic

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

    fm0 = (norm_i(end)-norm_i(1))/(1-norm_i(1));
    fm0_lb = 0.9*fm0;
    fm0_ub = 1.1*fm0;  
    if fm0_ub > 1
        fm0_ub = 1;
    else
    end

    % using this as a first guess point for the coefficients
    ft = fittype('fun_FRAPfit2(x,f,k,tau_n)');
        options = fitoptions(ft);
        options.StartPoint = [fm0,1.3, 3];
        options.DiffMinChange = 0.001;
        options.TolFun = 1e-8;
        options.Algorithm = 'Levenberg-Marquardt';
        options.Robust = 'off';  
    % perform the fit to the data
    [f,gof,output] = fit(x,norm_i, ft,options);
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
    dr = fits.(position).dFWHM .*0.5 .*fits.(position).pixel_size;
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
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(round(f.tau_n))]);
    disp(output)
end


disp(["all data fit in " + string(round(toc)) + " s."])
%
% adding variables to id that can be useful for plotting 
% need for loop to unpack
i = 0;
pd = struct(); %plot data structure
for field = fieldnames(id)'
    position = field{1};
    i = i+1;
    pd.('r')(i) = fits.(position).FWHM .* 0.5 .*fits.(position).pixel_size;
    pd.('c')(i) = id.(position).plwt;
    pd.('D')(i) = fits.(position).D/55.1; % from cheng to normalize by free soln BSA D
    pd.('Dneg')(i) = fits.(position).errD(1)/55.1;
    pd.('Dpos')(i) = fits.(position).errD(2)/55.1;
    pd.('fm0')(i) = fits.(position).fm0;   
    pd.('fm')(i) = fits.(position).fit_info.f;
    pd.('fmlb')(i) = fits.(position).fit_info.f - fits.(position).ci(1,1); 
    pd.('fmub')(i) =  fits.(position).ci(2,1) - fits.(position).fit_info.f;
    
end

% plotting the analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array

    c_temp = [num2str(id.(position).plwt), ' wt%'];
    fig = figure('name',c_temp,'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    
    % plotting data from image intensities, without the fits 
    t = [alldata.(position).dt]';
    t = t(2:end);
    I = [alldata.(position).I_ti]';
    I = I(2:end);
    IN = [alldata.(position).IN_ti]';
    IN = IN(2:end);
    ref_ratio = [alldata.(position).ref_ratio]';
    ref_ratio = ref_ratio(2:end);
    Ifitparam = fits.(position).('Ifitparam');

    % plotting the data
    hold on
%     plot(t,I,'x','color',C(40,:))
    plot(t,IN,'d','color',C(40,:),'markerfacecolor',C(40,:))
    % if not troubleshooting, plot the fits as well
    if trouble == 'n'
    t_fit = [fits.(position).t_fit]';
    I_fit = [fits.(position).I_fit]';
    plot(t_fit,I_fit,'--')
    xc = fits.(position).xc;
    shapec = fits.(position).shapec;
    fill(xc',shapec,'r','FaceAlpha',0.1,'LineStyle','--','linewidth',1);
    plot([0,t_fit(end)],[Ifitparam,Ifitparam])
    else 
    end
    title(position);
    axis([0,max(t),0,1.2]);
    xlabel('Time [s]');
    ylabel('Normalized Intensity');

% save images if selected above
if save_im == 'y'
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 6];             % define location to save the images 
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_intensity_fit'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png  
else    
end

end

% diffusivity data vs conc
%pd = plot data structure
if trouble == 'n'
fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 

subplot(1,3,1)
    set(gca,'YScale','log');
    hold all
    errorbar(pd.c,pd.D,pd.Dneg,pd.Dpos,'db')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.7*min(pd.D) 2*max(pd.D)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('D/D_0 [\mum^2s^{-1}]')

subplot(1,3,2)
    hold on
    plot(pd.c,pd.r,'or')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.8*min(pd.r) 1.2*max(pd.r)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('laser radius [um]')
    
subplot(1,3,3)
    hold on
    plot(pd.c,pd.fm,'om')
%     errorbar(pd.c,pd.fm,pd.fmlb,pd.fmub,'om')
    xlabel([id.(position).plur,' wt%'])
    ylabel('mobile fraction [fm]')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.5 1.2])
else
end

% save images if selected above
if save_im == 'y'
        a = fieldnames(id);
        posn = a{1};
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 10 4];             % define location to save the images 
        struct_name = [id.(posn).plur,'_',id.(posn).prot,'_',id.(posn).temp,'_',folder2];
        plot_name = [struct_name,'_DvC'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else    
end
%% Exporting the data structures for future use.
% can also save individual fields if we want to get crazy
struct_name = [folder1,'_',folder2]; % specify file name prefix
cd(outputfolder);% first cd to output folder (where we want to save the data)

alldata_name = [struct_name,'_','data','.mat'];% save alldata struct
save(alldata_name,'alldata'); 
fits_name = [struct_name,'_','fits','.mat']; % save fits struct
save(fits_name,'fits');
id_name = [struct_name,'_','id','.mat']; % save id struct
save(id_name,'id');
id_name = [struct_name,'_','pd','.mat']; % save pd struct (plot data)
save(id_name,'pd');

cd(mainfolder);% cd back to main folder
%% complete
disp('Analysis Complete')

% %% testing frapfit equation
% figure
% x = linspace(0,500,1000)';
% tau = 10;
% syms k
% y = 0.8*symsum((((-2).^k)./(1+k.*(1+ (2*x./tau))*factorial(k))), k, 0,12);
% y = double(y);
% plot(x,y)
% min(y)
% axis([0 500 0 1])
