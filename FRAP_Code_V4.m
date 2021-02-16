% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - add gui to flag circles that are not right
% - recomment/clean code
% - save thge D and c vectors for plotting into a new structure to export
%%%%%%%%%%% primary
% - multiple reference regions to help ignore irregularities
% - why are Ifit values so low?
% - error analysis of r_i vs tau
% - Ri should be the same more or less for every one- flag for this error
% - error propagation of std_dev of radius?
% - why is ifit taking so long
%%%%%%%%%%% thoughts
% - uniformity metrics for the bleached spot and the reference region
% - data cleaning of outliers 
% - reference sample in each array? HDFL gel?
% - recomment/clean code

%%%%%%%%%%% Diffusivity fitting %%%%%%%%%%%%%%%%%%%%%5
% - make fitting go faster
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
    
%% choose the folders
% prebleach folder must be named 'prebleach'
% frap folder must be named 'frap'

% choose the folders
global folder1 folder2
    folder1 = 'P123_BSA_25C';
    folder2 = 'trial_1';

% save images? 
    save_im = 'y'; 
    % save_im = 'y';
% parameters that may change from experiment to experiment during troubleshooting    
    dy = 350;   
    dx = 350;   
    lim1 = 50; % minimum pixel radius to look for
% troubleshooting? 
    trouble = 'n'; % does not fit the data 
    % trouble = 'n'; % working as normal
% frame to analyze as first bleached frame (frame 1)
    t1 = 2;
    
%% Initialize Data Folder, structures, and sample ID data
    global boxfolder outputfolder plotfolder datafolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    datafolder = fullfile(boxfolder,folder1,folder2); % full path where struct will be saved
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,folder1,folder2,'z_plots');
    
%initialize the structure to store all data in 
alldata = struct();
fits = struct();
images = struct();

% load in the info structure.
id_name = [folder1,'_',folder2,'_','info','.mat'];
temp_struct = load(fullfile(datafolder,id_name));
id = temp_struct.id;

%% Reading in data V2 - First Frap image to find bleached circle
    tic; % timing
    close all
for field = fieldnames(id)'
    position = field{1};
%     position = 'pos10'
    % Make a list of each image name in our position folder
    folder3 = 'frap';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);

    % define metadata filename for the current position
    global metadata_filename 
    txt = dir(fullfile(data(1).folder,'*.txt'));% lists all files with .txt ending
    metadata_filename = fullfile(txt.folder,txt.name);
    
     % Analyze first bleached frame
    f1 = fullfile(data(t1).folder,data(t1).name);
    im = imread(f1);

    % now smooth image
    imsmooth = imgaussfilt(im,8);
    
    % for plotting images
    imlb = min(im(:));
    imub = mean2(im);
    T = graythresh(imsmooth); 

    % First we do some pre-analysis on the images to find the bleached circle
    % bw1 = imbinarize(im,T); % make image black and white
    bw1 = im2bw(imsmooth,T); % old way to do it seems to work better
    bwf = imcomplement(bw1);% flip black and white    
    bw2 = bwareaopen(bwf,300); % remove specks and dots smaller than lim2 pixels in area
    bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area
    % filling in the boundaries in black and white image to make solid shapes
    [B,L] = bwboundaries(bw3,'noholes');

    % These objects are now analyzed by regionprops(). 
    % The output is a table, called stats, that has the information for each object
    stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','PixelList'); 
    
    if height(stats) == 0 % if no objects meet our criteria
        disp(['no good droplet found at ', position]) 
        % remove field name from id list
        id = rmfield(id,position);
        center = size(image)/2; % center is the middle of the image
        radius = 10; % radius set to 10 pixels (too small to be a real object)
        
    else % if the drop_info table is not empty
        % This is the center and radii of every "good" object found. The
        % code works even if multiple "good" drops are found. 
        main_idx = find(stats.Solidity == max(stats.Solidity));
    end
    
    % other errors
    if stats.EquivDiameter(main_idx)/2 < lim1
        disp(['no good droplet found at ', position]) 
        id = rmfield(id,position);
        center = size(image)/2; % center is the middle of the image
        radius = 10; % radius set to 10 pixels (too small to be a real object)
    else
        circularity = stats.Circularity(main_idx);
        idx = stats.PixelList(main_idx);
        idx = cell2mat(idx);
        idx_bleach = sub2ind(size(im),idx(:,2),idx(:,1));
        center = stats.Centroid(main_idx,:); 
        radius = (stats.EquivDiameter(main_idx)/2);
        idx_ref = sub2ind(size(im),idx(:,2)+dy,idx(:,1)+dx);  % reference region dx and dy defined above, convert to pixel # list instead of coords
        
        % save the circle and reference regions into fits 
        fits.(position).('radius') = radius;
        fits.(position).('center') = center;
        fits.(position).('idx_ref') = idx_ref;
        fits.(position).('idx_bleach') = idx_bleach;
        fits.(position).('imbds') = [imlb,imub]; 
        fits.(position).('images') = images;
    end
           
        fig = figure('name',position,'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
        
        subplot(2,2,1);
        %show the first image after photobleaching
            hold on 
            imshow(im,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            hold off 
        subplot(2,2,2)
        %show the laser region and reference region as black circles
            % modify bw3 to show reference region 
            bw4 = bw1;
            bw4(idx_ref) = 0;
            hold on
            % first plot the black and white image
            imshow(bw4,'Border','tight','InitialMagnification', 'fit')
            h= viscircles(center,radius,'linewidth',0.2,'color','g');
            hold off 
        subplot(2,2,3)
        % line profile across the radius
            x = [0 size(im,2)];
            y = [center(2) center(2)];
            c = improfile(im,x,y); % line profiles 
            c1 = improfile(imsmooth,x,y); % line profiles 
            hold on 
            plot(c(:,1,1),'r')
            plot(c1(:,1,1),'b')
            plot([0,2024],[imub,imub]) 
            % add radius and center
            plot([center(1)-radius,center(1)+radius],[1.2*mean([imub,imlb]),1.2*mean([imub,imlb])])
            axis([0 2024 imlb 1.2*imub])
end  
disp("bleached circles found in " + string(round(toc)) +" s");

%% save images if selected above
if save_im == 'y'
    for field = fieldnames(id)'
        position = field{1};
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 6];             % define location to save the images 
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_frames'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png 
    end    
else    
end

%% Readin the pre-bleach images next
for field = fieldnames(id)'
    position = field{1};
    % Make a list of each image name in our position folder
    folder3 = 'prebleach';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    % Analyze the pre-bleach image
    f1 = fullfile(data(1).folder,data(1).name);
    im0=double(imread(f1));    
    idx_ref = fits.(position).idx_ref;
    idx_bleach = fits.(position).idx_bleach;
    radius = fits.(position).radius;
    center = fits.(position).center; 
    
    % intensity of the reference region before bleaching
    refI_t0 = mean(im0(idx_ref));
    % intensity of the spot before bleaching 
    I_t0 = mean(im0(idx_bleach));
    
    % add info to fits.
    fits.(position).('refI_t0') = refI_t0;
    fits.(position).('I_t0') = I_t0;
end  
%% Reading in the rest of the frap data 
   tic % timing
for field = fieldnames(id)'
    position = field{1};
 
    % Make a list of each image name in our position folder
    folder3 = 'frap';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    
    % define metadata filename for the current position
    global metadata_filename 
    txt = dir(fullfile(data(1).folder,'*.txt'));% lists all files with .txt ending
    metadata_filename = fullfile(txt.folder,txt.name);
    
    % pull out the prebleach image info
    refI_t0 = fits.(position).refI_t0;
    I_t0 = fits.(position).I_t0;
    idx_ref = fits.(position).idx_ref;
    idx_bleach = fits.(position).idx_bleach;

    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time_V4(data,t); 
        fits.(position).('pixel_size') = id.(position).pixel_size_manual;
        
%         if data(1).pixel_size == 0
%             fits.(position).('pixel_size') = 0.323; %for 20x in case it was not found 
%             fits.(position).('pixel_size') = 0.645; % 10x case 
%         else
%             fits.(position).('pixel_size') = data(1).pixel_size;
%         end

        % add pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        im = double(imread(f1));
        I_ti = mean(im(idx_bleach))./I_t0; % normalized by initial Intensity
        refI_ti = mean(im(idx_ref));
        if refI_ti ==0
            refI_ti = 0.001;
        else
        end
        normalized_i = (I_ti).*(refI_t0/refI_ti);
       % normalized_i = (I_ti./I_t0);

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('refI_ti') = refI_ti;
        data(t).('ref_ratio') = (refI_t0/refI_ti);
        data(t).('IN_ti') = normalized_i;
        alldata.(position) = data;
    end
    
    % add end frame to plots 
    radius = fits.(position).radius;
    center = fits.(position).center;
    imbds = fits.(position).imbds;
    disp([position,' loaded'])
end
a = round(toc);
disp(['Data loaded in ',num2str(a),' s.'])


%% Adding some additional things to the data structure
% 1. elapsed time from first "initial" frame
% 2. fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam

for field = fieldnames(alldata)'
    tic % timing start
    position = field{1}; % use{} bcuz field is a cell array
    % elapsed time from received time of first frap frame.
    time1 = alldata.(position)(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(alldata.(position))
        t = alldata.(position)(i).r_time;
        t = datevec(t);
        dt = etime(t,time1); % calculate elapsed time from the first image in seconds
        alldata.(position)(i).('dt') = dt; %save to data structure
    end
    
    if trouble == 'n' % ignores the fit process if trouble == 'y'
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
        options.StartPoint = [0.75 1 1000];
        options.Lower = [0.15,0.1,1];
        options.Upper = [1,10,10000];
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
    % output values from fit 
    b = round(toc); % timing
    disp([position,' data fit in ',num2str(b),' s.'])
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(round(f.tau))]);
    else
    end
    % timing 
    b = round(toc);
end
disp('time added')

%% Exporting the data structures for future use.
% can also save individual fields if we want to get crazy
% first specify output folder (where we want to save the data
% maybe would be good to do conc instead of pos# here? 

% fix this section to just use folder1 and folder2 for naming
% then can also get temperature to be a number 
cd(outputfolder);
struct_name = [folder1,'_',folder2]; % specify file name prefix
alldata_name = [struct_name,'_','data','.mat'];
fits_name = [struct_name,'_','fits','.mat'];
id_name = [struct_name,'_','id','.mat']; % fwd id struct along with it.
save(alldata_name,'alldata'); 
save(fits_name,'fits');
save(id_name,'id');
cd(mainfolder);
%% plotting the analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    c_temp = id.(position).plwt;
    c_temp = [num2str(c_temp), ' wt%'];
    fig = figure('name',c_temp,'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    
    % plotting data from image intensities, without the fits 
    t = [alldata.(position).dt]';
    t = t(t1:end);
    I = [alldata.(position).I_ti]';
    I = I(t1:end);
    IN = [alldata.(position).IN_ti]';
    IN = IN(t1:end);
    ref_ratio = [alldata.(position).ref_ratio]';
    ref_ratio = ref_ratio(t1:end);

    % plotting the data
    hold on
    plot(t,IN,'x','color',C(40,:))

    % if not troubleshooting, plot the fits as well
    if trouble == 'n'
    t_fit = [fits.(position).t_fit]';
    I_fit = [fits.(position).I_fit]';
    plot(t_fit,I_fit,'--')
    xc = fits.(position).xc;
    shapec = fits.(position).shapec;
    fill(xc',shapec,'r','FaceAlpha',0.1,'LineStyle','--','linewidth',1);
    else 
    end
    title(position);
    axis([0,max(t),0,1.2]);
    xlabel('Time [s]');
    ylabel('Normalized Intensity');
end

% save images if selected above
if save_im == 'y'
    for field = fieldnames(id)'
        position = field{1};
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 6];             % define location to save the images 
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_intensity_fit'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png 
    end    
else    
end

%% diffusivity data vs conc
if trouble == 'n'
fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 

% need for loop to unpack
i = 0;
for field = fieldnames(alldata)'
    position = field{1};
    i = i+1;
    D = fits.(position).D/55.1;
    c = id.(position).plwt;
    pxl = fits.(position).pixel_size;
    r = fits.(position).radius*pxl;
    rall(i) =r;
    conc(i) = c;
    Dall(i) = D; % from cheng to normalize by free soln BSA D
    Dneg = fits.(position).errD(1)/55.1;
    Dpos = fits.(position).errD(2)/55.1;
    Dnegall(i) = Dneg;
    Dposall(i) = Dpos;
end

subplot(1,2,1)
    hold on 
    errorbar(conc,Dall,Dnegall,Dposall,'db')
    axis([0.9*min(conc) 1.1*max(conc) 0.7*min(Dall) 2*max(Dall)])
    set(gca,'YScale','log');
    xlabel([id.(position).plur,' wt%'])
    ylabel('D/D_0 [\mum^2s^{-1}]')

subplot(1,2,2)
    hold on
    plot(conc,rall,'or')
    axis([0.9*min(conc) 1.1*max(conc) 0.8*min(rall) 1.2*max(rall)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('laser radius [um]')
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