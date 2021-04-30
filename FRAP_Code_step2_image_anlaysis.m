% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - good idea to save another matlab structure that eliminates certain
% positions from being plotted when importing into the analysis main filer
% - make movies that have the circle highlighted to show drift tracking is
% working 
%- set radius constant for intensity.
% - account for reference region breaching edges
% - comment and clean up entire code
% - create function that loads an image given the folders, and frame#, etc?
% - speed up the code somehow
% - find a way to just have an empty entry when no circle is found

% - add data structure to input if there was drift or no
% - look up what the laser focus knob is actually doing

%%%%%%%%%%% primary
% - is photobleaching a permanent phenomena? irreversible?
% - why is radius so different than cheng??
% - cap fm off at one for the fits?

%%%%%%%%%%% Experiments to run %%%%%%%%%%%%%%%%%%%%%
% - energy balance on laser to make sure it isnt heating up

%%%%%%%%%%% Diffusivity fitting %%%%%%%%%%%%%%%%%%%%%
% - other models for FRAP and normalization?
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
    folder1 = 'F87_BSA_55C';
    folder2 = 'trial_6';

% Part 2: Would you like to save the plots generated? 
% Note:         Select 'y' or 'n'
save_im = 'n'; 

% Part 3: Are you troubleshooting the fits? 
% Note:         Select 'y' or 'n'
trouble = 'n'; 

% do you want to make movies? (yes is first time analyzing)
global makemovies npoints
makemovies = 'n';

% number of points per wt%. is usually 3 now.
npoints = 3;
%% %%%%%%%%%%%%%%%% Parameters Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Data Folder, structures, and sample ID data
global boxfolder outputfolder plotfolder datafolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    datafolder = fullfile(boxfolder,folder1,folder2); % full path where outputs structures will be saved
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,'z_outputs','z_plots'); % plot folder within the data folder for the defined experiment
    
%initialize the structures to store all data in 
alldata = struct();
fits = struct();
images = struct();

% load in the info structure. This strucutre is made by FRAP_data_cleaner
id_name = [folder1,'_',folder2,'_','info','.mat'];
temp_struct = load(fullfile(datafolder,id_name));
id = temp_struct.id;

%% Reading in the prebleached image, and the first image after bleaching
% radius is found using improfile and fitting a gaussian curve to it
tic; % timing 

[id,fits,fig,stats] = fun_radius_finder_0(id,fits,save_im);

disp(["First Circles found in " + string(round(toc)) + " s."])
% [id,fits,fig,stats] = fun_radius_finder_i(id,fits,save_im,1000*[1.1733 0.9340],0.9,60);
%% Reading in the rest of the frap data 

% [alldata] = fun_frap_frames(id,fits,alldata);

% disp(["Data loaded in " + string(round(toc)) + " s."])
for field = fieldnames(id)' % iterate through the position list in id structure
    tic
    position = field{1};
    folder3 = 'frap';

    % make one figure that will be updated at every 10th timepoint just to
    % watch it run
    fig2 = figure('name',position,'visible','on');
    set(fig2, 'WindowStyle', 'Docked');  %figure will dock instead of free float
    
    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    
    % iterate through the images in data structure
    for t = 1:n_images 
        try
            % add time information to the data structure
            [data] = fun_time(data,t); 
           % add elapsed time from received time of first frap frame to data struct.
            time1 = data(1).r_time;
            time1 = datevec(time1);
            time_i = data(t).r_time;
            time_i  = datevec(time_i);
            dt = etime(time_i,time1); % calculate elapsed time from the first image in seconds
            data(t).('dt') = dt; %save to data structure

            % calculate pixel intensity information
            f1 = fullfile(data(t).folder,data(t).name);
            imi = double(imread(f1));
            if t == 1
                previous_center = fits.(position).center;
                previous_radius = fits.(position).radius;
                ref_0_i = fits.(position).ref_0;
                ref_i = fits.(position).ref_1;
                I_ti = fits.(position).I_t1;
            else 
                [circle_mask, reference_mask,center,ref_i,I_ti,ref_0_i] = fun_radius_finder_i(position,imi,fits,previous_center,previous_norm_ratio,t);
                previous_center = center;
            end
            
            if ref_i ==0
                ref_i = 0.001;
                disp(position)
                disp('reference region mistake at frame #' +string(t))
            else
                norm_ratio =  ref_0_i/ref_i;
                previous_norm_ratio = norm_ratio; % for feeding into next loop
            end
            
            % normalized by prebleach Intensity(I_t0)
            % double normalized by the intensity ratio in the reference region
            % (refI_t0) divdided by the ref region ROI at time ti refI_ti
            normalized_i = (I_ti./fits.(position).I_t0).*(norm_ratio); 

            % adding to data structure
            data(t).('I_ti') = I_ti;
            data(t).('ref_i') = ref_i;
            data(t).('ref_ratio') = norm_ratio;
            data(t).('IN_ti') = normalized_i;
        catch e
            warning('werid thing happened ')
            fprintf(1,'There was an error! The message was:\n%s',e.message);
            data(t).('I_ti') = data(t-1).I_ti;
            data(t).('ref_i') = data(t-1).ref_i;
            data(t).('ref_ratio') = data(t-1).ref_ratio;
            data(t).('IN_ti') = data(t-1).IN_ti;
        end
            
    end

    
    % add the completed data structure to alldata
    alldata.(position) = data;
    disp([position, 'loaded in ' + string(round(toc)) + " s."])
end
%% testing boundaries of image masking

    fields = fieldnames(id)';
    position = field{1};
    folder3 = 'frap';
    t = 95;

    % make one figure that will be updated at every 10th timepoint just to
    % watch it run
    fig = figure('name',position,'visible','on');
    set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
    
    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure

        try
            % calculate pixel intensity information
            f1 = fullfile(data(t).folder,data(t).name);
            imi = double(imread(f1));
 
            % make a circle mask defined by the FWHM of the bleached area
            r = 100;
            circleCenterX = 200;
            circleCenterY = 200; % might have these flipped
            circle_mask = false(2048,2048);
            [x_im, y_im] = meshgrid(1:2048,1:2048);
            circle_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (r/2).^2) = true; 
    
            % find mean intensity inside of the circle
            masked_image_i = imi; % Initialize with the entire image.
            masked_image_i(~circle_mask) = 0; % Zero image outside the circle mask.
            I_ti = mean(masked_image_i(masked_image_i > 0))
            
%             imshow(imi,[min(imi(:)),max(imi(:))])
            imshow(masked_image_i,[min(imi(:)),max(imi(:))],'Border','tight','InitialMagnification', 'fit')
        catch 
            warning('werid thing happened ')
        end

%% Performing fits of the normalized data Cheng style
tic 


    [fits] = fun_intensity_fits(fits,alldata);


disp(["all data fit in " + string(round(toc)) + " s."])
%%
% adding variables to id that can be useful for plotting 
% need for loop to unpack
i = 0;
pd = struct(); %plot data structure
for field = fieldnames(id)'
    position = field{1};
    i = i+1;
    if fits.(position).GoodFit == 'y'
        pd.('r')(i) = fits.(position).radius .*fits.(position).pixel_size;
        pd.('rlb')(i) = (fits.(position).radius - fits.(position).err_radius(1)).*fits.(position).pixel_size;
        pd.('rub')(i) = (fits.(position).err_radius(2) - fits.(position).radius).*fits.(position).pixel_size;
        pd.('c')(i) = id.(position).plwt;
        pd.('D')(i) = fits.(position).D/55.1; % from cheng to normalize by free soln BSA D
        pd.('Dneg')(i) = fits.(position).errD(1)/55.1;
        pd.('Dpos')(i) = fits.(position).errD(2)/55.1;
        pd.('fm0')(i) = fits.(position).fm0;   
        pd.('fm')(i) = fits.(position).fit_info.f;
        pd.('fmlb')(i) = fits.(position).fit_info.f - fits.(position).ci(1,1); 
        pd.('fmub')(i) =  fits.(position).ci(2,1) - fits.(position).fit_info.f;
    else
        a = [];
        pd.('r')(i) = 1;
        pd.('rlb')(i) = 1;
        pd.('rub')(i) = 1;
        pd.('c')(i) = id.(position).plwt;
        pd.('D')(i) = 1;
        pd.('Dneg')(i) = 1;
        pd.('Dpos')(i) = 1;
        pd.('fm0')(i) = 0;
        pd.('fm')(i) = 0;
        pd.('fmlb')(i) = NaN;
        pd.('fmub')(i) = NaN;
    end
end
%% 
% plotting the analyzed data 
C = jet;
counter = 0;
for field = fieldnames(alldata)'
    if counter == npoints
        counter = 0;
    else
    end
    counter = counter +1;
    
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
    if fits.(position).GoodFit == 'y'
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
        %% check if this wt% has been analyzed yet 
        
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_intensity_fit_',num2str(counter)];
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
    errorbar(pd.c,pd.r,pd.rlb,pd.rub,'or')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.8*min(pd.r) 1.2*max(pd.r)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('laser radius [um]')
    
subplot(1,3,3)
    hold on
%     plot(pd.c,pd.fm,'om')
    errorbar(pd.c,pd.fm,pd.fmlb,pd.fmub,'om')
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

