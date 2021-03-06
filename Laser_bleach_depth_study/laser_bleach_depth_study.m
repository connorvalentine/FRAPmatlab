% FRAP laser bleach depth studies analysis
% Connor Valentine
%% to do list
% normalize intensity by the pre-bleach frame entirely??
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
    set(groot,'DefaultLineLineWidth',2)
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
    folder1 = 'z_Laser_Bleach_Study';
    folder2 = '3s timepoints trial 1 no ND';

% Part 2: Would you like to save the plots generated? 
% Note:         Select 'y' or 'n'
save_im = 'n'; 

%% %%%%%%%%%%%%%%%% Parameters Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Data Folder, structures, and sample ID data
global boxfolder plotfolder datafolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    datafolder = fullfile(boxfolder,folder1,folder2); % full path where outputs structures will be saved
    plotfolder = fullfile(boxfolder,folder1,folder2,'z_plots'); % plot folder within the data folder for the defined experiment
    
%initialize the structures to store all data in 
pb = struct();
data = struct();
%% reading in prebleach image
folder3 = 'prebleach';
% Make a list of each image name in our position folder
list = dir(fullfile(boxfolder,folder1,folder2,folder3,'Default','*.tif')); % lists all files with .tif ending in the position folder
data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
f1 = fullfile(data(1).folder,data(1).name);
im = imread(f1);
imsmooth = imgaussfilt(im,8); % now smooth image

pb.('imbds') = [min(im(:)),mean2(im)]; % parameters for plotting images
pb.('imsmoot') = imsmooth;
pb.('im') = im;
fig = figure('name','prebleach','visible','on');
    set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
    imlb = pb.imbds(1);
    imub = pb.imbds(2);
    % line profile across the radius
    x = [0 size(im,1)];
    y = [size(im,2)/2 size(im,2)/2];
    
    subplot(1,2,1);
    %show the first image after photobleaching
        hold on 
        imshow(im,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
        plot(x,y)
        hold off 
    subplot(1,2,2)
        c = improfile(im,x,y); % line profiles 
        c_smooth = improfile(imsmooth,x,y); % line profiles 
        pb.('profile') = c;
        pb.('sprofile') = c_smooth;
        hold on 
        plot(c(:,1,1),'r')
        plot(c_smooth(:,1,1),'b')
        % add radius and center
        axis([0 2048 1000 7500])
%% reading in each image to find the circle and make a line profile 
tic; % timing
% define other parameters 
folder3 = 'focused laser, no ND on Laser, 3s laser timepoints with intermediate pictures_2';
lim1 = 50; % minimum pixel radius to look for   

% Make a list of each image name in our position folder
list = dir(fullfile(boxfolder,folder1,folder2,folder3,'Default','*.tif')); % lists all files with .tif ending in the position folder
data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure

% Load first bleached frame
for i = 1:length(data)
    f1 = fullfile(data(i).folder,data(i).name);
    im = imread(f1);
    imsmooth = imgaussfilt(im,3.5); % sigma = 2 works for peak finding very well, but the center is all over the place


    % First we do some pre-analysis on the images 
    T = graythresh(imsmooth); 
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
    % now we sort through the objects to make sure the bleached circle is found
    if height(stats) == 0 % if no objects meet our criteria
        disp(['no good droplet found at ', i]) 
        % remove field name from id list
        id = rmfield(id,i);
        center = size(image)/2; % center is the middle of the image
        radius = 10; % radius set to 10 pixels (too small to be a real object)
    else 
        % the main object is the most solid one usually
        main_idx = find(stats.Solidity == max(stats.Solidity));
    end

    % Add data to fits structure if it doesn't meet the other errors
    if stats.EquivDiameter(main_idx)/2 < lim1
        % object is too small to be real 
        disp(['no good droplet found at ', num2str(i)]) 
        data(i).('center') = size(image)/2; % center is the middle of the image
        data(i).('radius') = 10; % radius set to 10 pixels (too small to be a real object)
        data(i).('imbds') = [min(im(:)),mean2(im)]; % parameters for plotting images
        data(i).('pixel_size') = 0.65;
    else % the main object chosen is fine. save the circle and reference regions into fits 
        data(i).('imbds') = [min(im(:)),mean2(im)]; % parameters for plotting images
        data(i).('radius') = (stats.EquivDiameter(main_idx)/2);
        data(i).('center') = stats.Centroid(main_idx,:); 
        % reference region dx and dy defined above, convert to pixel # list instead of coords
        idx = stats.PixelList(main_idx);
        idx = cell2mat(idx);
        data(i).('idx_bleach') = sub2ind(size(im),idx(:,2),idx(:,1));
        data(i).('pixel_size') = 0.65;
    end

    % adding time to the dataset
    % each timepoint represenents 3 seconds of bleaching
    data(i).('t') = i*3;
    
    % image profile for comparison to center finding
    center = data(i).center;
    x = [0 size(im,2)];
    y = [center(2) center(2)];
    c = improfile(im,x,y); % line profiles 
    c_smooth = improfile(imsmooth,x,y); % line profiles 
    
    ref_I = pb.sprofile(500);
%     data(i).('norm_factor') = mean(pb.sprofile(2:200))./mean(c_smooth(2:200));
    data(i).('norm_factor') = mean(pb.sprofile(2:200)./c_smooth(2:200));
    ind = find(c_smooth == min(c_smooth));
    prebleachint = mean(pb.sprofile(ind));
    data(i).('cs') = c_smooth; %cs = smoothed profile
    data(i).('cns') = c_smooth(:,1,1).*data(i).norm_factor;
    data(i).('bleachmin') = (min(c_smooth)./prebleachint)*data(i).norm_factor;
    
    x = linspace(1,length(data(i).cns),length(data(i).cns));
    y = 1-data(i).cns./pb.sprofile;
    [pk,lk,w,p] = findpeaks(y,x,'MinPeakProminence',0.1,'WidthReference','halfprom','MinPeakDistance',100);
    [pk2,lk2,w2,p2] = findpeaks(y,x,'MinPeakProminence',0.1,'WidthReference','halfheight','MinPeakDistance',100);
    if isempty(pk)
        disp(['empty peak find ', i])
        data(i).('rpk') = 0;
        data(i).('pk') = 0;
        data(i).('lk') = 0;
        data(i).('rpkh') = 0;
        data(i).('pkh') = 0;
        data(i).('lkh') = 0;
        pk =0; lk = 0; w = 0;
        pk2 = 0; lk2 = 0; w2 = 0;
    else
        data(i).('rpk') = w;
        data(i).('pk') = 1-pk;
        data(i).('lk') = lk;
        data(i).('rpkh') = w2;
        data(i).('pkh') = 1-pk2;
        data(i).('lkh') = lk2;
    end
    %plot figure
    fig = figure('name',num2str(data(i).t),'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
        imlb = data(i).imbds(1);
        imub = data(i).imbds(2);
        radius = data(i).radius;

        subplot(2,2,1);
        %show the first image after photobleaching
            hold on 
            imshow(im,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            hold off 
        subplot(2,2,2)
        %show the laser region and reference region as black circles
            hold on
            % first plot the black and white image
            imshow(bw1,'Border','tight','InitialMagnification', 'fit')
            h= viscircles(center,radius,'linewidth',0.2,'color','g');
            hold off 
        subplot(2,2,3)
        % line profile across the radius
            hold on 
            % plot(c(:,1,1),'r') % line profile (not smoothed)
            plot(data(i).cns,'b')
            plot([pb.sprofile],'g')
            % add radius and center
            plot([center(1)-radius,center(1)+radius],[1.2*mean([imub,imlb]),1.2*mean([imub,imlb])])
            axis([0 2048 3500 7500])
       subplot(2,2,4)
            hold on
            plot(data(i).lk,data(i).pk,'o','MarkerSize',12)
            errorbar(lk,1-(pk./2),w/2,'horizontal','color','k') % half prom
            errorbar(lk,1.05*(1-(pk./2)),radius,'horizontal','color','r') %imfind technique half height
            plot(data(i).cns./pb.sprofile,'b')
            axis([0 2048 0.4 1.1])
end
% save images if selected above
if save_im == 'y'
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 6];             % define location to save the images 
    a = fieldnames(id);
    pp = a{1};
    struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
    plot_name = [struct_name,'_',num2str(round(id.(i).plwt)),'wtp','_frames'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png  
else    
end
time2calculate = toc
%% plotting the results
fig = figure('name','radius','visible','on');
    set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
    subplot(2,2,1)
        x = [data.t];
        y = [data.radius];
        y2 = [data.rpk];
        y3 = [data.rpkh];
        hold on
        plot(x,y,'kd')
        plot(x,y2,'gd')
        plot(x,y3,'rd')
        xlabel('bleach time')
        ylabel('radius')
        legend('im analysis','peak find','peak find height','location','southeast')
    subplot(2,2,2)
        x = [data.t];
        y = [data.bleachmin];
        plot(x,y,'rd')
        xlabel('bleach time')
        ylabel('min bleach intensity')
    subplot(2,2,3)
        x = [data.t];
        y = [data.norm_factor];
        plot(x,y,'gd')
        xlabel('bleach time')
        ylabel('normalization factor intensity')  
    subplot(2,2,4)
        x = [data.t];
        y = [data.lk];
        plot(x,y,'gd')
        xlabel('bleach time')
        ylabel('center location by peakfinder')  
%% testing environments 
    