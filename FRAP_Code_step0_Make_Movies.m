% Making movies from the FRAP data
% Connor Valentine
%% to do list
% - have a more rigorous way of choosing a smaller radius for pixel analysis
% - remove speckles?
% - save plots into a results folder with names
% - move movie script to new matlab script
%%%%%%%%%%%%%%%%%future future future goals %%%%%%%%%%%%%%%%%%%%%%%%%
% - split this into two matlab scripts, one to load the pics and one to
% analyze 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - first script should output for each position: 
% - a data structure saved with:
% 1. Identifying information: what day this data was taken
% 2. temperature, concentration of pluronic, conc of BSA
% 2.1 necessary intensity data, etc for each frame
% 3. save them all to a master folder? could make fitting the data easier
% later, perhaps with filenames: Plur_conc_temp_protein?

% - second script for fitting the data
% 1. takes in the strctures from 1) and fits the data

% - make fitting go faster
% - confirm that the diffusivity is calculated correctly
% - add concentration too
% - apply Ionic hopping mechanism model to frap data to see if that is a
% good metric?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%% Part 1: Initialize some basic parameters and defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Arial')
    set(0,'DefaultAxesFontSize',9.5)
% Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial','DefaultTextFontSize',9.5)
% default axis settings for better plots
    set(groot,'defaultaxeslinewidth',1)
    set(groot,'DefaultLineLineWidth',1)
    set(groot, 'DefaultAxesBox', 'on')
    set(groot, 'DefaultAxesYGrid', 'off')
    set(groot, 'DefaultAxesXGrid', 'off')
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

% % F127 dataset
    pluronic = 'F127';
    folder1 = 'F127_BSA_25C';
    folder2 ='trial_1';
    conc = [17.5;20;22.5;25;27.5;30]; 
    temperature = '25C'; % later could have this pickup from the filename?
    dThresh = 0.7;
    dy = 300;   
    dx = 300;

% F87 dataset
%     pluronic = 'F87';
%     folder1 = 'F87_BSA_25C';
%     folder2 ='trial_1';
%     conc = [25;30;35;37.5;40;42.5]; 
%     temperature = '25C'; % later could have this pickup from the filename?
%     dThresh = 0.9; %0.7 didnt work
%     dy = 350;   
%     dx = 350;
    
% % P123 dataset trial 1
%     pluronic = 'P123';
%     folder1 = 'P123_BSA_25C';
%     folder2 ='trial_1';
%     conc = [20;25;27.5;30;32.5;35]; % later could have this pickup from a text file?
%     temperature = '25C'; % later could have this pickup from the filename?
%     dThresh = 1; %0.7 didnt work
%     dy = 350;   
%     dx = 350;
    
% number of positions 
% later this can be auto picked up by length of pos folders
    
% frame to analyze as first bleached frame (frame 1)
    t1 = 2;
  
% choose to make movies or not
  makemovies = 'y'; % change to 'y' for movies to be made
  
% folder where data is sorted data
    global boxfolder outputfolder plotfolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    outputfolder = fullfile(boxfolder,'z_outputs');
    
%initialize the structure to store all data in 
alldata = struct();
fits = struct();
id = struct();
images = struct();

% get the subfolder names from the folders instead of naming by number
    d = dir(fullfile(boxfolder,folder1,folder2,'frap'));
% remove all files (isdir property is 0)
    dfolders = d([d(:).isdir]);
% remove '.' and '..'from the list
	dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    pos_struct = rmfield(dfolders,{'bytes','date','isdir','datenum','folder'});  % cleaning up the data structure
    npos = length(pos_struct);
    
% re-order positions so they're arranged from lowest number to highest
% instead of alphabetical
positions = zeros(npos,1);
for p = 1:npos
    A = regexp(pos_struct(p).name,'\d*','Match');
    positions(p) = str2double(A{1});
end 
    %resort the position order
    [sortedpos,sortorder] = sort(positions);

% form position data structure with ID info 
for p = 1:npos
    % First we make the folder name for the position we want
    % add the positions in numerical order:
    idx = sortorder(p);
    position = pos_struct(idx).name;
    id.(position).('plur') = pluronic;
    id.(position).('plwt') = conc(p);
    id.(position).('prot') = 'BSA';
    id.(position).('protc') = '4p2';
    id.(position).('temp') = temperature;
end

%% Reading in data V2
close all
clc

% analyze frap data first.
    % timing
    tic;
    count = 1;
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
    
     % Analyze first bleached frame
    f1 = fullfile(data(t1).folder,data(t1).name);
    im = imread(f1);

    im1=double(im);
    % for plotting images
    imlb = min(im(:));
    imub = mean2(im);
    %threshold = (max(im1(:)) - min(im1(:)));
    T = graythresh(im)*dThresh; % works better with 0.7 for some

    % First we do some pre-analysis on the images to find the bleached circle
    bw1 = imbinarize(im,T); % make image black and white
    bwf = imcomplement(bw1);% flip black and white    
    bw2 = bwareaopen(bwf,200); % remove specks and dots smaller than lim2 pixels in area
    bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area
    % filling in the boundaries in black and white image to make solid shapes
    [B,L] = bwboundaries(bw3,'noholes');

    % These objects are now analyzed by regionprops(). 
    % The output is a table, called stats, that has the information for each object
    stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','PixelList'); 
    main_idx = find(stats.Solidity == max(stats.Solidity));

    circularity = stats.Circularity(main_idx);
    idx = stats.PixelList(main_idx);
    idx = cell2mat(idx);
    idx_bleach = sub2ind(size(im1),idx(:,2),idx(:,1));
    center = stats.Centroid(main_idx,:); 
    radius = (stats.EquivDiameter(main_idx)/2);

    % reference region dx and dy defined above 
    
    % convert to pixel # list instead of coords
    idx_ref = sub2ind(size(im1),idx(:,2)+dy,idx(:,1)+dx);  
        
       % %% show the first image after photobleaching
        % modify bw3 to show reference region 
        bw4 = bw1;
        bw4(idx_ref) = 0;
        
        %figure(f1);
        a = figure('name',position,'visible','on');
        set(a, 'WindowStyle', 'Docked');  %figure will dock instead of free float
        
        subplot(1,2,1);
            hold on 
            imshow(im,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            % Increase the size and position of the image 
            pos1 = get(gca,'Position'); pos1(1) = 0; pos1(2) = 0.2; pos1(3) = 0.5; 
            set(gca,'Position',pos1); 
            hold off 
        subplot(1,2,2)
            hold on
            % first plot the black and white image
            imshow(bw4,'Border','tight','InitialMagnification', 'fit')
            h= viscircles(center,radius,'linewidth',0.2,'color','g');
            % Increase the size and position of the image 
            pos2 = get(gca,'Position'); pos2(1)= 0.5; pos2(2) = 0.2; pos2(3) = 0.5;
            set(gca,'Position',pos2);
            hold off 
            
        % save the circle and reference regions into fits 
            fits.(position).('radius') = radius;
            fits.(position).('center') = center;
            fits.(position).('imbds') = [imlb,imub]; 
end  
disp("bleached circles found in " + string(round(toc)) +" s");
% now repeating for pre-bleached images


%% 
for field = fieldnames(id)'
    position = field{1};
    % timing
    tic
    % Make a list of each image name in our position folder
    folder3 = 'prebleach';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    % Analyze the pre-bleach image
    f1 = fullfile(data(1).folder,data(1).name);
    im0=double(imread(f1));
    if makemovies == 'y'
        images.(position).images{1} = im0;
    else
    end
end  
%%
for field = fieldnames(id)'
    position = field{1};
    tic % timing
    % Make a list of each image name in our position folder
    folder3 = 'frap';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    
    % define metadata filename for the current position
    global metadata_filename 
    txt = dir(fullfile(data(1).folder,'*.txt'));% lists all files with .txt ending
    metadata_filename = fullfile(txt.folder,txt.name);
    
    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time_V4(data,t); 

        if data(1).pixel_size == 0
            fits.(position).('pixel_size') = 0.323; %for 20x in case it was not found 
        else
            fits.(position).('pixel_size') = data(1).pixel_size;
        end
        % add pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        im = double(imread(f1));
        if t >1 && makemovies == 'y'
            images.(position).images{t} = im;
        else
        end
    end
    alldata.(position) = data;
    a = round(toc);
    disp([position,' loaded in ',num2str(a),' s.'])
end

disp("Data Loaded and intensity measured")

%% Adding some additional things to the data structure
% 1. elapsed time from first "initial" frame
% 2. fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam

for field = fieldnames(alldata)'
    tic
    position = field{1}; % use{} bcuz field is a cell array
    % 1. adding the elapsed time for each frame, based on the first recieved
    % time in the data structure. In this case, the first received time is the
    % time of the first frame in "initial" folder
    time1 = alldata.(position)(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(alldata.(position))
        t = alldata.(position)(i).r_time;
        t = datevec(t);
        dt = etime(t,time1)/60; % calculate elapsed time from the first image in minutes
        % add to the data structure
        alldata.(position)(i).('dt') = dt;
    end
end

%% making the videos for each position
% todo with videos: 
% -) change save location to a video folder?
% -) names should have ID info.
a = fieldnames(id);
pp = a{1};
struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
moviefolder = fullfile(outputfolder,'movies');
cd(moviefolder);
% 
if makemovies == 'y'
    for field = fieldnames(alldata)'
    position = field{1};

    disp(['Making video for ', position])
    movie_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp'];
    figure(14);
    set(gcf, 'Color','white')
    set(gca, 'nextplot','replacechildren', 'Visible','off');

    %# create AVI object
    nFrames =  length(images.(position).images);
    vidObj = VideoWriter([movie_name, '.avi']);
    vidObj.Quality = 100;
    vidObj.FrameRate = 10;
    open(vidObj);
    %# create movie

        for k=1:nFrames
        time = round(alldata.(position)(k).dt,1);
        timetext = [num2str(time),' min'];
        im = images.(position).images{k};
        imbds = fits.(position).imbds;
        imshow(im,imbds,'Border','tight','InitialMagnification', 'fit');
        hold on
        viscircles(center,radius,'linewidth',0.2,'color','g'); % circle for bleached reason
        viscircles(center+dy,radius,'linewidth',0.2,'color','r'); % circle for reference region
        text(1000,200,timetext,'color','m');
        % Increase the size and position of the image 
        writeVideo(vidObj, getframe(gca));
        hold off
        end

    close(gcf)

    %# save as AVI file, and open it using system video player
    close(vidObj);
    end
else
end
cd(mainfolder)