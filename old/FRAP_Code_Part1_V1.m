% FRAP analysis V1 
% Connor Valentine
%% to do list
% - have a more rigorous way of choosing a smaller radius for pixel analysis
% - remove speckles?
% - save plots into a results folder with names
% - make text file, or matlab file, with the starting inputs to run the
% script in each subfolder. for instance the file name, concentration in
% order, etc.
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
%     pluronic = 'F127';
%     folder1 = 'F127_BSA_25C';
%     folder2 ='trial_1';
%     conc = [17.5;20;22.5;25;27.5;30]; 
%     temperature = '25C'; % later could have this pickup from the filename?
%     dThresh = 0.7;
%     dy = 300;   
%     dx = 300;

% % F87 dataset
%     pluronic = 'F87';
%     folder1 = 'F87_BSA_25C';
%     folder2 ='trial_1';
%     conc = [25;30;35;37.5;40;42.5]; 
%     temperature = '25C'; % later could have this pickup from the filename?
%     dThresh = 0.9; %0.7 didnt work
%     dy = 350;   
%     dx = 350;
    
% P123 dataset trial 1
    pluronic = 'P123';
    folder1 = 'P123_BSA_25C';
    folder2 ='trial_1';
    conc = [20;25;27.5;30;32.5;35]; % later could have this pickup from a text file?
    temperature = '25C'; % later could have this pickup from the filename?
    dThresh = 1; %0.7 didnt work
    dy = 350;   
    dx = 350;
    
% frame to analyze as first bleached frame (frame 1)
    t1 = 2;
   
% folder where data is sorted data
    global boxfolder outputfolder plotfolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,folder1,folder2,'z_plots');
    
%initialize the structure to store all data in 
alldata = struct();
fits = struct();
id = struct();
images = struct();

% get the subfolder names from the folders instead of naming by number
    folder3 = 'frap';
    d = dir(fullfile(boxfolder,folder1,folder2,folder3));
% remove all files (isdir property is 0)
    dfolders = d([d(:).isdir]);
% remove '.' and '..'from the list
	dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    pos_struct = rmfield(dfolders,{'bytes','date','isdir','datenum','folder'});  % cleaning up the data structure
    npos = length(pos_struct); % number of positions
    npos = 1; % for testing
    
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

%% Reading in data V2 - First Frap image to find bleached circle
    tic; % timing
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
            
        % plotting
        fig = figure('name',position,'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
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
        
            % saving the figure as a high quality png 
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8 6];             % define location to save the images 
            a = fieldnames(id);
            pp = a{1};
            struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
            plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_frames'];
            plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
            print(fig,plot_path, '-painters', '-dpng', '-r600')

            % save the circle and reference regions into fits 
            fits.(position).('radius') = radius;
            fits.(position).('center') = center;
            fits.(position).('idx_ref') = idx_ref;
            fits.(position).('idx_bleach') = idx_bleach;
            fits.(position).('imbds') = [imlb,imub]; 
            fits.(position).('images') = images;
end  
disp("bleached circles found in " + string(round(toc)) +" s");

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
   
    % save the circle and reference regions into fits 
    fits.(position).('refI_t0') = refI_t0;
    fits.(position).('I_t0') = I_t0;
end  
%% Reading in the rest of the frap data 
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
    
    % pull out the prebleach image info
    refI_t0 = fits.(position).refI_t0;
    I_t0 = fits.(position).I_t0;
    idx_ref = fits.(position).idx_ref;
    idx_bleach = fits.(position).idx_bleach;

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
        I_ti = mean(im(idx_bleach));
        refI_ti = mean(im(idx_ref));
        normalized_i = (I_ti./I_t0).*refI_t0/(refI_ti);
        normalized_i = (I_ti./I_t0);

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('refI_ti') = refI_ti;
        data(t).('IN_ti') = normalized_i;
        alldata.(position) = data;
    end
    % add end frame to plots 
    radius = fits.(position).radius;
    center = fits.(position).center;
    imbds = fits.(position).imbds;
    a = round(toc);
    disp([position,' loaded in ',num2str(a),' s.'])
end
disp("Data Loaded and intensity measured")

%% Adding some additional things to the data structure
% 1. elapsed time from first "initial" frame
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    % The elapsed time for each frame is based on the first recieved
    % time (bleaching frame) in the data structure
    time1 = alldata.(position)(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(alldata.(position))
        t = alldata.(position)(i).r_time;
        t = datevec(t);
        dt = etime(t,time1); % calculate elapsed time from the first image in minutes
        % add to the data structure
        alldata.(position)(i).('dt') = dt;
    end
end

%% Exporting the data structures for future use.
% can also save individual fields if we want to get crazy
% first specify output folder (where we want to save the data
% maybe would be good to do conc instead of pos# here? 
cd(outputfolder);

% save the data structures for later plotting/analysis
a = fieldnames(id);
pp = a{1};
struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
alldata_name = [struct_name,'_','data','.mat'];
id_name = [struct_name,'_','id','.mat'];
save(alldata_name,'alldata');
save(id_name,'id');
cd(mainfolder);
%% plotting the un-analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    fig = figure('name',[position,'data'],'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    title(position)
    % choosing data to plot
    t = [alldata.(position).dt]';
    t = t(t1:end);
    I = [alldata.(position).IN_ti]';
    I = I(t1:end);
    I_raw = [alldata.(position).I_ti]';
    I_raw = I_raw(t1:end);
    % adding photobleaching effect
    ratioI_ti = [alldata.(position).refI_ti]'./[fits.(position).refI_t0]';
    ratioI_ti = ratioI_ti(t1:end);   
    % plotting and marker choice, linestyle, etc
    hold on
    plot(t,I,'x','color',C(40,:))
    plot(t,I_raw,'o','color',C(20,:))
    plot(t,ratioI_ti,'d')
    axis([0,max(t),0,max(I_raw)])
    xlabel('time');
end
