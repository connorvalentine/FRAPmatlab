ti%% micromanager tif folder analysis of the metadata
% Connor Valentine and Deyu Yang
% 2020 Walker Group

%% notes for things to add

% 3. button to flag a droplet result as bad so that it is not included in
%    the analysis

% 4. Centroid tracking between frames (connecting data between times)
% 5. Check to remind the user that pixel calibration might not be accurate
%    depends on the user making sure

% 9. Line scan of the drop? would have a profile across the drop

% 11. is it possible to go backwards in the analysis

% IDEA: what about a numbering system during troubleshooting? Because
% drop_info is sometimes more than one droplet long,we could have the user
% give the option to choose which drop it is

% 12. Right now if fun_droplets finds more than 1 droplet, it just chooses
% the first one in the list. Can we make this better? (see IDEA)

% now the problem with a numbering system, is how do we save the users
% droplet choices?

%% First section to reset the script, variables, and plot defaults
% Reset the environment
    clear all
    clc

% turn the error message beeps off
    beep off
    
% This just changes the default font, font size, and line width for better
% plots
    set(0,'DefaultAxesFontSize',16);
    set(findall(gcf,'type','text'),'FontSize',16, 'FontWeight','bold',...
    'fontName', 'arial');
    set(0,'defaultaxeslinewidth','default');
    set(0,'DefaultLineLineWidth',1);
    set(0,'DefaultAxesBox','on');
    set(0,'DefaultAxesYGrid','off');
    set(0,'DefaultAxesXGrid','off');

    close all;
%% inputs for the user
% device parameters
    h = 93;% device height in um
    C0 = 0.464; % initial solute concentration (in M)
    
% number of positions 
    n = 3;
    
% image analysis parameters
    level1 = 0.2; % black and white cutoff for the images
    level2 = 0.02;
    lim1 = 0.1; % cutoff for "circularity" parameter. Will ignore drops with circularity less than lim1
    lim2 = 0.1; % cutoff for "solidity" parameter
    lim3 = 100; % ignore "droplets" smaller than 100 pixels in diameter
    troubleshoot = 0; % make this = 1 if you want to output plots to click through to confirm circle finding

%% Before we do any analysis, specify what folders we are using
    % this line of code lists the filepath for the folder that this matlab
    % script is currently in. The user should not have to input the
    % filepath by hand
    mainfolder = cd; 
    subfolder1 = 'initial';
    subfolder2 = 'dehydration';
    
%% Analysis - first read in initial frames, and then dehydration frames
%initialize the structure to store all data in 
alldata = struct();

% iterate through each folder (initial, dehydration, etc)
for j=1:2 
    if j == 1
        subfolder = subfolder1;
        level = level1 ;
        n_images = 2;
    else
        subfolder = subfolder2;
        level = level2;  
        n_images = 10; %# of images
    end

    %iterating over each position and time 
    for p = 1:n
        % First we make the folder name for the position we want
        prefix = 'Pos';
        position = [prefix,num2str(p-1)];

        % Make a list of each image name in our position folder
        list = dir(fullfile(mainfolder,subfolder,position,'*.tif')); % lists all files with .tif ending in the position folder
        data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
        % if u want to analyze all images
        %       n_images = length(data); %# of images

        % Information about the metadata file
        % global variable means we can use it inside of other functions
        global metadata_filename
        txt = dir(fullfile(mainfolder,subfolder,position,'*.txt'));% lists all files with .txt ending
        metadata_filename = fullfile(txt.folder,txt.name);

        for t = 1:n_images % iterating over each time point in the position folder
            % add time information to the data structure
            [data] = fun_time_V4(data,t); 
            
            % run the function to analyze the droplets 
            [data] = fun_droplets_V4(data,t,level,lim1,lim2,lim3,troubleshoot);        
        end
        
        % add results for all times "t" in position "p" to alldata
        if j == 1
            alldata.(position) = data;
        else
            alldata.(position) = [alldata.(position);data];
        end
    end
end
% Say that the analysis is compete
    disp('Analysis Complete')
    
%% Adding some additional things to the data structure
% 1. elapsed time from first "initial" frame
% 2. Diameter corrected for pixel size
% 3. Volume and Concentration at each frame

for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    
    % 1. adding the elapsed time for each frame, based on the first recieved
    % time in the data structure. In this case, the first received time is the
    % time of the first frame in "initial" folder
    t0 = alldata.(position)(1).r_time;
    t0 = datevec(t0);
    for i = 1:length(alldata.(position))
        t = alldata.(position)(i).r_time;
        t = datevec(t);
        dt = etime(t,t0)/60; % calculate elapsed time from the first image in minutes
        % add to the data structure
        alldata.(position)(i).('dt') = dt;
    end
    
    % 2. Adding the corrected diameter (with correct pixel size)
    for i = 1:length(alldata.(position))
        diameter = alldata.(position)(i).diameter;
        pixel_size = alldata.(position)(i).pixel_size;
        real_diameter = diameter.*pixel_size;
        % add to the data structure
        alldata.(position)(i).('real_diameter') = real_diameter;
    end
    
    % 3. calculating the volume and concentration
    for i = 1:length(alldata.(position))
        D0 = alldata.(position)(i).real_diameter; %dropdiameter in um 
        V0=(1+3/2./h.^2.*((D0-h).*(D0-0.5.*(2-pi).*h)));
        V10 = V0.*(4/3*pi.*(h/2).^3); %drop volume in um^3
        % add to the data structure
        alldata.(position)(i).('volume') = V10;
        % calculating the concentration change
        initial_volume = alldata.(position)(1).('volume');
        alldata.(position)(i).('concentration') = C0*initial_volume./V10;
    end

end

%% finding errors in the data, flagging bad bubbles, etc
    % iterate through all of alldata and if the errorcode entry is not "no
    % error" then that row of data is removed.

%% plotting 
fig = figure;
set(fig, 'WindowStyle', 'Docked');  

hold on

for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array

    % choosing data to plot
    t = [alldata.(position).dt];
    c = [alldata.(position).concentration];
    
    % plotting and marker choice, linestyle, etc
    plot(t,c,'o')
   
end
% add legend
lgd = legend(fieldnames(alldata),'location','best');

% formatting the plot
xlabel('Time [min]')
ylabel('Concentration [M]')


