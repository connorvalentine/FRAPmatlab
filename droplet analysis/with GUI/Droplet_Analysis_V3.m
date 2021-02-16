%% micromanager tif folder analysis of the metadata
% Connor Valentine and Deyu Yang
% 2020 Walker Group

%% notes for things to add
% 1. button to abort the script during troubleshootin
% 2. explain circularity and solidity
% 3. button to flad a droplet result as bad so that it is not included in
% the analysis
% show image of intersted drop



%% First section to reset the script, variables, and plot defaults
% Reset the environment
    clear all;
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
    n = 11;% number of positions 
    level = 0.2; % black and white cutoff for the bright field images
    level1 = 0.02; % black and white cutoff for the polarized images
    lim1 = 0.1; % cutoff for "circularity" parameter. Will ignore drops with circularity less than lim1
    lim2 = 0.5; % cutoff for "solidity" parameter
    lim3 = 100; % ignore "droplets" smaller than 100 pixels in diameter
    lim4 = 0.3; % cutoff for "Eccentricity" parameter
    h = 93;% device height in um
    C0 = 0.464; % initial solute concentration (in M)
    troubleshoot = 0; % make this = 1 if you want to output plots to click through to confirm circle finding

%% Before we do any analysis, specify what folders we are using
    % this line of code lists the filepath for the folder that this matlab
    % script is currently in. The user should not have to input the
    % filepath by hand
    mainfolder = cd; 
    subfolder1 = 'initial';
    subfolder2 = 'dehydration';

%% Analysis - iterating over the 'initial' and 'dehydration folder

for j=1:2 
    if j == 1
        subfolder = subfolder1;
        level=level;
        n_images = 1;
    else
        subfolder = subfolder2;
        level = level1;  
        n_images = 10;
        %n_images = length(data); %# of images
    end

    %iterating over each position and time 

    for p = 1:n
        % First we make the folder name for the position we want
        prefix = 'Pos';
        position = [prefix,num2str(p-1)];

        % Make a list of each image name in our position folder
        list = dir(fullfile(mainfolder,subfolder,position,'*.tif')); % lists all files with .tif ending in the position folder
        data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
        %n_images = length(data); %# of images

        % Information about the metadata file
        txt = dir(fullfile(mainfolder,subfolder,position,'*.txt'));% lists all files with .txt ending
        metadata_filename = fullfile(txt.folder,txt.name);

        for t = 1:n_images 
            % Specify file name for the image we want to analyze
            frame_filename = fullfile(data(t).folder,data(t).name);
            image = imread(frame_filename); % load image

            % run the function to analyze the droplets 
            [stats,drop_info,center,radius,error,figure_out] = fun_dropletsV3(image,level,lim1,lim2,lim3,lim4,troubleshoot);

            % this block of code will open up each figure when troubleshooting.
            % You will have to click the button to confirm each frame is good
            if troubleshoot == 1
                disp(position + ", " + "time point = " + string(t))
                disp(drop_info)
                Button = uicontrol('Style','PushButton','String','Next Image',...
                                         'Position',[20 20 200 60],'fontsize',20,...
                                         'Callback', 'uiresume(figure(1))');
                uiwait(figure(1));
            else
            end

            % find the time for each frame 
            [e_time,r_time,pixel_size] = fun_time(metadata_filename,data(t).name);
                % e_time is elapsed time from start of experiment
                % r_time is the time the picture was taken

            % add the results to our data structure
            data(t).('dt') = e_time/1000; % elapsted time from start of experiment in seconds
            data(t).('r_time') = r_time; % recieved time according to the computer clock
            data(t).('pixel_size') = pixel_size; % microns and assumes pixel-size calibration was run
            data(t).('diameter')= radius*2;   
            data(t).('error') = error;
            
            % get time between current frame to the initial frame
            if j == 1
                data(t).('time')=datevec(r_time);
            else
                data(t).('time')=(etime(datevec(r_time),maindata.initial.(position)(1).time))/60;
            end
            
            % Calculate droplet volume
                B0 = radius*2*pixel_size;  %dropdiameter in um 
                V0=(1+3/2./h.^2.*((B0-h).*(B0-0.5.*(2-pi).*h)));
                V10 = V0.*(4/3*pi.*(h/2).^3); %drop volume in um^3
                data(t).('volume')=V10;
             % Calculate Concentration
             if j==1
                 data(t).('concentration')= C0;
             else
                 data(t).('concentration')= C0*maindata.initial.(position)(1).volume/V10;
             end
                
                
        end
        % add the results for all times "t" in position "p" to alldata
             alldata.(position) = data;
    end
    % add "initial'' and ''dehyration'' results to maindata
    maindata.(subfolder)=alldata;
end

% Say that the analysis is compete
    disp('Analysis Complete')
    

%% Example for calculating elapsed time using the recieved time
% look up datevec() function in matlab help. Turns a data+time into a 6
% column vector. These date vectors can be subtracted using etime() to
% provide accurate time-spans

t0 = datevec(alldata.Pos0(1).r_time);
t2 = datevec(alldata.Pos0(5).r_time);
dt = alldata.Pos0(5).dt;

% etime(t2,t0) %elapsed time in seconds

%% Example of plotting 
%% Plot d vs. time 
figure(2)
set(figure(2), 'WindowStyle', 'Docked');  
hold on

% this gets the names of each position in the data structure, and convert to a string
Positions = string(fieldnames(alldata)); 

% plotting by iteration through each position
for p = 1:n
    % First specifiy what position we are plotting
    Position = Positions(p);
    
    % choosing data to plot
    x = [0 alldata.(Position).time];
    y = [[maindata.initial.(Position)(1).diameter.*maindata.initial.(Position)(1).pixel_size] [alldata.(Position).diameter].*[alldata.(Position).pixel_size]];
    
    % plotting and marker choice, linestyle, etc
    plot(x,y,'o')
end

% formatting the plot
    xlabel('time [min]')
    ylabel('Drop Diameter [\mum]')
    % automatically label the information according to position name 
    % *** only works because Positions vector is in the same order that we plot the
    % data
    legend(Positions,'location','best')

 %% Plot concentration vs.time
figure(3)
set(figure(3), 'WindowStyle', 'Docked');  
hold on
% this gets the names of each position in the data structure, and convert to a string
Positions = string(fieldnames(alldata)); 

% plotting by iteration through each position
for p = 1:n
    % First specifiy what position we are plotting
    Position = Positions(p);
    
    % choosing data to plot
    x1 = [0 alldata.(Position).time];
    y1= [C0 [alldata.(Position).concentration]];
    
    % plotting and marker choice, linestyle, etc
    plot(x1,y1,'o')
end

% formatting the plot
    xlabel('time [min]')
    ylabel('Concentration [M]')
    % automatically label the information according to position name 
    % *** only works because Positions vector is in the same order that we plot the
    % data
    legend(Positions,'location','best')  
    
  
%% Cleaning the data structure
% going to clean bad drops from alldata based on the value of error, remove drops with error not equals to 0  
Dropidx=zeros(1,n);
for p = 1:n
    Position = Positions(p);
    
    % choosing good drop based on the value of error for both
    y2 = [maindata.initial.(Position)(1).error alldata.(Position)(1:10).error]; 
    %if the first 1/3 of the image frames in a position have error 
    % not equal to zero, then keep that position ; This can be changed based
    % on the data set
    err= find(y2~=0);
    if length(err)==0
        Dropidx(1,p)=p;
    end
end

% Plot good drops
O = find(Dropidx~=0); % Positions where drop is good 
figure(4)
set(figure(4), 'WindowStyle', 'Docked');  
subplot (1,2,1)
hold on
subplot (1,2,2)
hold on
% plotting by iteration through each position
for p = 1:length(O)
    Position = Positions(O(p));
    
    % choosing data to plot
    x3 = [0 alldata.(Position).time];
    y3 = [[maindata.initial.(Position)(1).diameter.*maindata.initial.(Position)(1).pixel_size] [alldata.(Position).diameter].*[alldata.(Position).pixel_size]];
    y4 = [C0 [alldata.(Position).concentration]];
    % plotting and marker choice, linestyle, etc
    subplot (1,2,1)
    plot(x3,y3,'o')
    ylim([200 600])
    xlabel('time [min]')
    ylabel('Drop Diameter [\mum]')
    legend(Positions(O),'location','best','Fontsize',10,'NumColumns',2)
    
    subplot (1,2,2)
    plot(x3,y4,'o')
    xlabel('time [min]')
    ylabel('Concentration [M]')
    
end


