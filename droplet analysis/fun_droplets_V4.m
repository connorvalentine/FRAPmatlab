%% function description
% this function takes in an image of a droplet in a microfluidic device,
% and outputs information about the droplet. Lim1,Lim2,Lim3, and level can
% all be adjusted to tailor how the image is analyzed, and how objects in
% the image are classified as circles or not.

%% inputs
% data - the "data" structure for the position you want to analyze
% t - the frame number of the image within the position folder
% lim1 - the circularity cut-off of the dropplets
% lim2 - the solidity cut-off for the droplets
% lim3 - the minimum area cut-off to ignore small objects
% level - the black and white cutoff when binarizing the image
% troubleshoot - to make the plots with bubble images or no (1 == yes, 0 == no)

%% outputs
% data - we add a lot of information about the droplets, and every object
%        found. also add these to the data structure:
%      - Two tables, drop_info and stats, are added to index, t, in data 
%      - figure_out - a plot with the original droplet images along with overlaid
%           circles. This will only be output if troubleshoot ==1.
%% main part of the function
function [data] = fun_droplets_V4(data,t,level,lim1,lim2,lim3,troubleshoot)
% loading in the image 
frame_filename = fullfile(data(t).folder,data(t).name);
image = imread(frame_filename);

% First we do some pre-analysis on the images
bw1 = im2bw(image,level); % make image black and white
bw1 = imcomplement(bw1);% flip black and white    
bw2 = bwareaopen(bw1,200); % remove specks and dots smaller than lim2 pixels in area
bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area

% filling in the boundaries in black and white image to make solid shapes
[B,L] = bwboundaries(bw3,'noholes');

% These objects are now analyzed by regionprops(). 
% The output is a table, called stats, that has the information for each object
stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area');
    
%% Next we sort the objects into things we think are droplets, and those that are not
drop_info = table(); % initialize the table to contain information for objects we think are droplets
count = 0;
    
% Objects in stats table are classified as a droplet if they meet all of these
% critera:
    % their "circularity" is > lim 1 
    % and their solidity is > lim 2
    % and their diameter is larger than a minimum value (lim3)
for d = 1:height(stats) %height(stats) is # of objects found by regionprops
    if stats.Circularity(d) > lim1 && stats.Solidity(d) > lim2 && stats.EquivDiameter(d) > lim3
        count = count+1; % count is # of objects that meet all 3 criteria
        drop_info(count,:) = stats(d,:); % take the row of the stats table and add it to drop_info table
    else
    end
end

% for some images we will not be able to find any objects that are droplets  
if height(drop_info) == 0 % if no objects meet our 3 criteria, the drop_info table will be empty
    disp('no good droplet found') %output to command window
    data(t).('diameter')= [];
    data(t).('circularity')= [];
    data(t).('solidity')= [];
    data(t).('drop_center')= [];
    center = size(image)/2; % center is the middle of the image
    radius = 10; % radius set to 10 pixels (too small to be a real droplet)
else % if the drop_info table is not empty
    
% This is the center and radii of every "good" object found. The
% code works even if multiple "good" drops are found. This way
% during troubleshooting we can tune the parameters until "good
% drops" are actually good.
% adding what we found to the data structure 
data(t).('diameter')= drop_info.EquivDiameter(1);
data(t).('circularity')= drop_info.Circularity(1);
data(t).('solidity')= drop_info.Solidity(1);
data(t).('drop_center')= drop_info.Centroid(1,:);

center = drop_info.Centroid(:,:); 
radius = drop_info.EquivDiameter/2;
end

% including the drop_info table and stats table in case we need it later
data(t).('drop_info') = drop_info;
data(t).('object_stats') = stats;

%% output plots for troubleshooting mode
% labels: a vector that is 1 --> length of drop info.
% Used to number each circle that is plotted.
labels = string(linspace(1,height(drop_info),height(drop_info)));
        
% if you want to output the plots, troubleshoot == 1
if troubleshoot == 1                
    fig = figure('visible','off');
    set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
%     set(fig, 'Number',1);
    subplot(1,2,1);
        % this subplot plot is just the un-anlayzed image from
        % micromanager with the final droplet circles overlaid 
        hold on % allows multiple plots on same figure

        % first plot the image 
        I=double(image);
        i_min = min(I(:)); % rescale brightness by min 
        i_max = max(I(:)); % rescale brightness by max
        imshow(I,[i_min,i_max*0.3],'Border','tight','InitialMagnification', 'fit');

        % plotting the results to confirm. 
        % Viscircles() plots circles on top of the image for each center and radii
        viscircles(center,radius,'linewidth',0.2,'color','g');

        % lable the circles by their row number in drop_info
        text(center(:,1),center(:,2),labels,'Color','m','FontSize',18,...
            'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'middle');

        % Increase the size and position of the image 
        pos1 = get(gca, 'Position');
        pos1(1) = 0; % set x position of the image to left side
        pos1(2) = 0.2; % set y to a little bit higher
        pos1(3) = 0.5; % set width to 1/2 of available space
        set(gca,'Position',pos1);

        % reset hold 
        hold off 

    subplot(1,2,2)
        hold on
        % this plot is the image and binarizaiton, and iamge processing

        % first plot the black and white image
        imshow(bw3,'Border','tight','InitialMagnification', 'fit')

        % add the boundaries from regionprops, with colors
        imshow(label2rgb(L,@jet,[.5 .5 .5]),'Border','tight','InitialMagnification', 'fit');

        % plot the drop_info results to show where we think the
        % "good drops" are located
        viscircles(center,radius,'linewidth',0.2,'color','g');
        text(center(:,1),center(:,2),labels,'Color','m','FontSize',18,...
            'FontWeight','bold','FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'middle'); 
        hold off   

        % Increase the size and position of the image 
        pos2 = get(gca, 'Position');
        pos2(1) = 0.5;
        pos2(2) = 0.2; % set y to a little bit higher
        pos2(3) = 0.5;
        set(gca,'Position',pos2);
else % if troubleshoot is set to 0, no plot is made
    fig = 'no plot';
end

% add the figure to the data structure
data(t).('figure') = fig;
end