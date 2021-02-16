%% function description
% this function takes in an image of a droplet in a microfluidic device,
% and outputs information about the droplet. Lim1,Lim2,Lim3, and level can
% all be adjusted to tailor how the image is analyzed, and how objects in
% the image are classified as circles or not.

%% inputs
% lim1 is the circularity cut-off of the dropplets
% lim2 is the solidity cut-off for the droplets
% lim3 is the minimum area cut-off to ignore small objects
% level is the black and white cutoff when binarizing the image
% troubleshoot is to make the plots with bubble images or no (1 == yes, 0 == no)

%% outputs
% stats is a table with information about every object in the image found by
%           regionprops().
% drop_info is a table with the rows from stats that meet the criteria in
%           lim1,lim2 and lim3.
% center is a list of centers for objects that meet the criteria (also included in drop_info)
% radius is a list of redii for objects that meet criteria lim1-lim3 (also
%           included in drop_info
% figure_out is a plot with the original droplet images along with overlaid
%           circles. This will only be output if troubleshoot ==1.
%% main part of the function
function [stats,drop_info,center,radius,error,figure_out] = fun_dropletsV2(image,level,lim1,lim2,lim3,lim4,troubleshoot)
%% First we do some pre-analysis on the images
% first we Binarize the image into black and white 
    % image pre-processing steps
    %     bw1 = imbinarize(image,level);
    bw1 = im2bw(image,level); % make image black and white
    bw1 = imcomplement(bw1);% flip black and white    
    bw2 = bwareaopen(bw1,200); % remove specks and dots smaller than lim2 pixels in area
    bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area

% filling in the boundaries in the black and white image to make solid
% shapes
    [B,L] = bwboundaries(bw3,'noholes');

% These objects are now analyzed by regionprops(). The output is a table,
% called stats, that has the information for each object
    stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','Eccentricity');
    
%% Next we sort the objects into things we think are droplets, and those that are not
    drop_info = table(); % initialize the table to contain information for objects we think are droplets
    
    % Objects in stats table are classified as a droplet if they meet all of these
    % critera:
    % their "circularity" is > lim 1 
    % and their solidity is > lim 2
    % and their diameter is larger than a minimum value (lim3)
    

        A = find(stats.EquivDiameter > lim3 );

        if length(A)== 0
           Error = 1;  % 'no drop is in trap'

        else C= find (stats.Circularity(A) > lim1 & stats.Eccentricity(A) < lim4 & stats.Solidity(A) > lim2);
        end
        if length (C) == 0 
            Error =2; %'weird shape or wetting to the channel';
        else drop_info (1,:)= stats(A(C),:); %'good' object found and store its info into drop_info
             Error = 0; %'good' object found
        end
    
    % for some images we will not be able to find any objects that are
    % droplets 
    
    if height(drop_info) == 0 % if no objects meet our 3 criteria, the drop_info table will be empty
        disp('no good droplet found') %output to command window
            center = size(image)/2; % center is the middle of the image
            radius = 10; % radius set to 10 pixels (too small to be a real droplet)
    else % if the drop_info table is not empty
        % This is the center and radii of every "good" object found. The
        % code works even if multiple "good" drops are found. This way
        % during troubleshooting we can tune the parameters until "good
        % drops" are actually good.
        center = drop_info.Centroid(:,:); 
        radius = drop_info.EquivDiameter/2;
    end
        
%% output plots for troubleshooting mode
% labels is a vector that is 1 --> length of drop info. This will be used
% to number each circle that is plotted.

        labels = string(linspace(1,height(drop_info),height(drop_info)));
        
% if you want to output the plots, troubleshoot == 1
    if troubleshoot == 1                
        
            F = figure(1); % start the figure where we will plot the images
                set(F, 'WindowStyle', 'Docked');  %figure will dock instead of free float

            subplot(1,2,1);
                % this plot is just the un-anlayzed image, straight from
                % micromanager with the final droplet circles overlaid to
                % make sure that the algorithm is working 
               
                % first plot the image 
                I=double(image);
                i_min = min(I(:)); % rescale brightness by min 
                i_max = max(I(:)); % rescale brightness by max
                imshow(I,[i_min,i_max*0.3],'Border','tight','InitialMagnification', 'fit');
    
                % plotting the results to confirm. Viscircles() plots
                % circles on top of the image for each center and radii
                hold on % lets us keep plotting on same figure
                viscircles(center,radius,'linewidth',0.2,'color','g');
                
                % lable the circles by their row number in drop_info
                text(center(:,1),center(:,2),labels,'Color','m','FontSize',18,...
                    'FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'middle');
                
                %lable error type 
                    txt = ['Error: ' num2str(Error)];
                    text (100,100,txt,'Color','red','FontSize',18);
                
                % Increase the size and position of the image so that is is
                % larger
                pos1 = get(gca, 'Position');
                pos1(1) = 0;
                pos1(3) = 0.5;
                set(gca,'Position',pos1);
                
                % reset hold 
                hold off 

            subplot(1,2,2)
                % this plot is the image and binarizaiton, and iamge processing
               
                % first plot the black and white image
                imshow(bw3,'Border','tight','InitialMagnification', 'fit')
                
                % add the boundaries from regionprops, with colors
                imshow(label2rgb(L,@jet,[.5 .5 .5]),'Border','tight','InitialMagnification', 'fit');
                hold on 
                
                % plot the drop_info results to show where we think the
                % "good drops" are located
                
                viscircles(center,radius,'linewidth',0.2,'color','g');
                text(center(:,1),center(:,2),labels,'Color','m','FontSize',18,...
                    'FontWeight','bold','FontWeight','bold','HorizontalAlignment', 'center','VerticalAlignment', 'middle'); 
                hold off   
                
                % Increase the size and position of the image so that is is
                % larger
                pos2 = get(gca, 'Position');
                pos2(1) = 0.5;
                pos2(3) = 0.5;
                set(gca,'Position',pos2);
    else % if troubleshoot is set to 0, no plot is made
        F = 'no plot';
    end
    
%% list our function outputs 
drop_info = drop_info;
stats = stats;
center = center;
radius = radius;
error = Error;
bw_image = bw3;
figure_out = F;
end