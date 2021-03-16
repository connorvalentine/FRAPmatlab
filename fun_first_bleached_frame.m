%% Function description: 
% readins in the first image from the frap dataset, finds the bleached
% circle and produces some plots showing the line profile across the found
% cicle

%% inputs and outputs
% inputs: id structure and fits structure. Save_im is 'y' or 'no'

function [struct_out1,struct_out2,fig_out] = fun_first_bleached_frame(id,fits,save_im)
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

    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure

    % Load first bleached frame
    f1 = fullfile(data(t1).folder,data(t1).name);
    im = imread(f1);
    imsmooth = imgaussfilt(im,8); % now smooth image
    
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
        disp(['no good droplet found at ', position]) 
        % remove field name from id list
        id = rmfield(id,position);
        center = size(image)/2; % center is the middle of the image
        radius = 10; % radius set to 10 pixels (too small to be a real object)
    else 
        % the main object is the most solid one usually
        main_idx = find(stats.Solidity == max(stats.Solidity));
    end
    
    % Add data to fits structure if it doesn't meet the other errors
    if stats.EquivDiameter(main_idx)/2 < lim1
        % object is too small to be real 
        disp(['no good droplet found at ', position]) 
        id = rmfield(id,position);
        center = size(image)/2; % center is the middle of the image
        radius = 10; % radius set to 10 pixels (too small to be a real object)
    else % the main object chosen is fine. save the circle and reference regions into fits 
        fits.(position).('imbds') = [min(im(:)),mean2(im)]; % parameters for plotting images
        fits.(position).('radius') = (stats.EquivDiameter(main_idx)/2);
        fits.(position).('center') = stats.Centroid(main_idx,:); 
        % reference region dx and dy defined above, convert to pixel # list instead of coords
        idx = stats.PixelList(main_idx);
        idx = cell2mat(idx);
        fits.(position).('idx_bleach') = sub2ind(size(im),idx(:,2),idx(:,1));
        fits.(position).('idx_ref') = sub2ind(size(im),idx(:,2)+dy,idx(:,1)+dx); 
        fits.(position).('pixel_size') = id.(position).pixel_size_manual;
        
        center = fits.(position).center;
        x = [center(1) center(1)];
        y = [0 size(im,1)];
        prof = improfile(imsmooth,x,y); % line profiles 
        prof = prof(2:end);
        [pk,lk,w,p] = findpeaks(prof,x,'MinPeakProminence',0.3,'WidthReference','halfprom','MinPeakDistance',100)
    end


end  

% plotting loop
for field = fieldnames(id)' % iterate through the position list in id structure
    position = field{1};
    %plot figure
    fig = figure('name',position,'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
        imlb = fits.(position).imbds(1);
        imub = fits.(position).imbds(2);
        center = fits.(position).center;
        radius = fits.(position).radius;
        idx_ref = fits.(position).idx_ref;
        % vertical line profile across the radius
        x = [center(1) center(1)];
        y = [0 size(im,1)];

        subplot(2,2,1);
        %show the first image after photobleaching
            hold on 
            imshow(im,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            plot(x,y,'m')
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
        c = improfile(im,x,y); % line profiles 
        c1 = improfile(imsmooth,x,y); % line profiles 
        
            hold on 
            plot(c1(:,1,1),'b')
            % add radius and center
            plot([center(2)-radius,center(2)+radius],[1.2*mean([imub,imlb]),1.2*mean([imub,imlb])])
            axis([0 2024 imlb 1.2*imub])

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
struct_out1 = id;
struct_out2 = fits;
fig_out = fig;
end