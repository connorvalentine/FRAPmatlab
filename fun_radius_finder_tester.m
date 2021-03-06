%% Function description: 
% readins in the first image from the frap dataset, finds the bleached
% circle and produces some plots showing the line profile across the found
% cicle

%% inputs and outputs
% inputs: id structure and fits structure. Save_im is 'y' or 'no'

function [struct_out1,struct_out2,fig_out,stats] = fun_radius_finder(id,fits,save_im)
% define the global variables
global boxfolder folder1 folder2 plotfolder npoints

% define other parameters 


for field = fieldnames(id)' % iterate through the position list in id structure
    position = field{1};
    plotflag = 'y';
    
    % First we load the prebleach image
    folder3 = 'prebleach';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    f1 = fullfile(data(1).folder,data(1).name); % loads the pre-bleach image
    im0 = imread(f1);    
    im0sm = imgaussfilt(im0,2); % smoothed prebleach image
    fits.(position).('pbim') = im0; % add info to fits.
    
    % Second we load the first image after bleaching
    folder3 = 'frap';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    f1 = fullfile(data(2).folder,data(2).name); % Load first bleached frame - frame 1 is the image of the laser.
    im1 = imread(f1);
    im1sm = imgaussfilt(im1,2); % now smooth image
    fits.(position).('imbds') = [min(im1(:)),mean2(im1)]; % parameters for plotting images
    
    imTEST = 10000*(1- im1sm./im0sm); % because each image is an int, the imTest is all int, rounded to either 0 or 10,000
    imlbtest = min(imTEST(:));
    imubtest = mean2(imTEST);
    imlb = min(im1(:));
    imub = mean2(im1);
     % First we do some pre-analysis on the images 
    T = graythresh(imTEST); 
    bw1 = im2bw(imTEST,T); % old way to do it seems to work better
    bw2 = bwareaopen(bw1,300); % remove specks and dots smaller than lim2 pixels in area
    bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area
    
    % filling in the boundaries in black and white image to make solid shapes
    
    % filling in the boundaries in black and white image to make solid shapes
    [B,L] = bwboundaries(bw3,'noholes');
    
    % These objects are now analyzed by regionprops(). 
    % The output is a table, called stats, that has the information for each object
    stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','PixelList');
    fits.(position).('objects_in_image') = stats;        % the main object is the most solid one usually
    
    % now we iterate through the ojects found and choose the best one. This
    % could be improved some what I think
    lim1 = 50; % minimum pixel radius to look for on most solid droplet 
    if height(stats) == 0 | stats.EquivDiameter(find(stats.Solidity == max(stats.Solidity)))/2 < lim1  
        fits.(position).('center') = [1024,1024]; % center is the middle of the image
        fits.(position).('radius') = 10;% radius set to 10 pixels (too small to be a real object)
    else
        main_idx = find(stats.Solidity == max(stats.Solidity));
        fits.(position).('radius') = (stats.EquivDiameter(main_idx)/2);
        fits.(position).('center') = stats.Centroid(main_idx,:); 
    end
    fits.(position).('pixel_size') = id.(position).pixel_size_manual;

end

% now second loop to correct for weird centers by using other centers
counter = 0; % position #(1-3) at current sample, is tracked to save images, etc properly.
for field = fieldnames(id)' % iterate through the position list in id structure
    position = field{1};
    stats = fits.(position).objects_in_image;
    
    % now we sort through the objects to make sure the bleached circle is found
    % if no droplets are found (height of table == 0) or the most solid
    % droplet is very small (<lim1) then we advance with artificial values
    % for the center of the droplet.
    lim1 = 50; % minimum pixel radius to look for on most solid droplet 
    if height(stats) == 0 | stats.EquivDiameter(find(stats.Solidity == max(stats.Solidity)))/2 < lim1
        if height(stats) == 0
            disp(['no droplets found', position])
        else
            disp(['no droplets found at ', position]) 
        end 

        temp_pos_list = fieldnames(id)';
        for i = 1:length(fieldnames(id))
            temp_pos = temp_pos_list{i};
            temp_stats = fits.(temp_pos).objects_in_image;
            
            if height(temp_stats) == 0 | temp_stats.EquivDiameter(find(temp_stats.Solidity == max(temp_stats.Solidity)))/2 < lim1
                continue
            else 
                main_idx = find(temp_stats.Solidity == max(temp_stats.Solidity));
                fits.(position).('radius') = (temp_stats.EquivDiameter(main_idx)/2);
                fits.(position).('center') = temp_stats.Centroid(main_idx,:); 
            end
        end
        
        disp(['using center of ',temp_pos, ' instead']);
    else 
    end    
    
    % now we take some image profiles based on the center of the laser hole
    % First we make a vertical line profile across the radius
    center = fits.(position).center;
    radius = fits.(position).radius;
    profile_line_x = [center(1) center(1)];
    profile_line_y = [0 size(im1,1)];
    
    % Now we record the profile on the first bleached image
    prof = improfile(im1sm,profile_line_x,profile_line_y); % line profiles 
    profile1 = prof(2:end);
    
    % Now we record the profile of the same line, but on the pre-bleached
    % image
    prof = improfile(im0sm,profile_line_x,profile_line_y); % line profiles 
    profile0 = prof(2:end);  
    
    % temporary photobleaching normalization of profile 1 by profile 0 for the fit
    norm_ratio = mean(profile0(250:500))/mean(profile1(250:500));
    
    % now we normalize profile1 by profile0 and flip 
    y_norm = 1- (profile1./profile0)*norm_ratio;
    x_norm = linspace(1,2048,2048)'; 
    
    f = fit(x_norm,y_norm,'gauss1'); %number after gauss changes number of peaks in guass
    % a = height of peak
    % b = position of the center of the peak
    % c = sigma*sqrt(2)
    FWHM = 2*sqrt(log(2))*f.c1;
    peak_height = f.a1;
    peak_center = f.b1;
    
    ci = confint(f,0.95);
    FWHM_err = 2*sqrt(log(2))*abs(ci(:,3));
    
    % make a circle mask defined by the FWHM of the bleached area
    circleCenterX = center(1);
    circleCenterY = center(2); % might have these flipped
    circle_mask = false(2048,2048);
    [x_im, y_im] = meshgrid(1:2048,1:2048);
    circle_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (FWHM/2).^2) = true; 
    
    % now make a concentric circle mask to define the reference region
    r_out = 4;
    r_in = 2.5;
    reference_mask = false(2048,2048);
    [x_im, y_im] = meshgrid(1:2048,1:2048);
    reference_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (r_out*FWHM/2).^2 ...
        & (x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 >= (r_in*FWHM/2).^2) = true; 
    
    % now we mask the original image to only see pixels we want in the
    % prebleach image to find the intensity in the circle and the reference
    % region
    masked_image0 = im0; % Initialize with the entire image.
	masked_image0(~circle_mask) = 0; % Zero image outside the circle mask.
    I_t0 = mean(masked_image0(masked_image0 > 0));    
    
    % mask the prebleach image once more to find the reference region
    masked_reference_image0 = im0;
    masked_reference_image0(~ reference_mask) = 0;
    ref_0 = mean(masked_reference_image0(masked_reference_image0 > 0)); % find mean of all pixels that are not == 0
    
    % masking image 1 just for the plot
    masked_image1 = im1; % Initialize with the entire image.
	masked_image1(~circle_mask) = 0; % Zero image outside the circle mask.
    masked_reference_image1 = im1;
    masked_reference_image1(~ reference_mask) = 0;    
    
    % find the mean value of the pixels in this circle
    % add peak info to the fits 
    fits.(position).('FWHM') = FWHM;
    fits.(position).('dFWHM') = FWHM_err;
    fits.(position).('peak_height') = peak_height;
    fits.(position).('peak_center') = peak_center;
    fits.(position).('profile0') = profile0;
    fits.(position).('profile1') = profile1;
    fits.(position).('profile_line_x') = profile_line_x;
    fits.(position).('profile_line_y') = profile_line_y;
    fits.(position).('circle_mask') = circle_mask;
    fits.(position).('reference_mask') = reference_mask;
    fits.(position).('I_t0') = I_t0;
    fits.(position).('ref_0') = ref_0;

    if counter == npoints
        counter = 0;
    else
        counter = counter +1;
    end
    
    
    fig = figure('name',position,'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
%         imlb = fits.(position).imbds(1);
%         imub = fits.(position).imbds(2);

        subplot(2,2,1);
        %show the first image after photobleaching
            hold on 
            imshow(imTEST,[imlbtest,imubtest],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            plot(profile_line_x,profile_line_y,'m')
            title('image1 profile line')
            hold off 
        if plotflag == 'y'
        subplot(2,2,2)
        %show the laser region k circles
            hold on
            % first plot the black and white image
%             imshow(masked_image1,[imlb,imub],'Border','tight','InitialMagnification', 'fit')
            imshow(masked_reference_image1+masked_image1,[imlb,imub],'Border','tight','InitialMagnification', 'fit')
            h= viscircles(center,radius,'linewidth',0.2,'color','g');
            plot(profile_line_x,profile_line_y,'m')
            title('Gr imfind V. regions')
            hold off 
        subplot(2,2,3)
            hold on 
            plot(profile1,'b')
            plot(profile0,'g')
            title('prebleach v. im1')
            axis([0 2024 min(profile1)*0.8 max(profile0)*1.2])
        subplot(2,2,4)
            hold on 
            plot(f,x_norm,y_norm)
            errorbar(peak_center,peak_height/2,FWHM/2,'horizontal')
            title('norm im1 & Guass')
            axis([0 2024 0  1])
        else
        end
    % save images if selected above
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 6];             % define location to save the images 
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_frames_',num2str(counter)];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png  
    else    
    end
end  
    
%designate the outputs properly
struct_out1 = id;
struct_out2 = fits;
fig_out = fig;
end