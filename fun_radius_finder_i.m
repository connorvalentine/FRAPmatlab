%% inputs and outputs
% inputs: id structure and fits structure. Save_im is 'y' or 'no'
% need to import position as well
function [center,ref_i,I_ti,ref_0_i,drift_flag] = fun_radius_finder_i(position,imi,fits,previous_center_coords,previous_norm_ratio,t)
    drift_flag = 0;
    im0 = fits.(position).im0;
    im0sm = fits.(position).im0sm2;
    radius = fits.(position).radius; % using same radius at all time points, just changing the center
    %Next we load the image we want to use after bleaching
    imism = imgaussfilt(imi,69).*previous_norm_ratio; % now smooth image, sorta intensity adjusted for photobleaching
    
    imTEST = 10000*abs(1- round(double(imism)./double(im0sm),1));% because each image is an int, the imTest is all int, rounded to either 0 or 10,000
    % First we do some pre-analysis on the images 
    T = graythresh(imTEST);
%     T2 = graythresh(imTEST2);
    bw1 = imbinarize(imTEST,'adaptive','ForegroundPolarity','dark','Sensitivity',T);
%     bw2 = bwareaopen(bw2,10000); 
%     bw1 = im2bw(imTEST,T); % old way 
    bw3 = bwareaopen(bw1,10000); % (may need to be reverseD?)remove specks and dots smaller than 10000 pixels in area
%     [B,L] = bwboundaries(bw3,'holes'); % filling in the boundaries in black and white image to make solid shapes
    
    % These objects are now analyzed by regionprops(). 
    % The output is a table, called stats, that has the information for each object
    stats = regionprops('table', bw3, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','PixelList');

    % finding the distance of each object in the image from the previous
    % image
    dist = ones(height(stats),1);
    for k = 1:height(stats)
        cent = stats.Centroid(k,:);     % compute the distance from center to center of each shape
        dist(k) = sqrt((previous_center_coords(1) - cent(1)).^2 + (previous_center_coords(2) - cent(2)).^2);
    end
    stats.('distance') = dist;
    
    stats2 = stats(find(stats.Area < 500000),:);
    
    main_idx = find(stats2.distance == min(stats2.distance));

    if height(stats2) == 0 | stats2(main_idx,:).distance > 200
        center = previous_center_coords; % center is the middle of the image
        disp('center moved too far from previous frame, skipping to next')
        disp(t)
        disp(stats2(main_idx,:))
    else
        center = stats2.Centroid(main_idx,:);
    end
    
    % make a circle mask defined by the FWHM of the bleached area
    circleCenterX = center(1);
    circleCenterY = center(2); % might have these flipped
    circle_mask = false(2048,2048);
    [x_im, y_im] = meshgrid(1:2048,1:2048);
    circle_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (radius).^2) = true; 
        
    % now make a concentric circle mask to define the reference region
    r_out = 4;
    r_in = 2.5;
    reference_mask = false(2048,2048);
    [x_im, y_im] = meshgrid(1:2048,1:2048);
    reference_mask((x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 <= (r_out*radius).^2 ...
        & (x_im - circleCenterX).^2 + (y_im - circleCenterY).^2 >= (r_in*radius).^2) = true; 
    
    % check to make sure it is not approaching the edge of the circle
    if abs(1 - circleCenterX) < (r_out*radius) || abs(2048 - circleCenterX) < (r_out*radius) 
        disp('center too close to edge of frame for analysis')
        drift_flag = 1;
    elseif abs(1 - circleCenterY) < (r_out*radius) || abs(2048 - circleCenterY) < (r_out*radius) 
        disp('center too close to edge of frame for analysis')
        drift_flag = 1;
    else
    end 
    % masking image 1 just for the plot
    masked_image_i = imi; % Initialize with the entire image.
	masked_image_i(~circle_mask) = 0; % Zero image outside the circle mask.
    I_ti = mean(masked_image_i(masked_image_i > 0));
    
    masked_reference_image1 = imi;
    masked_reference_image1(~ reference_mask) = 0;    
    ref_i = mean(masked_reference_image1(masked_reference_image1 > 0));

    masked_reference_image0 = im0;
    masked_reference_image0(~ reference_mask) = 0;
    ref_0_i = mean(masked_reference_image0(masked_reference_image0 > 0)); % find mean of all pixels that are not == 0
     
    if t == 2 || mod(t,2) == 0
         
        subplot(1,2,1) %show the first image after photobleaching 
            imshow(imi,[min(imi(:)),max(imi(:))],'Border','tight','InitialMagnification', 'fit');
            hold on
            viscircles(center,radius,'linewidth',0.2,'color','g');
            title(["image at time " + string(t)])
            hold off 
        subplot(1,2,2) % plot the black and white image
            imshow(bw3,[0,1],'Border','tight','InitialMagnification', 'fit');
            hold on
            viscircles(center,radius,'linewidth',0.2,'color','g','LineStyle',':');
            for k = 1:height(stats2)
                viscircles(stats2.Centroid(k,:),stats2.EquivDiameter(k,1)/2,'linewidth',0.2,'color','b');
            end
            hold off
            title(["analyzed at time " + string(t)])

        drawnow % necessary or it wont plot
     else
     end
            
end