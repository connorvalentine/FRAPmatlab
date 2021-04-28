function [alldata] = fun_frap_frames(id,fits,alldata)
global boxfolder folder1 folder2 makemovies plotfolder npoints
folder3 = 'frap';
counter = 0;

for field = fieldnames(id)'
    position = field{1};
    
    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
  
    % pull out the prebleach image info from fits structure
    ref_0 = fits.(position).ref_0;
    I_t0 = fits.(position).I_t0;
    circle_mask = fits.(position).circle_mask;
    reference_mask = fits.(position).reference_mask;
    
    if counter == npoints
        counter = 0;
    else
    end
    counter = counter +1;
    
    if makemovies == 'y'
        disp(['Making video for ', position])
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        
        movie_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_movie_',num2str(counter)];
        movie_path = fullfile(plotfolder,movie_name); 
        figure(69);
        set(gcf, 'Color','white')
        set(gca, 'nextplot','replacechildren', 'Visible','off');

        %# create AVI object
        vidObj = VideoWriter([movie_path '.avi']);
        vidObj.Quality = 100;
        vidObj.FrameRate = 5;
        open(vidObj);
        %# create movie
        video_mask = circle_mask + reference_mask;
    else 
    end

    % iterate through the images in data structure
    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time(data,t); 
       % add elapsed time from received time of first frap frame to data struct.
        time1 = data(1).r_time;
        time1 = datevec(time1);
        time_i = data(t).r_time;
        time_i  = datevec(time_i);
        dt = etime(time_i,time1); % calculate elapsed time from the first image in seconds
        data(t).('dt') = dt; %save to data structure
        
        % calculate pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        imi = double(imread(f1));
        if makemovies == 'y'
            if mod(t,5) == 0 || t == 1 
                imbds = fits.(position).imbds;
                vid_im = imi;
                lb = min(vid_im(:));
                ub = mean2(vid_im);
                vid_im(~video_mask) = vid_im(~video_mask)*0.95; % darken image outside the circle and reference mask.
                imshow(vid_im,[lb,ub],'Border','tight','InitialMagnification', 'fit');
                hold on
                time = round(data(t).dt/60,1);
                timetext = [num2str(time),' min'];
                text(1000,200,timetext,'color','m');
                writeVideo(vidObj, getframe(gca));
                hold off
            else 
            end
        else
        end
        
        % find mean intensity inside of the circle
        masked_image_i = imi; % Initialize with the entire image.
        masked_image_i(~circle_mask) = 0; % Zero image outside the circle mask.
        I_ti = mean(masked_image_i(masked_image_i > 0));

        % calculate mean intensity inside the reference region
        masked_reference_imagei = imi;
        masked_reference_imagei(~ reference_mask) = 0;
        ref_i = mean(masked_reference_imagei(masked_reference_imagei > 0));

        norm_ratio = ref_0/ref_i;
        
        if ref_i ==0
            ref_i = 0.001;
            disp(position)
            disp('reference region mistake at frame #' +string(t))
        else
        end
        % normalized by prebleach Intensity(I_t0)
        % double normalized by the intensity ratio in the reference region
        % (refI_t0) divdided by the ref region ROI at time ti refI_ti
        normalized_i = (I_ti./I_t0).*(norm_ratio); 

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('ref_i') = ref_i;
        data(t).('ref_ratio') = norm_ratio;
        data(t).('IN_ti') = normalized_i;
    end
    % add elapsed time from received time of first frap frame to data struct.
%     time1 = data(1).r_time;
%     time1 = datevec(time1);
%     for i = 1:length(data)
%         t = data(i).r_time;
%         t = datevec(t);
%         dt = etime(t,time1); % calculate elapsed time from the first image in seconds
%         data(i).('dt') = dt; %save to data structure
%     end
        
    if makemovies == 'y'
        close(gcf)
        %# save as AVI file, and open it using system video player
        close(vidObj);
    else 
    end
    
    
    % add the completed data structure to alldata
    alldata.(position) = data;
    disp([position,' loaded'])
end
end
