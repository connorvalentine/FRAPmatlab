function [alldata] = fun_frap_frames(id,fits,alldata)
global boxfolder folder1 folder2
folder3 = 'frap';

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

    % iterate through the images in data structure
    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time(data,t); 
   
        % calculate pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        imi = double(imread(f1));
        
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
    time1 = data(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(data)
        t = data(i).r_time;
        t = datevec(t);
        dt = etime(t,time1); % calculate elapsed time from the first image in seconds
        data(i).('dt') = dt; %save to data structure
    end
    
    % add the completed data structure to alldata
    alldata.(position) = data;
    disp([position,' loaded'])
end
end
