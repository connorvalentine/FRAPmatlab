%% Function description: 
% readins in the prebleach images and adds information to data structures
%% inputs and outputs
% inputs: id structure and fits structure. 

function [struct_out] = fun_frap_frames(id,fits,alldata)
% define the global variables
global boxfolder folder1 folder2
folder3 = 'frap';

for field = fieldnames(id)'
    position = field{1};
    
    % pull out the prebleach image info from fits structure
    
    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    
    % iterate through the images in data structure
    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time(data,t); 
   
        % calculate pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        im = double(imread(f1));
        I_ti = mean(im(idx_bleach)); 
        refI_ti = mean(im(idx_ref));
        if refI_ti ==0
            refI_ti = 0.001;
            disp(position)
            disp('reference region mistake at frame #' +string(t))
        else
        end
        % normalized by prebleach Intensity(I_t0)
        % double normalized by the intensity ratio in the reference region
        % (refI_t0) divdided by the ref region ROI at time ti refI_ti
        normalized_i = (I_ti./I_t0).*(refI_t0/refI_ti); 

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('refI_ti') = refI_ti;
        data(t).('ref_ratio') = (refI_t0/refI_ti);
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

% function outputs 
struct_out = alldata;
end 