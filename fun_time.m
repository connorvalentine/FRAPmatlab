%% function purpose
% reads in the metadata from micromanager, and outputs the time, and pixel
% size for the specific image that you input

%% inputs
% metadata_filename =  the full filename of the metadata .txt file, which
%                        includes the path
% frame_filename = the .tif name of the picture that we are trying to find
%                        metadata for (without the full file path)

%% outputs
% e_time is elapsed time since the start of the experiment in micromanager
% r_time is the calendar time that the image was taken 
% pixel_size is the size of the pixel in microns 

% PIXEL SIZE ONLY WORKS IF YOU RAN A PIXEL SIZE CALIBRATION AT THE SAME
% % MAGNIFICATION BEFORE YOU STARTED THE EXPERIMENT
% % Make a list of each image name in our position folder
% global boxfolder folder1 folder2 
% folder3 = 'frap';
% counter = 0;
% position = 'Pos9';
% list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
% data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
% n_images = length(list);
%    
% for t = 1:n_images 
%     [data] = fun_time(data, t);
% end

function [data] = fun_time(data, t)
%define the metadata filename
    txt = dir(fullfile(data(1).folder,'*.txt'));% lists all files with .txt ending
    metadata_filename = fullfile(txt.folder,txt.name);
    
% define some other variables 
    frame_filename = data(t).name;
    
% read in the metadata text 
    metadata = fileread(metadata_filename); 
    idx = strfind(metadata,frame_filename); % indices where the picture filename is found in metadata

% the filename is mentioned 3 times. We know that the time information
% is in between the 2nd and 3rd mentions
    lb = idx(1,1);
    ub = idx(1,3);
    info = metadata(lb:ub); % isolate the portion of metadata where we know the time is.

% find the index where "Recieved Time" is mentioned.
    p1 = strfind(info,'"ElapsedTime-ms":');
    p2 = strfind(info,'"ReceivedTime":');
    p3 = strfind(info,'"PixelSizeUm":');
    p4 = strfind(info,'"PixelSizeAffine"');

% pull out Recieved time
    temp2 = info(p2:p3);
    lb2 = strfind(temp2,':');
    ub2 = strfind(temp2,'-');
    r_time = temp2(lb2(1)+3:ub2(3)-2)% keep as string

% pull out pixel size
    temp3 = info(p3:p4);
    lb3 = strfind(temp3,':');
    ub3 = strfind(temp3,',');
    pixel_size = temp3(lb3(1)+2:ub3(1)-1); 
    pixel_size = str2num(pixel_size); % convert to a number
    
% function outputs
    % add the metadata time information to our data structure
    data(t).('r_time') = r_time; % recieved time according to the computer clock 
    data(t).('pixel_size') = pixel_size; % microns and assumes pixel-size calibration was run
end
