%% Function description: 
% readins in the prebleach images and adds information to data structures
%% inputs and outputs
% inputs: id structure and fits structure. 

function [struct_out] = fun_prebleach_frame(id,fits)
% define the global variables
global boxfolder folder1 folder2
folder3 = 'prebleach';

for field = fieldnames(id)'
    position = field{1};
    
    % pull parameters from fits (made by fun_first_bleached_frame)
    idx_ref = fits.(position).idx_ref;
    idx_bleach = fits.(position).idx_bleach;
    radius = fits.(position).radius;
    center = fits.(position).center; 
    
    % Make a list of each image file name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
   
    % loads the pre-bleach image
    f1 = fullfile(data(1).folder,data(1).name);
    im0 = double(imread(f1));    

    % calculate intensity of the reference region before bleaching
    refI_t0 = mean(im0(idx_ref));
    
    % calculate intensity of the spot before bleaching 
    I_t0 = mean(im0(idx_bleach));
    
    % add info to fits.
    fits.(position).('refI_t0') = refI_t0;
    fits.(position).('I_t0') = I_t0;
end  

% list the function outputs
struct_out = fits;
end 