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
    
    % Make a list of each image file name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
   
    % loads the pre-bleach image
    f1 = fullfile(data(1).folder,data(1).name);
    im0 = imread(f1);    
    
    % add info to fits.
    fits.(position).('pbim') = im0;
end  

% list the function outputs
struct_out = fits;
end 