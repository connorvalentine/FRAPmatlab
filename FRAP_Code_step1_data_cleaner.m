% FRAP parameter generator 
% Connor Valentine

% this script is used to list the parameters required to run FRAP_Code_V3.
% The aim of including this extra step is so we just need to run each
% experiment analysis by choosing the folder where the data is.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%% Part 1: Initialize some basic parameters and defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Arial')
    set(0,'DefaultAxesFontSize',9.5)
% Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial','DefaultTextFontSize',9.5)
% default axis settings for better plots
    set(groot,'defaultaxeslinewidth',1)
    set(groot,'DefaultLineLineWidth',1)
    set(groot, 'DefaultAxesBox', 'on')
    set(groot, 'DefaultAxesYGrid', 'off')
    set(groot, 'DefaultAxesXGrid', 'off')
% Reset the environment
    clear all;
    clear global all;
    clc
    close all;
% Establish the main directory, 
    global mainfolder 
    mainfolder = cd;
    
%% define the folder/experiment we are going to make a parameter structure for.
    folder1 = 'F127_HDF_25C';
    folder2 ='trial_2';

%% Specify the parameters you want to save that do not change. Experimental info
% 
%     pluronic = 'F127';
%     conc = [20;20;20;25;25;25;30;30;30];
%     conc = [17.5;20;22.5;25;27.5;30];
%     conc = [17.5;17.5;17.5;20;20;20;22.5;22.5;22.5;25;25;25;27.5;27.5;27.5;30;30;30]; 

    pluronic = 'F127';
    conc = [25;25;25];
%     pluronic = 'F127';
%     conc = [25;25;25;27.5;27.5;27.5;30;30;30];
%     
    prot = 'HDFL'; % must be 3 letter str
    protc = '1p0'; %must be 3 letter str
    temperature = '25C'; % needs to be string in this format. can pull num out later if needed
   

%% pixel size calibrations are prone to changing by other users... damn
 magnification = '20X';
if magnification == '10X'
    pixel_size_manual = 0.645; % 0.645 um/pixel for 10x 
elseif magnification == '20X'
    pixel_size_manual = 0.323; % 0.323 um/pixel for 20x
else 
    disp('error in inputs');
end
%% manually input the parameters here that will be saved into a structure.
% this way to run the experiment in the future you just need to specify the
% folders you want to run. (folder 1 and folder 2)
%% Initialize Data Folder, structures, and sample ID data
    global boxfolder 
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    id = struct();    

% get the subfolder names from the folders instead of naming by number;
    d = dir(fullfile(boxfolder,folder1,folder2,'frap'));
% remove all files from the list to just get folders (isdir property is 0)
    dfolders = d([d(:).isdir]);
% remove '.' and '..'from the list
	dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    pos_struct = rmfield(dfolders,{'bytes','date','isdir','datenum','folder'});  % cleaning up the data structure
    npos = length(pos_struct); % number of positions
    
% re-order positions so they're arranged from lowest number to highest
% instead of alphabetical
positions = zeros(npos,1);
for p = 1:npos
    A = regexp(pos_struct(p).name,'\d*','Match');
    positions(p) = str2double(A{1});
end 
% resort the position order
    [sortedpos,sortorder] = sort(positions);

% write position data structure with ID info. This will be saved as a .mat
% file inside folder2.

for p = 1:npos
    % First we make the folder name for the position we want
    % add the positions in numerical order:
    idx = sortorder(p);
    position = pos_struct(idx).name;
    id.(position).('plur') = pluronic;
    id.(position).('plwt') = conc(p);
    id.(position).('prot') = prot;
    id.(position).('protc') = protc;
    id.(position).('temp') = temperature;
    id.(position).('magnification') = magnification;
    id.(position).('pixel_size_manual') = pixel_size_manual;
end
%% saving ID structure as a .mat file we can load into other script
datafolder = fullfile(boxfolder,folder1,folder2); % full path where struct will be saved
cd(datafolder); % change over to this folder to save here 

id_name = [folder1,'_',folder2,'_','info','.mat'];
save(id_name,'id');
cd(mainfolder);

%% say ur done
disp('finished')