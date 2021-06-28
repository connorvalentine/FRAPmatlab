%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%% Part 1: Initialize some basic parameters and defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Arial')
    set(0,'DefaultAxesFontSize',14)
% Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial','DefaultTextFontSize',14)
% default axis settings for better plots
    set(groot,'defaultaxeslinewidth',1)
    set(groot,'DefaultLineLineWidth',1)
    set(groot, 'DefaultAxesBox', 'on')
    set(groot, 'DefaultAxesYGrid', 'off')
    set(groot, 'DefaultAxesXGrid', 'off') 
% turn the error message beeps off
    beep off
% Reset the environment
    clear all;
    clear global all;
    clc
    close all;
% Establish the main directory, 
    global mainfolder 
    mainfolder = cd; 
    
%% %%%%%%%%%%%%%%%% Inputs Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Everything that must be selected when code is running smoothly 
% Part 1: choose the folders
% NOTE:         prebleach folder must be named 'prebleach'
%               frap folder must be named 'frap'

global folder1 folder2
    folder1 = 'F87_BSA_25C';
    folder2 = 'trial_6';

% Part 2: Would you like to save the plots generated? 
% Note:         Select 'y' or 'n'
save_im = 'y'; 

% Part 3: Are you troubleshooting the fits? 
% Note:         Select 'y' or 'n'
trouble = 'n'; 

% do you want to make movies? (yes is first time analyzing)
global makemovies npoints
makemovies = 'y';

% number of points per wt%. is usually 3 now.
npoints = 3;
%% %%%%%%%%%%%%%%%% Parameters Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Data Folder, structures, and sample ID data
global boxfolder outputfolder plotfolder datafolder moviefolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    datafolder = fullfile(boxfolder,folder1,folder2); % full path where outputs structures will be saved
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,'z_outputs','z_plots'); % plot folder within the data folder for the defined experiment
    moviefolder = fullfile(boxfolder,'z_outputs','z_movies'); % plot folder within the data folder for the defined experiment
%initialize the structures to store all data in 
alldata = struct();
fits = struct();
images = struct();

% load in the info structure. This strucutre is made by FRAP_data_cleaner
id_name = [folder1,'_',folder2,'_','info','.mat'];
temp_struct = load(fullfile(datafolder,id_name));
id = temp_struct.id;

%% Reading in the prebleached image, and the first image after bleaching
% radius is found using improfile and fitting a gaussian curve to it
tic; % timing 

[id,fits,fig,stats] = fun_radius_finder_0(id,fits,save_im);

disp(["First Circles found in " + string(round(toc)) + " s."])
%% Reading in the rest of the frap data 
tic;
close all
[alldata] = fun_frap_frames_drift_tracking(id,fits,alldata);

disp(["Data loaded in " + string(round(toc)) + " s."])
