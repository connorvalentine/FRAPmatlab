%% script to analyze the data structures

% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - fix fun_data_average so the error propagation is solid

% - Go through the dataset and remove/fix outliers (45C etc)
% - go through the movies, find way to flag an experiment to ignore it.

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
    
%% Initialize Data Folder, structures, and sample ID data
    global boxfolder outputfolder 
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,'z_outputs','main_analysis');
    
    av = struct(); % averaged data structure
    all = struct(); % all data we want to plot maybe
    save_im = 'y';


%% importing the D0 data for BSA in free solution
id_BSA = struct();
data_BSA = struct();
avBSA = struct(); % averaged data structure
allBSA = struct(); % all data we want to plot maybe
    
bsa_folder = fullfile(boxfolder,'z_BSA_in_Free_Solution');
cd(bsa_folder)
    
    id_BSA.T25 = load('z_BSA_in_Free_Solution_25C_id').id;
    data_BSA.T25 = load('z_BSA_in_Free_Solution_25C_pd').pd;
    id_BSA.T35 = load('z_BSA_in_Free_Solution_35C_id').id;
    data_BSA.T35 = load('z_BSA_in_Free_Solution_35C_pd').pd;
    id_BSA.T45 = load('z_BSA_in_Free_Solution_45C_id').id;
    data_BSA.T45 = load('z_BSA_in_Free_Solution_45C_pd').pd;
    id_BSA.T55 = load('z_BSA_in_Free_Solution_55C_id').id;
    data_BSA.T55 = load('z_BSA_in_Free_Solution_55C_pd').pd;
cd(mainfolder)

[avBSA,allBSA] = fun_data_average(data_BSA,id_BSA,avBSA);

%% importing the 25C data 
% temporary structures
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_BSA_25C_trial_6_id').id;
    data.F87 = load('F87_BSA_25C_trial_6_pd').pd;
    id.F127 = load('F127_BSA_25C_trial_6_id').id;
    data.F127 = load('F127_BSA_25C_trial_6_pd').pd;
    id.P123 = load('P123_BSA_25C_trial_6_id').id;
    data.P123 = load('P123_BSA_25C_trial_6_pd').pd;    
    cd(mainfolder)

[av,all] = fun_data_average(data,id,av);
%% importing the 35C data 
% temporary structures
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_BSA_35C_trial_6_id').id;
    data.F87 = load('F87_BSA_35C_trial_6_pd').pd;
    id.F127 = load('F127_BSA_35C_trial_6_id').id;
    data.F127 = load('F127_BSA_35C_trial_6_pd').pd;
    id.P123 = load('P123_BSA_35C_trial_6_id').id;
    data.P123 = load('P123_BSA_35C_trial_6_pd').pd;    
    cd(mainfolder)

[av,all] = fun_data_average(data,id,av);
%% importing the 45C data 
% temporary structures
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_BSA_45C_trial_6_id').id;
    data.F87 = load('F87_BSA_45C_trial_6_pd').pd;
    id.F127 = load('F127_BSA_45C_trial_6_id').id;
    data.F127 = load('F127_BSA_45C_trial_6_pd').pd;
    id.P123 = load('P123_BSA_45C_trial_6_id').id;
    data.P123 = load('P123_BSA_45C_trial_6_pd').pd;    
    cd(mainfolder)

[av,all] = fun_data_average(data,id,av);

%% importing the 55C data 
% temporary structures
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_BSA_55C_trial_6_id').id;
    data.F87 = load('F87_BSA_55C_trial_6_pd').pd;
    id.F127 = load('F127_BSA_55C_trial_6_id').id;
    data.F127 = load('F127_BSA_55C_trial_6_pd').pd;
    id.P123 = load('P123_BSA_55C_trial_6_id').id;
    data.P123 = load('P123_BSA_55C_trial_6_pd').pd;    
    cd(mainfolder)

[av,all] = fun_data_average(data,id,av);

%% Normalizing the data for use down below
%normalizing the data by the bsa in solution D0 data. Using Light
%scattering data for D0 in solution not the FRAP
D_mean = [56.0;75.7;94.8;99.8];
D_err = [0.75;0.80;0.40;3.18];
k = 0;
for topfield = fieldnames(av)'
    k = k+1;
    temperature_field = topfield{1}
    for field = fieldnames(id)'
        pluronic_field = field{1}
        BSA_fields = fieldnames(avBSA.(temperature_field));
        BSA_field = BSA_fields{1};
        
        D_temp = av.(temperature_field).(pluronic_field).D;
        Dneg_temp = av.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = av.(temperature_field).(pluronic_field).Dpos;
        
%         D0 = avBSA.(temperature_field).(BSA_field).D
%         Dneg0 = avBSA.(temperature_field).(BSA_field).Dneg;
%         Dpos0 = avBSA.(temperature_field).(BSA_field).Dpos;
        D0 = D_mean(k);
        Dneg0 = D_err(k);
        Dpos0 = D_err(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        av.(temperature_field).(pluronic_field).D =  D_normalized;
        av.(temperature_field).(pluronic_field).Dneg = D_normalized_neg;
        av.(temperature_field).(pluronic_field).Dpos =  D_normalized_pos;
    end
end

%% exporting to Excel
c = av.all25C.F87.c
% T = table(LastName,Age,Weight,Smoker);
% T(1:5,:)

for topfield = fieldnames(av)'
    k = k+1;
    temperature_field = topfield{1};
    file_name = [temperature_field,'.xlsx'];
    for field = fieldnames(id)'
        pluronic_field = field{1}
        concentration = av.(temperature_field).(pluronic_field).c';
        D_D0 = av.(temperature_field).(pluronic_field).D'; %normalized 
        D_D0_negative_error = av.(temperature_field).(pluronic_field).Dneg';
        D_D0_positive_error = av.(temperature_field).(pluronic_field).Dpos';
        fm0 = av.(temperature_field).(pluronic_field).fm0';
        
        T = table(concentration,D_D0,D_D0_negative_error,D_D0_positive_error,fm0);
        writetable(T,file_name,'Sheet',pluronic_field)
    end
end