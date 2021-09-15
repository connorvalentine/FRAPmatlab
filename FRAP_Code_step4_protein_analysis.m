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
    
%% pick your colors 
all_colors = parula(69);
colors6(1,:) = all_colors(60,:);
colors6(2,:) = all_colors(50,:);
colors6(3,:) = all_colors(40,:);
colors6(4,:) = all_colors(30,:);
colors6(5,:) = all_colors(20,:);
colors6(6,:) = all_colors(10,:);

all_colors = parula(69);
colors4(1,:) = all_colors(50,:);
colors4(2,:) = all_colors(40,:);
colors4(3,:) = all_colors(20,:);
colors4(4,:) = all_colors(1,:);

all_colors = parula(69);
colors3(1,:) = all_colors(10,:);
colors3(2,:) = all_colors(30,:);
colors3(3,:) = all_colors(50,:);

%% Initialize Data Folder, structures, and sample ID data
    global boxfolder outputfolder 
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,'z_outputs','main_analysis');
    
    BSA = struct(); % averaged data structure
    allBSA = struct(); % all data we want to plot 
    LYS = struct();
    allLYS = struct();
    CHA = struct();
    allCHA = struct();   
    HSA = struct();
    allHSA = struct();
    
    save_im = 'y';

%% protein diffusivity data from light scattering
% Normalizing the data for use down below
%normalizing the data by the bsa in solution D0 data. Using Light
%scattering data for D0 in solution not the FRAP
dls_path = 'C:\Users\user\Desktop\mathub\APS 2021';
filename = 'Protein DLS.xlsx';
dls_data = readtable(fullfile(dls_path,filename),'Range','A1:F17');

% define ranges of the table
bsa = 1:4;
lys = 5:8;
cha = 9:12;
hsa = 13:16;

rh_BSA = dls_data.rh(bsa);
rh_err_BSA = dls_data.stdRh(bsa);
D_mean_BSA = dls_data.D(bsa);
D_err_BSA = dls_data.stdD(bsa);

rh_LYS = dls_data.rh(lys);
rh_err_LYS = dls_data.stdRh(lys);
D_mean_LYS = dls_data.D(lys);
D_err_LYS = dls_data.stdD(lys);

rh_CHA = dls_data.rh(cha);
rh_err_CHA = dls_data.stdRh(cha);
D_mean_CHA = dls_data.D(cha);
D_err_CHA = dls_data.stdD(cha);

rh_HSA = dls_data.rh(hsa);
rh_err_HSA = dls_data.stdRh(hsa);
D_mean_HSA = dls_data.D(hsa);
D_err_HSA = dls_data.stdD(hsa);



%% importing the BSA data
% importing the 25C data 
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
[BSA,allBSA] = fun_data_average(data,id,BSA);

%importing the 35C data 
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
[BSA,allBSA] = fun_data_average(data,id,BSA);

% importing the 45C data 
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
[BSA,allBSA] = fun_data_average(data,id,BSA);

% importing the 55C data 
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
[BSA,allBSA] = fun_data_average(data,id,BSA);


k = 0;
for topfield = fieldnames(BSA)'
    k = k+1;
    temperature_field = topfield{1};
    rh = rh_BSA(k); %hydrodynamic radius at this temperature from DLS
    rh_std = rh_err_BSA(k);
    for temperature_name = fieldnames(id)'
        pluronic_field = temperature_name{1};
        
        D_temp = BSA.(temperature_field).(pluronic_field).D;
        Dneg_temp = BSA.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = BSA.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_BSA(k);
        Dneg0 = D_err_BSA(k);
        Dpos0 = D_err_BSA(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        D_raw = 55.1* D_temp;
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        % adding the std of hydrodynamic radius of the protein to this data
        rh_prot_std = zeros(length(D_normalized),1);
        rh_prot_std = rh_prot_std + rh_std;
        
        BSA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
        BSA.(temperature_field).(pluronic_field).rh_prot_std = rh_prot_std';
        BSA.(temperature_field).(pluronic_field).D_raw =  D_raw;
        BSA.(temperature_field).(pluronic_field).D =  D_normalized;
        BSA.(temperature_field).(pluronic_field).Dneg = D_normalized_neg;
        BSA.(temperature_field).(pluronic_field).Dpos =  D_normalized_pos;
    end
end
%% importing the LYS data
% importing the 25C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_LYS_25C_trial_1_id').id;
    data.F87 = load('F87_LYS_25C_trial_1_pd').pd;
    id.F127 = load('F127_LYS_25C_trial_1_id').id;
    data.F127 = load('F127_LYS_25C_trial_1_pd').pd;

    cd(mainfolder)
[LYS,allLYS] = fun_data_average(data,id,LYS);

%importing the 35C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_LYS_35C_trial_1_id').id;
    data.F87 = load('F87_LYS_35C_trial_1_pd').pd;
    id.F127 = load('F127_LYS_35C_trial_1_id').id;
    data.F127 = load('F127_LYS_35C_trial_1_pd').pd;

    cd(mainfolder)
[LYS,allLYS] = fun_data_average(data,id,LYS);

% importing the 45C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_LYS_45C_trial_1_id').id;
    data.F87 = load('F87_LYS_45C_trial_1_pd').pd;
    id.F127 = load('F127_LYS_45C_trial_1_id').id;
    data.F127 = load('F127_LYS_45C_trial_1_pd').pd;
 
    cd(mainfolder)
[LYS,allLYS] = fun_data_average(data,id,LYS);

% importing the 55C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_LYS_55C_trial_1_id').id;
    data.F87 = load('F87_LYS_55C_trial_1_pd').pd;
    id.F127 = load('F127_LYS_55C_trial_1_id').id;
    data.F127 = load('F127_LYS_55C_trial_1_pd').pd;
   
    cd(mainfolder)
[LYS,allLYS] = fun_data_average(data,id,LYS);

% Normalizing the data for use down below
%normalizing the data by the bsa in solution D0 data. Using Light
%scattering data for D0 in solution not the FRAP

k = 0;
for topfield = fieldnames(LYS)'
    k = k+1;
    temperature_field = topfield{1};
    rh = rh_LYS(k); %hydrodynamic radius at this temperature from DLS
    rh_std = rh_err_LYS(k);
    
    for temperature_name = fieldnames(id)'
        pluronic_field = temperature_name{1};
        
        D_temp = LYS.(temperature_field).(pluronic_field).D;
        Dneg_temp = LYS.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = LYS.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_LYS(k);
        Dneg0 = D_err_LYS(k);
        Dpos0 = D_err_LYS(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        D_raw = 55.1* D_temp;
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        % adding the std of hydrodynamic radius of the protein to this data
        rh_prot_std = zeros(length(D_normalized),1);
        rh_prot_std = rh_prot_std + rh_std;
        
        LYS.(temperature_field).(pluronic_field).rh_prot = rh_prot';
        LYS.(temperature_field).(pluronic_field).rh_prot_std = rh_prot_std';
         LYS.(temperature_field).(pluronic_field).D_raw =  D_raw;
        LYS.(temperature_field).(pluronic_field).D =  D_normalized;
        LYS.(temperature_field).(pluronic_field).Dneg = D_normalized_neg;
        LYS.(temperature_field).(pluronic_field).Dpos =  D_normalized_pos;
    end
end
%% importing the CHA data
% importing the 25C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_CHA_25C_trial_1_id').id;
    data.F87 = load('F87_CHA_25C_trial_1_pd').pd;
    id.F127 = load('F127_CHA_25C_trial_1_id').id;
    data.F127 = load('F127_CHA_25C_trial_1_pd').pd;

    cd(mainfolder)
[CHA,allCHA] = fun_data_average(data,id,CHA);

%importing the 35C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_CHA_35C_trial_1_id').id;
    data.F87 = load('F87_CHA_35C_trial_1_pd').pd;
    id.F127 = load('F127_CHA_35C_trial_1_id').id;
    data.F127 = load('F127_CHA_35C_trial_1_pd').pd;

    cd(mainfolder)
[CHA,allCHA] = fun_data_average(data,id,CHA);

% importing the 45C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_CHA_45C_trial_1_id').id;
    data.F87 = load('F87_CHA_45C_trial_1_pd').pd;
    id.F127 = load('F127_CHA_45C_trial_1_id').id;
    data.F127 = load('F127_CHA_45C_trial_1_pd').pd;
 
    cd(mainfolder)
[CHA,allCHA] = fun_data_average(data,id,CHA);

% importing the 55C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_CHA_55C_trial_1_id').id;
    data.F87 = load('F87_CHA_55C_trial_1_pd').pd;
    id.F127 = load('F127_CHA_55C_trial_1_id').id;
    data.F127 = load('F127_CHA_55C_trial_1_pd').pd;
   
    cd(mainfolder)
[CHA,allCHA] = fun_data_average(data,id,CHA);

% Normalizing the data for use down below
%normalizing the data by the bsa in solution D0 data. Using Light
%scattering data for D0 in solution not the FRAP

k = 0;
for topfield = fieldnames(CHA)'
    k = k+1;
    temperature_field = topfield{1};
    rh = rh_CHA(k); %hydrodynamic radius at this temperature from DLS
    rh_std = rh_err_CHA(k);
    for temperature_name = fieldnames(id)'
        pluronic_field = temperature_name{1};
        
        D_temp = CHA.(temperature_field).(pluronic_field).D;
        Dneg_temp = CHA.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = CHA.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_CHA(k);
        Dneg0 = D_err_CHA(k);
        Dpos0 = D_err_CHA(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        D_raw = 55.1* D_temp;
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        % adding the std of hydrodynamic radius of the protein to this data
        rh_prot_std = zeros(length(D_normalized),1);
        rh_prot_std = rh_prot_std + rh_std;
        
        CHA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
        CHA.(temperature_field).(pluronic_field).rh_prot_std = rh_prot_std';
        CHA.(temperature_field).(pluronic_field).D_raw =  D_raw;
        CHA.(temperature_field).(pluronic_field).D =  D_normalized;
        CHA.(temperature_field).(pluronic_field).Dneg = D_normalized_neg;
        CHA.(temperature_field).(pluronic_field).Dpos =  D_normalized_pos;
    end
end
%% importing the HSA data
% importing the 25C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_HSA_25C_trial_1_id').id;
    data.F87 = load('F87_HSA_25C_trial_1_pd').pd;
    id.F127 = load('F127_HSA_25C_trial_1_id').id;
    data.F127 = load('F127_HSA_25C_trial_1_pd').pd;

    cd(mainfolder)
[HSA,allHSA] = fun_data_average(data,id,HSA);

%importing the 35C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_HSA_35C_trial_1_id').id;
    data.F87 = load('F87_HSA_35C_trial_1_pd').pd;
    id.F127 = load('F127_HSA_35C_trial_1_id').id;
    data.F127 = load('F127_HSA_35C_trial_1_pd').pd;

    cd(mainfolder)
[HSA,allHSA] = fun_data_average(data,id,HSA);

% importing the 45C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_HSA_45C_trial_1_id').id;
    data.F87 = load('F87_HSA_45C_trial_1_pd').pd;
    id.F127 = load('F127_HSA_45C_trial_1_id').id;
    data.F127 = load('F127_HSA_45C_trial_1_pd').pd;
 
    cd(mainfolder)
[HSA,allHSA] = fun_data_average(data,id,HSA);

% importing the 55C data 
    id = struct();
    data = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_HSA_55C_trial_1_id').id;
    data.F87 = load('F87_HSA_55C_trial_1_pd').pd;
    id.F127 = load('F127_HSA_55C_trial_1_id').id;
    data.F127 = load('F127_HSA_55C_trial_1_pd').pd;
   
    cd(mainfolder)
[HSA,allHSA] = fun_data_average(data,id,HSA);

% Normalizing the data for use down below
%normalizing the data by the bsa in solution D0 data. Using Light
%scattering data for D0 in solution not the FRAP

k = 0;
for topfield = fieldnames(HSA)'
    k = k+1;
    temperature_field = topfield{1};
    rh = rh_HSA(k); %hydrodynamic radius at this temperature from DLS
    rh_std = rh_err_HSA(k);
    for temperature_name = fieldnames(id)'
        pluronic_field = temperature_name{1};
        
        D_temp = HSA.(temperature_field).(pluronic_field).D;
        Dneg_temp = HSA.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = HSA.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_HSA(k);
        Dneg0 = D_err_HSA(k);
        Dpos0 = D_err_HSA(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        D_raw = 55.1* D_temp;
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        % adding the std of hydrodynamic radius of the protein to this data
        rh_prot_std = zeros(length(D_normalized),1);
        rh_prot_std = rh_prot_std + rh_std;
        
        
        % adding the temperature to the matrix (is dumb workaround)
        Tmatrix = zeros(length(D_normalized),1);
        T = Tmatrix+str2num(temperature_field(4:5));
        
        HSA.(temperature_field).(pluronic_field).T = T';
        HSA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
        HSA.(temperature_field).(pluronic_field).rh_prot_std =rh_prot_std';
        HSA.(temperature_field).(pluronic_field).D_raw = D_raw
        HSA.(temperature_field).(pluronic_field).D =  D_normalized;
        HSA.(temperature_field).(pluronic_field).Dneg = D_normalized_neg;
        HSA.(temperature_field).(pluronic_field).Dpos =  D_normalized_pos;
    end
end

%% creating the main data structure to use for plots below
main = struct();
main.BSA = BSA;
main.LYS = LYS;
main.CHA = CHA;
main.HSA = HSA;


%% Making better data structures for plotting vs tempe
% data we need (structure layers)
% 1. D vs C (prot,T,Plur)         DONE, see main
% 2. fm vs C (prot,T,Plur)        DONE, see main

% This makes a better data structure for plotting things vs Temperature
pluronics = struct();
pluronics.('F87') = struct();
pluronics.('F127') = struct();
pluronics.('P123') = struct();

pluronics2 = struct();
pluronics2.('F87') = struct();
pluronics2.('F127') = struct();

mainT = struct();
mainT.('BSA') = pluronics;
mainT.('LYS') = pluronics2;
mainT.('CHA') = pluronics2;
mainT.('HSA') = pluronics2;

for proteins = fieldnames(main)'
    protein_name = proteins{1};
    i = 0;   
    for temperatures = fieldnames(main.(protein_name))' %all25C, all35C etc
        i = i+1;
        temperature_name = temperatures{1};
        temperature = str2num(temperature_name(4:5));      
        
        for k = 1:length(fieldnames(main.(protein_name).(temperature_name))')
            pluronics = fieldnames(main.(protein_name).(temperature_name))';
            pluronic = pluronics{k};
            data = main.(protein_name).(temperature_name).(pluronic);
            
            for v = 1:length(data.c)
                if temperature == 25
                    concentration = ['wt',num2str(round(data.c(v),0))];
                    mainT.(protein_name).(pluronic).(concentration) = struct();
                    mainT.(protein_name).(pluronic).(concentration).('T') = temperature;
                    mainT.(protein_name).(pluronic).(concentration).('D') = data.D(v);
                    mainT.(protein_name).(pluronic).(concentration).('D_raw') = data.D_raw(v);
                    mainT.(protein_name).(pluronic).(concentration).('Dneg') = data.Dneg(v);
                    mainT.(protein_name).(pluronic).(concentration).('Dpos') = data.Dpos(v);
                    mainT.(protein_name).(pluronic).(concentration).('fm0') = data.fm0(v);
                    mainT.(protein_name).(pluronic).(concentration).('fmlb') = data.fmlb(v);
                    mainT.(protein_name).(pluronic).(concentration).('fmub') = data.fmub(v);
                    mainT.(protein_name).(pluronic).(concentration).('rh_prot') = data.rh_prot(v);
                else 
                    concentration = ['wt',num2str(round(data.c(v),0))];
                    mainT.(protein_name).(pluronic).(concentration).T(i) = temperature;
                    mainT.(protein_name).(pluronic).(concentration).D(i) = data.D(v);
                    mainT.(protein_name).(pluronic).(concentration).D_raw(i) = data.D_raw(v);
                    mainT.(protein_name).(pluronic).(concentration).Dneg(i) = data.Dneg(v);
                    mainT.(protein_name).(pluronic).(concentration).Dpos(i)  = data.Dpos(v);
                    mainT.(protein_name).(pluronic).(concentration).fm0(i)  = data.fm0(v);
                    mainT.(protein_name).(pluronic).(concentration).fmlb(i)  = data.fmlb(v);
                    mainT.(protein_name).(pluronic).(concentration).fmub(i)  = data.fmub(v);
                    mainT.(protein_name).(pluronic).(concentration).rh_prot(i)  = data.rh_prot(v);   
                end
            end
        end
    end

end

%% adding a structure for legend names bcuz cant start strucutre names with a number
lgd_names = struct();
pluronics = struct();
pluronics.('F87') = struct();
pluronics.('F127') = struct();
pluronics.('P123') = struct();

pluronics2 = struct();
pluronics2.('F87') = struct();
pluronics2.('F127') = struct();

lgd_names.('BSA') = pluronics;
lgd_names.('LYS') = pluronics2;
lgd_names.('CHA') = pluronics2;
lgd_names.('HSA') = pluronics2;

for proteins = fieldnames(lgd_names)'
    protein_name = proteins{1};
    
    for pluronics = fieldnames(lgd_names.(protein_name))' 
        pluronic = pluronics{1};
        
        if strcmp(protein_name, 'BSA')
            if strcmp(pluronic,'F87') 
                lgd_names.(protein_name).(pluronic) = {'25 wt%';'30 wt%';'35 wt%';'37.5 wt%';'40 wt%';'42.5 wt%'};
            elseif strcmp(pluronic,'F127') 
                lgd_names.(protein_name).(pluronic) = {'17.5 wt%';'20 wt%';'22.5 wt%';'25 wt%';'27.5 wt%';'30 wt%'};
            elseif strcmp(pluronic,'P123') 
                lgd_names.(protein_name).(pluronic) = {'20 wt%';'25 wt%';'27.5 wt%';'30 wt%';'32.5 wt%';'35 wt%'};
            else
            end
        else % other proteins have only 3 concentrations and all the same
            if strcmp(pluronic,'F87') 
                lgd_names.(protein_name).(pluronic) = {'37.5 wt%';'40 wt%';'42.5 wt%'};
            elseif strcmp(pluronic,'F127') 
                lgd_names.(protein_name).(pluronic) = {'25 wt%';'27.5 wt%';'30 wt%'};
            else
            end
        end
    end
end

lgd_names.('all') = {'37.5wt% F87';'40wt% F87';'42.5wt% F87';'25wt% F127';'27.5wt% F127';'30wt% F127'};

% legend names when only sampels with saxs data
lgd_names_saxs = struct();
pluronics = struct();
pluronics.('F87') = struct();
pluronics.('F127') = struct();
pluronics.('P123') = struct();

pluronics2 = struct();
pluronics2.('F87') = struct();
pluronics2.('F127') = struct();

lgd_names_saxs.('BSA') = pluronics;
lgd_names_saxs.('LYS') = pluronics2;
lgd_names_saxs.('CHA') = pluronics2;
lgd_names_saxs.('HSA') = pluronics2;

for proteins = fieldnames(lgd_names_saxs)'
    protein_name = proteins{1};
    
    for pluronics = fieldnames(lgd_names_saxs.(protein_name))' 
        pluronic = pluronics{1};
        
        if strcmp(protein_name, 'BSA')
            if strcmp(pluronic,'F87') 
                lgd_names_saxs.(protein_name).(pluronic) = {'37.5 wt%';'40 wt%';'42.5 wt%'};
            elseif strcmp(pluronic,'F127') 
                lgd_names_saxs.(protein_name).(pluronic) = {'22.5 wt%';'30 wt%'};
            elseif strcmp(pluronic,'P123') 
                lgd_names_saxs.(protein_name).(pluronic) = {'35 wt%'};
            else
            end
        else % other proteins have only 3 concentrations and all the same
            if strcmp(pluronic,'F87') 
                lgd_names_saxs.(protein_name).(pluronic) = {'37.5 wt%';'40 wt%';'42.5 wt%'};
            elseif strcmp(pluronic,'F127') 
                lgd_names_saxs.(protein_name).(pluronic) = {'25 wt%';'27.5 wt%';'30 wt%'};
            else
            end
        end
    end
end

%% crystal size data from APS SAXS
% This makes a better data structure for plotting things vs Temperature
mainSAXS = mainT;

SAXS_path = 'C:\Users\user\Box\Presentations\0_Papers\BSA Diffusion\Geometry of FCC and BCC and calculations\SAXS 2021 data and DLS';
filename = 'SAXS_data_clean.xlsx';
SAXS_data = readtable(fullfile(SAXS_path,filename),'Range','A1:L53');

pluronic_find = ismember(SAXS_data.pluronic,'F127');
wt_find =   ismember(SAXS_data.wt,22);
temperature_find =   ismember(SAXS_data.T,25);
   % comma separated list expansion 
test = logical(pluronic_find.*wt_find.*temperature_find);
   
test = SAXS_data(test,:);
for proteins = fieldnames(mainT)'
    protein_name = proteins{1};
    for pluronics = fieldnames(mainT.(protein_name))'
        pluronic = pluronics{1};
        pluronic_find = ismember(SAXS_data.pluronic,pluronic);
        temporary = struct();
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_string = concentrations{1};
            mainSAXS.(protein_name).(pluronic).(c_string) = struct();
            c = str2num(c_string(3:4));
            wt_find =   ismember(SAXS_data.wt,c);
            
            if sum(wt_find.*pluronic_find) == 0
                mainSAXS.(protein_name).(pluronic) = rmfield(mainSAXS.(protein_name).(pluronic),c_string);
                disp(['No SAXS data', protein_name, pluronic, c_string])
                mainSAXS.(protein_name).(pluronic)
                continue
            else
            end
            
            data = mainT.(protein_name).(pluronic).(c_string);
            rm = zeros(1,length(data.T));
            ro = zeros(1,length(data.T));
            extrapolated = ["","","",""];
            for i = 1:length(data.T)
                
                temperature_find = ismember(SAXS_data.T,data.T(i));
                
                SAXS_data_row = logical(pluronic_find.*wt_find.*temperature_find);
                
                test = SAXS_data(SAXS_data_row,:);
                if isempty(test)
                    % do nothing - there is no SAXS data for this example
                else
                    rm(1,i) = test.rm;
                    ro(1,i) = test.r_oct;
                    conc(1,i) = c;
                    a = test.extrapolated{1};
                    extrapolated(1,i) = a;
                    if c == 43
                        conc(1,i) = 42.5; % cant have decimals in name
                    elseif c == 38
                        conc(1,i) = 37.5; % cant have decimals in name
                    end
                end
            end
            mainSAXS.(protein_name).(pluronic).(c_string).('T') = data.T;
            mainSAXS.(protein_name).(pluronic).(c_string).('c') = conc;
            mainSAXS.(protein_name).(pluronic).(c_string).('D') = data.D;
            mainSAXS.(protein_name).(pluronic).(c_string).('Dneg') = data.Dneg;
            mainSAXS.(protein_name).(pluronic).(c_string).('Dpos') = data.Dpos;
            mainSAXS.(protein_name).(pluronic).(c_string).('fm0') = data.fm0;
            mainSAXS.(protein_name).(pluronic).(c_string).('fmlb') = data.fmlb;
            mainSAXS.(protein_name).(pluronic).(c_string).('fmub') = data.fmub;
            mainSAXS.(protein_name).(pluronic).(c_string).('rh_prot') = data.rh_prot;
            mainSAXS.(protein_name).(pluronic).(c_string).('rm') = rm;
            mainSAXS.(protein_name).(pluronic).(c_string).('ro') = ro;
            mainSAXS.(protein_name).(pluronic).(c_string).('extrapolated') = extrapolated;
            
        end
    end
end
%% Plotting all D v C (figure == protein),(subplots == pluronic)
close all
protein_count = 0;
for proteins = fieldnames(main)'
    protein_count = protein_count +1;
    protein_name = proteins{1};
    % % cheng comparison data
    % F127c = [7;10;12;15;17.5;20;23;25;27];
    % F127d = [0.2;0.1;0.08;0.04;0.025;0.012;0.009;0.006;0.003];
    % 
    % F87c = [15;20;25;30;37;40;43];
    % F87d = [0.06;0.035;0.017;0.005;0.0017;0.001;0.0004];
    % 
    % P123c = [10;15;20;30;33;35];
    % P123d = [0.18;0.12;0.07;0.004;0.0004;0.00025];
    
    fig = figure(protein_count);
    set(fig,'Position',[2600 300 1500 500]);
    t = 0;
    for temperatures = fieldnames(main.(protein_name))' %all25C, all35C etc
        t = t+1;
        temperature_name = temperatures{1};
        temperature = str2num(temperature_name(4:5));
        color = colors4(t,:);
        
        
        for k = 1:length(fieldnames(main.(protein_name).(temperature_name))')
            pluronics = fieldnames(main.(protein_name).(temperature_name))';
            pluronic = pluronics{k};
            data = main.(protein_name).(temperature_name).(pluronic);

            subplot(1,length(fieldnames(main.(protein_name).(temperature_name))),k);
                p1 = errorbar(data.c,data.D,data.Dneg,data.Dpos);
                hold on
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';

                ax = gca;
                ax.YScale = 'log';
                if  strcmp(pluronic, 'F127')
                    ax.XLim = [17 31];
                elseif strcmp(pluronic, 'F87')
                    ax.XLim = [24 44];
                elseif strcmp(pluronic, 'P123')
                     ax.XLim = [19 36];
                else
                end
                
                ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = [pluronic ,' wt%'];
                ax.YLabel.String = 'D/D_0 [\mum^2s^{-1}]';
                yticks([1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

            if temperature == 55
                lgd = legend({'25 \circC','35 \circC','45 \circC','55 \circC'});
                lgd.Orientation = 'vertical';
                lgd.NumColumns = 2;
                lgd.Location = 'best';
                lgd.Title.String = [protein_name,' in ',pluronic];
                lgd.LineWidth = 0.5; 
            else
            end        

            if strcmp(protein_name,'BSA')
                fig.PaperPosition = [0 0 16 5]; 

            else
                fig.PaperPosition = [0 0 12 5];  
            end
            if save_im == 'y' && temperature == 55
                fig.PaperUnits = 'inches';
                plot_name = [protein_name,' All Pluronics D vs C'];
                plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
                print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
            else
            end
        end

    end
end

%% Plotting all fm0 v C (figure == protein),(subplots == pluronic)
close all
protein_count = 0;
for proteins = fieldnames(main)'
    protein_count = protein_count +1;
    protein_name = proteins{1};
    
    fig = figure(protein_count);
    set(fig,'Position',[2600 300 1500 500]);
    t = 0;
    for temperatures = fieldnames(main.(protein_name))' %all25C, all35C etc
        t = t+1;
        temperature_name = temperatures{1};
        temperature = str2num(temperature_name(4:5));
        color = colors4(t,:);
        
        
        for k = 1:length(fieldnames(main.(protein_name).(temperature_name))')
            pluronics = fieldnames(main.(protein_name).(temperature_name))';
            pluronic = pluronics{k};
            data = main.(protein_name).(temperature_name).(pluronic);

            subplot(1,length(fieldnames(main.(protein_name).(temperature_name))),k);
                p1 = errorbar(data.c,1-data.fm0,data.fmlb,data.fmub);
                hold on
                p1.Marker = 'o';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';

                ax = gca;
                ax.YScale = 'linear';
                if  strcmp(pluronic, 'F127')
                    ax.XLim = [17 31];
                elseif strcmp(pluronic, 'F87')
                    ax.XLim = [24 44];
                elseif strcmp(pluronic, 'P123')
                     ax.XLim = [19 36];
                else
                end
                
                ax.YLim = [0 1]; 
                ax.XLabel.String = [pluronic ,' wt%'];
                ax.YLabel.String = 'Immobile protein fraction';

            if temperature == 55
                lgd = legend({'25 \circC','35 \circC','45 \circC','55 \circC'});
                lgd.Orientation = 'vertical';
                lgd.NumColumns = 2;
                lgd.Location = 'best';
                lgd.Title.String = [protein_name,' at T = '];
                lgd.LineWidth = 0.5; 
            else
            end        

            if strcmp(protein_name,'BSA')
                fig.PaperPosition = [0 0 16 5]; 

            else
                fig.PaperPosition = [0 0 12 5];  
            end
            if save_im == 'y' && temperature == 55 && k == length(fieldnames(main.(protein_name).(temperature_name))')
                fig.PaperUnits = 'inches';
                plot_name = [protein_name,' All Pluronics fm0 vs C'];
                plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
                print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
            else
            end
        end

    end
end
    

%% Plotting D vs T (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all
for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    
    fig = figure(proteinCount);
    set(fig,'Position',[2600 300 1500 500]);
    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,length(fieldnames(mainT.(protein_name))),i);
        length(fieldnames(mainT.(protein_name)));
        if strcmp(protein_name, 'BSA')
            c_count = 0;
        else
            c_count = 3;
        end
            for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
                c_count = c_count+1;
                c = concentrations{1};
                data = mainT.(protein_name).(pluronic).(c);

                p1 = errorbar(data.T,data.D,data.Dneg,data.Dpos);
                hold on
                p1.Marker = 'o';
                p1.Color = colors6(c_count,:);
                p1.MarkerFaceColor = colors6(c_count,:);
                p1.LineStyle = 'none';
            end
            
                ax = gca;
                ax.XLim = [20 60]; 
                ax.YScale = 'log';
                ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = 'Temperature [\circC]';
                ax.YLabel.String = 'D/D_0';
                
                if strcmp(protein_name,'BSA')
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 2;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                else
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                end
    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 12 5];  
    end
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = [protein_name,' All Pluronics D vs T'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end

%% Plotting fm0 vs T (figure == protein),(subplots == pluronic)

proteinCount = 0;
numFigures = 0;
close all
for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    
    fig = figure(proteinCount);
    set(fig,'Position',[2600 300 1500 500]);
    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,length(fieldnames(mainT.(protein_name))),i);
        length(fieldnames(mainT.(protein_name)));
        if strcmp(protein_name, 'BSA')
            c_count = 0;
        else
            c_count = 3;
        end
            for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
                c_count = c_count+1;
                c = concentrations{1};
                data = mainT.(protein_name).(pluronic).(c);

                p1 = errorbar(data.T,data.fm0,data.fmlb,data.fmub);
                hold on
                p1.Marker = 'o';
                p1.Color = colors6(c_count,:);
                p1.MarkerFaceColor = colors6(c_count,:);
                p1.LineStyle = 'none';
            end
            
                ax = gca;
                ax.XLim = [20 60]; 
                ax.YScale = 'linear';
                ax.YLim = [1e-2 1e0]; 
                ax.XLabel.String = 'Temperature [\circC]';
                ax.YLabel.String = 'Mobile protein fraction';
                
                if strcmp(protein_name,'BSA')
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'northwest';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                    lgd.Location = 'southwest';
                else
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                    lgd.Location = 'southwest';
                end
    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 12 5];  
    end
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = [protein_name,' All Pluronics fm0 vs T'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end


%% Plotting D vs rh from DLS (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    

    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,3,i);
        if strcmp(protein_name, 'BSA')
            c_count = 0;
        else
            c_count = 3;
        end
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_count = c_count+1;
            c = concentrations{1};
            data = mainT.(protein_name).(pluronic).(c);

            p1 = errorbar(data.rh_prot,data.D,data.Dneg,data.Dpos);
            hold on
            p1.Marker = 'o';
            p1.Color = colors6(c_count,:);
            p1.MarkerFaceColor = colors6(c_count,:);
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [1 6]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1e0]; 
        ax.XLabel.String = 'R_h Protein';
        ax.YLabel.String = 'D/D_0 [\mum^2s^{-1}]';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = [protein_name,' All Pluronics D vs RhProt'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end

%% Plotting fm0 vs rh from DLS (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,3,i);
        if strcmp(protein_name, 'BSA')
            c_count = 0;
        else
            c_count = 3;
        end
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_count = c_count+1;
            c = concentrations{1};
            data = mainT.(protein_name).(pluronic).(c);

            p1 = errorbar(data.rh_prot,1-data.fm0,data.fmlb,data.fmub);
            hold on
            p1.Marker = 'o';
            p1.Color = colors6(c_count,:);
            p1.MarkerFaceColor = colors6(c_count,:);
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [1 6]; 
        ax.YScale = 'linear';
        ax.YLim = [0 1]; 
        ax.XLabel.String = 'R_h Protein';
        ax.YLabel.String = 'Immobile protein fraction';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'best';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = [protein_name,' All Pluronics fm0 vs RhProt'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end

%% Plotting D vs rh/rm from DLS (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    

    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        if i ==3
            break
        else
        end
        pluronic = pluronics{1};
        
        subplot(1,2,i);
            c_count = 0;

        if protein_name == 'BSA'
            ind = 3;
        else
            ind = 0;
        end
        for k = 1:(length(fieldnames(mainSAXS.(protein_name).(pluronic))')-ind)
            n_conc = length(fieldnames(mainSAXS.(protein_name).(pluronic))');
            concentrations = fieldnames(mainSAXS.(protein_name).(pluronic))';
            c_count = c_count+1;
            
            
            c = concentrations{1+n_conc-k};
            disp([protein_name,' ',pluronic,' ',c])
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);

            xdata = data.rh_prot./SAXSdata.rm;
            
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            p1.Marker = 'o';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [0 1]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1e0]; 
        ax.XLabel.String = 'R_h Protein / R_{micelle}';
        ax.YLabel.String = 'D/D_0 [\mum^2s^{-1}]';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(flip([lgd_names_saxs.(protein_name).(pluronic)],1));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end

end
if save_im == 'y'
    fig.PaperUnits = 'inches';
    plot_name = [' All Pluronics All Prot D vs Rh_div_Rm'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end


%% Plotting fm0 vs rh/rm from DLS (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    

    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,3,i);
            c_count = 0;

        for k = 1:length(fieldnames(mainSAXS.(protein_name).(pluronic))')
            n_conc = length(fieldnames(mainSAXS.(protein_name).(pluronic))');
            concentrations = fieldnames(mainSAXS.(protein_name).(pluronic))';
            c_count = c_count+1;
            c = concentrations{1+n_conc-k};
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.rm;
            
            p1 = errorbar(xdata,1-data.fm0,data.fmlb,data.fmub);
            hold on
            p1.Marker = 'o';
            p1.Color = colors6(c_count,:);
            p1.MarkerFaceColor = colors6(c_count,:);
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [0 1.5]; 
        ax.YScale = 'linear';
        ax.YLim = [0 1]; 
        ax.XLabel.String = 'R_h Protein / R_{micelle}';
        ax.YLabel.String = 'Immobile Fraction (fm)';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(flip([lgd_names_saxs.(protein_name).(pluronic)],1));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end

end
if save_im == 'y'
    fig.PaperUnits = 'inches';
    plot_name = [' All Pluronics All Prot fm vs Rh_div_Rm'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end
%% Plotting D vs rh/ro from DLS (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    

    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,3,i);
            c_count = 0;

        for k = 1:length(fieldnames(mainSAXS.(protein_name).(pluronic))')
            n_conc = length(fieldnames(mainSAXS.(protein_name).(pluronic))');
            concentrations = fieldnames(mainSAXS.(protein_name).(pluronic))';
            c_count = c_count+1;
            c = concentrations{1+n_conc-k};
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.ro;
            
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            p1.Marker = 'o';
            p1.Color = colors6(c_count,:);
            p1.MarkerFaceColor = colors6(c_count,:);
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [0 2]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1e0]; 
        ax.XLabel.String = 'R_h Protein / R_{octahedral}';
        ax.YLabel.String = 'D/D_0';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(flip([lgd_names_saxs.(protein_name).(pluronic)],1));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end

end
if save_im == 'y'
    fig.PaperUnits = 'inches';
    plot_name = [' All Pluronics All Prot D vs Rh_div_Ro'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end
%% Plotting D vs rh/rm from DLS (figure == ALLDATA EVERYTHING ONE PLOT (Except p123)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

all_colors = parula(69);
colors3(1,:) = all_colors(10,:);
colors3(2,:) = all_colors(30,:);
colors3(3,:) = all_colors(50,:);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    

    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        if i ==3
            break
        else
        end
        pluronic = pluronics{1};
        
%         subplot(1,2,i);
            c_count = 0;

        if protein_name == 'BSA'
            ind = 3;
        else
            ind = 0;
        end
        for k = 1:(length(fieldnames(mainSAXS.(protein_name).(pluronic))')-ind)
            n_conc = length(fieldnames(mainSAXS.(protein_name).(pluronic))');
            concentrations = fieldnames(mainSAXS.(protein_name).(pluronic))';
            c_count = c_count+1;
            
            
            c = concentrations{1+n_conc-k};
            disp([protein_name,' ',pluronic,' ',c])
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);

            xdata = data.rh_prot./SAXSdata.rm;
            
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            if strcmp(pluronic,'F87')
                p1.Marker = 'o';
            else
                p1.Marker = 's';
            end
                
                
            p1.Color = colors3(c_count,:);
            if SAXSdata.extrapolated(1) == "y"
                p1.MarkerFaceColor = 'none';
            else
                p1.MarkerFaceColor = colors3(c_count,:);
            end
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [0 1]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1]; 
        ax.XLabel.String = 'R_h Protein / R_{micelle}';
        ax.YLabel.String = 'D/D_0 [\mum^2s^{-1}]';
        ax.Title.String = 'F127 and F87';
        
        if strcmp(protein_name,'HSA')
            lgd = legend(flip([lgd_names.all],1));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = ['wt% , pluronic'];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end

end
if save_im == 'y'
    fig.PaperUnits = 'inches';
    plot_name = [' EVERYTHING D vs Rh_div_Rm'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end
%% Plotting D vs rh/ro from DLS (figure == ALLDATA EVERYTHING ONE PLOT (Except p123)
% normal plots
proteinCount = 0;
numFigures = 0;
close all    
fig = figure(1);
set(fig,'Position',[2600 300 1500 500]);

all_colors = parula(69);
colors3(1,:) = all_colors(10,:);
colors3(2,:) = all_colors(30,:);
colors3(3,:) = all_colors(50,:);

for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    

    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        if i ==3
            break
        else
        end
        pluronic = pluronics{1};
        
%         subplot(1,2,i);
            c_count = 0;

        if protein_name == 'BSA'
            ind = 3;
        else
            ind = 0;
        end
        for k = 1:(length(fieldnames(mainSAXS.(protein_name).(pluronic))')-ind)
            n_conc = length(fieldnames(mainSAXS.(protein_name).(pluronic))');
            concentrations = fieldnames(mainSAXS.(protein_name).(pluronic))';
            c_count = c_count+1;
            
            
            c = concentrations{1+n_conc-k};
            disp([protein_name,' ',pluronic,' ',c])
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);

            xdata = data.rh_prot./SAXSdata.ro;
            
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            if strcmp(pluronic,'F87')
                p1.Marker = 'o';
            else
                p1.Marker = 's';
            end
                
                
            p1.Color = colors3(c_count,:);
            if SAXSdata.extrapolated(1) == "y"
                p1.MarkerFaceColor = 'none';
            else
                p1.MarkerFaceColor = colors3(c_count,:);
            end
            p1.LineStyle = 'none';
        end
            
        ax = gca;
        ax.XLim = [0 1.5]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1]; 
        ax.XLabel.String = 'R_h Protein / R_{octahedral}';
        ax.YLabel.String = 'D/D_0';
        ax.Title.String = 'F127 and F87';
        
        if strcmp(protein_name,'HSA')
            lgd = legend(flip([lgd_names.all],1));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = ['wt% , pluronic'];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 16 5];  
    end

end
if save_im == 'y'
    fig.PaperUnits = 'inches';
    plot_name = [' EVERYTHING D vs Rh_div_Ro'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end

%% Plotting ln(D/T) vs 1/T  (figure == protein),(subplots == pluronic)
% normal plots
proteinCount = 0;
numFigures = 0;
close all
for proteins = fieldnames(mainT)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1};
    
    fig = figure(proteinCount);
    set(fig,'Position',[2600 300 1500 500]);
    i = 0;
    for pluronics = fieldnames(mainT.(protein_name))' 
        i = i+1; % number of pluronics
        pluronic = pluronics{1};
        
        subplot(1,length(fieldnames(mainT.(protein_name))),i);
        length(fieldnames(mainT.(protein_name)));
        if strcmp(protein_name, 'BSA')
            c_count = 0;
        else
            c_count = 3;
        end
            for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
                c_count = c_count+1;
                c = concentrations{1};
                data = mainT.(protein_name).(pluronic).(c);
                ydata = log(data.D_raw);
                xdata = 1./(data.T+273);

                p1 = plot(xdata,ydata);
                hold on
                p1.Marker = 'o';
                p1.Color = colors6(c_count,:);
                p1.MarkerFaceColor = colors6(c_count,:);
                p1.LineStyle = 'none';
            end
            
                ax = gca;
%                 ax.XLim = [20 60]; 
%                 ax.YScale = 'log';
%                 ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = 'Temperature [\circC]';
                ax.YLabel.String = 'ln(D)';
                
                if strcmp(protein_name,'BSA')
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 2;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                else
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                end
    end
    if strcmp(protein_name,'BSA')
        fig.PaperPosition = [0 0 16 5]; 
        
    else
        fig.PaperPosition = [0 0 12 5];  
    end
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = [protein_name,' All Pluronics lnD v lnT'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end