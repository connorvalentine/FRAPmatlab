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
    for field = fieldnames(id)'
        pluronic_field = field{1};
        
        D_temp = BSA.(temperature_field).(pluronic_field).D;
        Dneg_temp = BSA.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = BSA.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_BSA(k);
        Dneg0 = D_err_BSA(k);
        Dpos0 = D_err_BSA(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        BSA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
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
    for field = fieldnames(id)'
        pluronic_field = field{1};
        
        D_temp = LYS.(temperature_field).(pluronic_field).D;
        Dneg_temp = LYS.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = LYS.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_LYS(k);
        Dneg0 = D_err_LYS(k);
        Dpos0 = D_err_LYS(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        LYS.(temperature_field).(pluronic_field).rh_prot = rh_prot';
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
    for field = fieldnames(id)'
        pluronic_field = field{1};
        
        D_temp = CHA.(temperature_field).(pluronic_field).D;
        Dneg_temp = CHA.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = CHA.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_CHA(k);
        Dneg0 = D_err_CHA(k);
        Dpos0 = D_err_CHA(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        CHA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
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
    for field = fieldnames(id)'
        pluronic_field = field{1};
        
        D_temp = HSA.(temperature_field).(pluronic_field).D;
        Dneg_temp = HSA.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = HSA.(temperature_field).(pluronic_field).Dpos;
        
        D0 = D_mean_HSA(k);
        Dneg0 = D_err_HSA(k);
        Dpos0 = D_err_HSA(k);

        D_normalized = 55.1* D_temp/D0; %% bug, its already normalized by 55.1 when it comes out of the data analysis
        
        D_normalized_neg =  D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos =  D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        % adding the hydrodynamic radius of the protein to this data
        rh_prot = zeros(length(D_normalized),1);
        rh_prot = rh_prot + rh;
        
        % adding the temperature to the matrix (is dumb workaround)
        Tmatrix = zeros(length(D_normalized),1);
        T = Tmatrix+str2num(temperature_field(4:5));
        
        HSA.(temperature_field).(pluronic_field).T = T';
        HSA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
        HSA.(temperature_field).(pluronic_field).rh_prot = rh_prot';
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

% 3. D vs T (prot,C,plur)
pluronics = struct();
pluronics.('F87') = struct();
pluronics.('F127') = struct();
pluronics.('P123') = struct();

mainT = struct();
mainT.('BSA') = pluronics;
mainT.('LYS') = pluronics;
mainT.('CHA') =pluronics;
mainT.('HSA') = pluronics



%
clc
proteinCount = 0;
numFigures = 0;

for proteins = fieldnames(main)'
    protein_name = proteins{1};
    protein_data = main.(proteins{1});
    i = 0;
   
   
    for temperatures = fieldnames(protein_data)' %all25C, all35C etc
        i = i+1;
        field = temperatures{1};
        temperature = str2num(field(4:5));      
        
        for k = 1:length(fieldnames(protein_data.(field))')
            pluronics = fieldnames(protein_data.(field))';
            pluronic = pluronics{k};
            data = protein_data.(field).(pluronic);
            
            for v = 1:length(data.c)
                if temperature == 25
                    concentration = ['wt',num2str(round(data.c(v),0))];
                    mainT.(protein_name).(pluronic).(concentration) = struct();
                    mainT.(protein_name).(pluronic).(concentration).('T') = temperature;
                    mainT.(protein_name).(pluronic).(concentration).('D') = data.D(v);
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
%% Plotting all D v C BSA data on separate  pluronic plots

% % cheng comparison data
% F127c = [7;10;12;15;17.5;20;23;25;27];
% F127d = [0.2;0.1;0.08;0.04;0.025;0.012;0.009;0.006;0.003];
% 
% F87c = [15;20;25;30;37;40;43];
% F87d = [0.06;0.035;0.017;0.005;0.0017;0.001;0.0004];
% 
% P123c = [10;15;20;30;33;35];
% P123d = [0.18;0.12;0.07;0.004;0.0004;0.00025];
close all

% normal plots
proteinCount = 0;
numFigures = 0;
for proteins = fieldnames(main)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1}
    protein_data = main.(proteins{1});
    
    all_colors = parula(69);

    colors(1,:) = all_colors(5,:);
    colors(2,:) = all_colors(30,:);
    colors(3,:) = all_colors(45,:);
    colors(4,:) = all_colors(60,:);
    i = 0;
   
    for temperatures = fieldnames(protein_data)' %all25C, all35C etc
        i = i+1;
        field = temperatures{1}
        temperature = str2num(field(4:5));
        color = colors(i,:);
        
        
        for k = 1:length(fieldnames(protein_data.(field))')
            
            pluronics = fieldnames(protein_data.(field))';
            pluronic = pluronics{k};
            temp = protein_data.(field).(pluronic);

            fig = figure(k+numFigures);
            set(k+numFigures, 'WindowStyle', 'Docked');

            p1 = errorbar(temp.c,temp.D,temp.Dneg,temp.Dpos);
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

            if i ==4
                lgd = legend('25C','35C','45C','55C');
                lgd.Orientation = 'vertical';
                lgd.NumColumns = 1;
                lgd.Location = 'eastoutside';
                lgd.Title.String = protein_name;
            else
            end        

            if save_im == 'y' && temperature == 55
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 6.5 4];             % define location to save the images 
    %             plot_name = [pluronic,'DvC_withChengData'];
                plot_name = [pluronic,'DvC_',protein_name];
                plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
                print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
            else
            end
        end

    end
     numFigures = numFigures + length(fieldnames(protein_data.(field)));
end
% 


%% Plotting all fm0vC data on separate  pluronic plot 
close all

% normal plots
proteinCount = 0;
numFigures = 0;
for proteins = fieldnames(main)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1}
    protein_data = main.(proteins{1});
    
    all_colors = parula(69);

    colors(1,:) = all_colors(5,:);
    colors(2,:) = all_colors(30,:);
    colors(3,:) = all_colors(45,:);
    colors(4,:) = all_colors(60,:);
    i = 0;
   
    for temperatures = fieldnames(protein_data)' %all25C, all35C etc
        i = i+1;
        field = temperatures{1}
        temperature = str2num(field(4:5));
        color = colors(i,:);
        
        for k = 1:length(fieldnames(protein_data.(field))')
            
            pluronics = fieldnames(protein_data.(field))';
            pluronic = pluronics{k};
            temp = protein_data.(field).(pluronic);

            fig = figure(k+numFigures);
            set(k+numFigures, 'WindowStyle', 'Docked');

            p1 = plot(temp.c,temp.fm0);
            hold on
            p1.Marker = 'o';
            p1.Color = color;
            p1.MarkerFaceColor = color;
            p1.LineStyle = 'none';
                
            ax = gca;
            if  strcmp(pluronic, 'F127')
                ax.XLim = [17 31];
            elseif strcmp(pluronic, 'F87')
                ax.XLim = [24 44];
            elseif strcmp(pluronic, 'P123')
                 ax.XLim = [19 36];
            else
            end
            ax.YScale = 'linear';
            ax.YLim = [0 1.05]; 
            ax.XLabel.String = [pluronic ,' wt%'];
            ax.YLabel.String = 'fraction of mobile proteins';
 
            lgd = legend('25C','35C','45C','55C');
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'southwest';
            lgd.Title.String = protein_name;
            lgd.LineWidth = 0.5;  

            if save_im == 'y' && temperature == 55
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 4.5 4];             % define location to save the images 
                plot_name = [pluronic,'fm0vC_',protein_name];
                plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
                print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
            else
            end

        end

    end
     numFigures = numFigures + length(fieldnames(protein_data.(field)));
end

%% Plotting D vs T for each concentration separate plots
close all


% normal plots
proteinCount = 0;
numFigures = 0;
for proteins = fieldnames(main)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1}
    protein_data = main.(proteins{1});
    
    colors = parula(6);
    all_colors = parula(69);

    colors(6,:) = all_colors(1,:);
    colors(5,:) = all_colors(20,:);
    colors(4,:) = all_colors(35,:);
    colors(3,:) = all_colors(45,:);
    colors(2,:) = all_colors(55,:);
    colors(1,:) = all_colors(69,:);
    i = 0;
   
    for temperatures = fieldnames(protein_data)' %all25C, all35C etc
        i = i+1;
        field = temperatures{1};
        temperature = str2num(field(4:5));
        
        k = 0;
        for pluronics = fieldnames(protein_data.(field))'
            k = k+1;
            pluronic = pluronics{1};
            temp = protein_data.(field).(pluronic);

            fig = figure(k+numFigures);
            set(k+numFigures, 'WindowStyle', 'Docked');

            for j = 1:length(temp.D)
                if length(temp.D) > 3
                    color = colors(j,:);
                else
                    color = colors(j+3,:);
                end
                p1 = errorbar(temperature,temp.D(j),temp.Dneg(j),temp.Dpos(j));
                hold on
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';
                
                ax = gca;
                ax.YScale = 'log';
                ax.XLim = [20 60];
                ax.YLim = [1e-4 1];
                ax.XLabel.String = 'Temperature \circC';
                ax.YLabel.String = [protein_name,' D/D_0'];
                yticks([ 1e-4 1e-3 1e-2 1e-1 1e0])
            end
            
            if strcmp(protein_name,'BSA')
                if  strcmp(pluronic, 'F127')
                     lgd = legend('17.5','20','22.5','25','27.5','30');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F127 wt%';
                elseif strcmp(pluronic, 'F87')
                     lgd = legend('25','30','35','37.5','40','42.5');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F87 wt%';
                elseif strcmp(pluronic, 'P123')
                     lgd = legend('20 ','25 ','27.5','30','32.5','35');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'P123 wt%';              
                else
                end
            else
                if  strcmp(pluronic, 'F127')
                     lgd = legend('25','27.5','30');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F127 wt%';
                elseif strcmp(pluronic, 'F87')
                     lgd = legend('37.5','40','42.5');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F87 wt%';        
                else
                end
            end

            if save_im == 'y' && temperature == 55
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 6 4];             % define location to save the images 
                plot_name = [pluronic,'DvT_',protein_name];
                plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
                print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
            else
            end

        end

    end
     numFigures = numFigures + length(fieldnames(protein_data.(field)));
end

%% Plotting fm0 vs T for each concentration separate plots
close all

% normal plots
proteinCount = 0;
numFigures = 0;
for proteins = fieldnames(main)'
    proteinCount = proteinCount+1;
    protein_name = proteins{1}
    protein_data = main.(proteins{1});
    
    colors = parula(6);
    all_colors = parula(69);

    colors(6,:) = all_colors(1,:);
    colors(5,:) = all_colors(20,:);
    colors(4,:) = all_colors(35,:);
    colors(3,:) = all_colors(45,:);
    colors(2,:) = all_colors(55,:);
    colors(1,:) = all_colors(69,:);
    i = 0;
   
    for temperatures = fieldnames(protein_data)' %all25C, all35C etc
        i = i+1;
        field = temperatures{1};
        temperature = str2num(field(4:5));
        color = colors(i,:);
        k = 0;
        for pluronics = fieldnames(protein_data.(field))'
            k = k+1;
            pluronic = pluronics{1};
            temp = protein_data.(field).(pluronic);

            fig = figure(k+numFigures);
            set(k+numFigures, 'WindowStyle', 'Docked');

            for j = 1:length(temp.D)  
                if length(temp.D) > 3
                    color = colors(j,:);
                else
                    color = colors(j+3,:);
                end
                p1 = plot(temperature,temp.fm0(j));
                hold on
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';
                
                ax = gca;
                ax.YScale = 'linear';
                ax.XLim = [20 60];
                ax.YLim = [0 1.1];
                ax.XLabel.String = 'Temperature \circC';
                ax.YLabel.String = [protein_name,' fm_0'];
            end
            
            if strcmp(protein_name,'BSA')
                if  strcmp(pluronic, 'F127')
                     lgd = legend('17.5','20','22.5','25','27.5','30');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F127 wt%';
                elseif strcmp(pluronic, 'F87')
                     lgd = legend('25','30','35','37.5','40','42.5');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F87 wt%';
                elseif strcmp(pluronic, 'P123')
                     lgd = legend('20 ','25 ','27.5','30','32.5','35');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'P123 wt%';              
                else
                end
            else
                if  strcmp(pluronic, 'F127')
                     lgd = legend('25','27.5','30');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F127 wt%';
                elseif strcmp(pluronic, 'F87')
                     lgd = legend('37.5','40','42.5');
                     lgd.Location = 'eastoutside';
                     lgd.Title.String = 'F87 wt%';        
                else
                end
            end

            if save_im == 'y' && temperature == 55
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 6 4];             % define location to save the images 
                plot_name = [pluronic,'fm0vT_',protein_name];
                plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
                print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
            else
            end

        end

    end
     numFigures = numFigures + length(fieldnames(protein_data.(field)));
end

%% Plotting D vs hydrodynamic radius of protein for each temperature
close all

% normal plots
proteinCount = 0;
numFigures = 0;

   
colors = parula(6);
all_colors = parula(69);

colors(6,:) = all_colors(1,:);
colors(5,:) = all_colors(20,:);
colors(4,:) = all_colors(35,:);
colors(3,:) = all_colors(45,:);
colors(2,:) = all_colors(55,:);
colors(1,:) = all_colors(69,:);
i = 0;

for temperatures = fieldnames(main.BSA)' %all25C, all35C etc
    i = i+1;
    field = temperatures{1};
    temperature = str2num(field(4:5))
    color = colors(i,:);
    k = 0;

    for pluronics = fieldnames(main.BSA.(field))'
        k = k+1;
        pluronic = pluronics{1};
        
        if strcmp(pluronic,'P123')
            continue
        else  
        end
        fig = figure(k+numFigures);
        set(k+numFigures, 'WindowStyle', 'Docked');

        for proteins = fieldnames(main)'
            
            proteinCount = proteinCount+1;
            protein_name = proteins{1};
            protein_data = main.(proteins{1});
            temp = protein_data.(field).(pluronic);
            
            for j = 1:length(temp.D)  
                if length(temp.D) > 3
                    color = colors(j,:);
                else
                    color = colors(j+3,:);
                end
                p1 = errorbar(temp.rh_prot(j),temp.D(j),temp.Dneg(j),temp.Dpos(j));
                hold on
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';

                ax = gca;
                ax.YScale = 'log';
                ax.XLim = [1.5 6];
                ax.YLim = [1e-4 1];
                ax.XLabel.String = 'R_h [nm]';
                ax.YLabel.String = ['D/D_0 at ', field(4:5),'\circC'];
                yticks([ 1e-4 1e-3 1e-2 1e-1 1e0])
            end
            
            if  strcmp(pluronic, 'F127')
                 lgd = legend('17.5','20','22.5','25','27.5','30');
                 lgd.Location = 'eastoutside';
                 lgd.Title.String = 'F127 wt%';
            elseif strcmp(pluronic, 'F87')
                 lgd = legend('25','30','35','37.5','40','42.5');
                 lgd.Location = 'eastoutside';
                 lgd.Title.String = 'F87 wt%';
            elseif strcmp(pluronic, 'P123')
                 lgd = legend('20 ','25 ','27.5','30','32.5','35');
                 lgd.Location = 'eastoutside';
                 lgd.Title.String = 'P123 wt%';              
            else
            end
        end  
        if save_im == 'y'
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 6 4];             % define location to save the images 
            plot_name = [pluronic,'DvRh_',field(4:5),'C'];
            plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
            print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
        else
        end

    end
numFigures = numFigures + length(fieldnames(protein_data.(field)));
end
     
%% Plotting fm0 vs hydrodynamic radius of protein for each temperature
close all

% normal plots
proteinCount = 0;
numFigures = 0;

   
colors = parula(6);
all_colors = parula(69);

colors(6,:) = all_colors(1,:);
colors(5,:) = all_colors(20,:);
colors(4,:) = all_colors(35,:);
colors(3,:) = all_colors(45,:);
colors(2,:) = all_colors(55,:);
colors(1,:) = all_colors(69,:);
i = 0;

for temperatures = fieldnames(main.BSA)' %all25C, all35C etc
    i = i+1;
    field = temperatures{1};
    temperature = str2num(field(4:5))
    color = colors(i,:);
    k = 0;

    for pluronics = fieldnames(main.BSA.(field))'
        k = k+1;
        pluronic = pluronics{1};
        
        if strcmp(pluronic,'P123')
            continue
        else  
        end
        fig = figure(k+numFigures);
        set(k+numFigures, 'WindowStyle', 'Docked');

        for proteins = fieldnames(main)'
            
            proteinCount = proteinCount+1;
            protein_name = proteins{1};
            protein_data = main.(proteins{1});
            temp = protein_data.(field).(pluronic);
            
            for j = 1:length(temp.D)  
                if length(temp.D) > 3
                    color = colors(j,:);
                else
                    color = colors(j+3,:);
                end
                p1 = plot(temp.rh_prot(j),temp.fm0(j));
                hold on
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';

                ax = gca;
                ax.YScale = 'linear';
                ax.XLim = [1 5];
                ax.YLim = [0 1.1];
                ax.XLabel.String = 'R_h [nm]';
                ax.YLabel.String = ['fm_0 at ', field(4:5),'\circC'];
            end
            
            if  strcmp(pluronic, 'F127')
                 lgd = legend('17.5','20','22.5','25','27.5','30');
                 lgd.Location = 'eastoutside';
                 lgd.Title.String = 'F127 wt%';
            elseif strcmp(pluronic, 'F87')
                 lgd = legend('25','30','35','37.5','40','42.5');
                 lgd.Location = 'eastoutside';
                 lgd.Title.String = 'F87 wt%';
            elseif strcmp(pluronic, 'P123')
                 lgd = legend('20 ','25 ','27.5','30','32.5','35');
                 lgd.Location = 'eastoutside';
                 lgd.Title.String = 'P123 wt%';              
            else
            end
        end  
        if save_im == 'y' && temperature == 55
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 6 4];             % define location to save the images 
            plot_name = [pluronic,'fm0vRh_',field(4:5),'C'];
            plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
            print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
        else
        end

    end
numFigures = numFigures + length(fieldnames(protein_data.(field)));
end
% %% Plotting ln(D) vs 1/T 
% close all
% colors = parula(7);
% for i = 1:length(fieldnames(BSA))
%     temperatures = fieldnames(BSA);
%     field = temperatures{i};
%     flag = str2num(field(4:5))
%     temperature = 1./(str2num(field(4:5)) + 273)
%     k = 0;
%     
%     for pluronics = fieldnames(BSA.(field))'
%         k = k+1;
%         fig = figure(2*k);
%         set(2*k, 'WindowStyle', 'Docked');
%         hold on
%         pluronic = pluronics{1};
%         temp = BSA.(field).(pluronic);
%         
%         figure(2*k);
%         for j = 1:6
%             hold on
%             set(gca,'YScale','linear');
%             plot(temperature,log(temp.D(j)),'d','color',colors(j,:),'markerfacecolor',colors(j,:))
%             axis([0.003 0.0034 -9 0])
%             xlabel('1/Temperature [1/k]')
%             ylabel('ln(D/D_0)')
% %             yticks([ 1e-4 1e-3 1e-2 1e-1 1e0])
%         end
%         if k == 1
%             lgd = legend('25','30','35','37.5','40','42.5','location','eastoutside');
%             lgd.Title.String = 'F87 wt%';
%         elseif k == 2
%             lgd = legend('17.5','20','22.5','25','27.5','30','location','eastoutside');
%             lgd.Title.String = 'F127 wt%';
%         else
%             lgd = legend('20 ','25 ','27.5','30','32.5','35','location','eastoutside');
%             lgd.Title.String = 'P123 wt%';
%         end
%         
%         if save_im == 'y' && flag == 55
%             fig.PaperUnits = 'inches';
%             fig.PaperPosition = [0 0 6 4];             % define location to save the images 
%             plot_name = [pluronic,'lnDv1_T'];
%             plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
%             print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
%         else
%         end
%     end
%     
% end
