%% data loader for main and protein analysis

%% initialize function
function [main, mainT, mainSAXS,DLS] = data_loader()
    global mainfolder boxfolder outputfolder plotfolder
%% Initialize Data Folder, structures, and sample ID data

    
    
    BSA = struct(); % averaged data structure
    allBSA = struct(); % all data we want to plot 
    LYS = struct();
    allLYS = struct();
    CHA = struct();
    allCHA = struct();   
    HSA = struct();
    allHSA = struct();
    
  
    
%% protein diffusivity data from light scattering
% Normalizing the data for use down below
%normalizing the data by the bsa in solution D0 data. Using Light
%scattering data for D0 in solution not the FRAP
dls_path = 'C:\Users\user\Box\Sorted Data FRAP\APS data';
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

DLS = struct();
DLS.('BSA').('rh') = rh_BSA;
DLS.('BSA').('rh_err') = rh_err_BSA;
DLS.('BSA').('D') = D_mean_BSA;
DLS.('BSA').('D_err') = D_err_BSA;

DLS.('LYS').('rh') = rh_LYS;
DLS.('LYS').('rh_err') = rh_err_LYS;
DLS.('LYS').('D') = D_mean_LYS;
DLS.('LYS').('D_err') = D_err_LYS;

DLS.('CHA').('rh') = rh_CHA;
DLS.('CHA').('rh_err') = rh_err_CHA;
DLS.('CHA').('D') = D_mean_CHA;
DLS.('CHA').('D_err') = D_err_CHA;

DLS.('HSA').('rh') = rh_HSA;
DLS.('HSA').('rh_err') = rh_err_HSA;
DLS.('HSA').('D') = D_mean_HSA;
DLS.('HSA').('D_err') = D_err_HSA;

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
        HSA.(temperature_field).(pluronic_field).D_raw = D_raw;
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

%% crystal size data from APS SAXS
% This makes a better data structure for plotting things vs Temperature
mainSAXS = mainT;

SAXS_path = 'C:\Users\user\Box\Sorted Data FRAP\APS data';
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
                mainSAXS.(protein_name).(pluronic);
                continue
            else
            end
            
            data = mainT.(protein_name).(pluronic).(c_string);
            a = zeros(1,length(data.T));
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
                    a(1,i) = test.a;
                    rm(1,i) = test.rm;
                    ro(1,i) = test.r_oct;
                    conc(1,i) = c;
                    k = test.extrapolated{1};
                    extrapolated(1,i) = k;
                    if c == 43
                        conc(1,i) = 42.5; % cant have decimals in name
                    elseif c == 38
                        conc(1,i) = 37.5; % cant have decimals in name
                    end
                    
                    if c == 18
                        conc(1,i) = 17.5; % cant have decimals in name
                    elseif c == 23
                        conc(1,i) = 22.5; % cant have decimals in name
                    elseif c == 28
                        conc(1,i) = 27.5; % cant have decimals in name
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
            mainSAXS.(protein_name).(pluronic).(c_string).('a') = a;
            mainSAXS.(protein_name).(pluronic).(c_string).('rm') = rm;
            mainSAXS.(protein_name).(pluronic).(c_string).('ro') = ro;
            mainSAXS.(protein_name).(pluronic).(c_string).('extrapolated') = extrapolated;
            
        end
    end
end
end
