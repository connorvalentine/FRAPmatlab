%% script to analyze the data structures

% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - Selection of the files that should be made into plots
% - add comparison to cheng data
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
    id = struct();
    fits = struct();
% load in the info structure.
    cd(outputfolder)
    id.F87 = load('F87_BSA_25C_trial_1_id').id;
    fits.F87 = load('F87_BSA_25C_trial_1_fits').fits;
    id.F127 = load('F127_BSA_25C_trial_3_id').id;
    fits.F127 = load('F127_BSA_25C_trial_3_fits').fits;
    id.P123 = load('P123_BSA_25C_trial_1_id').id;
    fits.P123 = load('P123_BSA_25C_trial_1_fits').fits;    
    
    cd(mainfolder)
%% cheng comparison data
F127c = [7;10;12;15;17.5;20;23;25;27];
% F127d = [0.2;
    

%% plotting    
close all
fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 


% need for loop to unpack
k = 0;
lgd = [];

for topfield = fieldnames(id)'
    pluronic = topfield{1};
    C = parula(3*length(fieldnames(id)));
    i = 0;
    k = k+1;
    conc = [];
    C1 = C(2+2*k,:);
    for field = fieldnames(id.(pluronic))'
        position = field{1};
        fits2 = fits.(pluronic);
        id2 = id.(pluronic);
        i = i+1;
        D = fits2.(position).D/55.1;
        c = id2.(position).plwt;
        pxl = fits2.(position).pixel_size;
        r = fits2.(position).radius*pxl;
        rall(i) =r;
        conc(i) = c;
        Dall(i) = D; % from cheng to normalize by free soln BSA D
        Dneg = fits2.(position).errD(1)/55.1;
        Dpos = fits2.(position).errD(2)/55.1;
        Dnegall(i) = Dneg;
        Dposall(i) = Dpos;
    end
    lgd{k} = pluronic;
    hold on 
    errorbar(conc,Dall,Dnegall,Dposall,'linestyle','none','Color',C1,'MarkerFaceColor',C1,'Marker','d')
end
axis([15 45 1e-3 1e2])
set(gca,'YScale','log');
xlabel(['wt%'])
ylabel('D/D_0 [\mum^2s^{-1}]')
title('Diffusivity of BSA in Pluronic at 25\circC')
legend(lgd)

