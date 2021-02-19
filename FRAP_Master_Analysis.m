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
    plotfolder = fullfile(boxfolder,'z_outputs','main_analysis');
    
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
F127d = [0.2;0.1;0.08;0.04;0.025;0.012;0.009;0.006;0.003];

F87c = [15;20;25;30;37;40;43];
F87d = [0.06;0.035;0.017;0.005;0.0017;0.001;0.0004];

P123c = [10;15;20;30;33;35];
P123d = [0.18;0.12;0.07;0.004;0.0004;0.00025];

%% plotting    
close all
fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 


% need for loop to unpack
k = 0;
lgd = [];
syms = ['d','s','o'];
for topfield = fieldnames(id)'
    pluronic = topfield{1};
    C = parula(3*length(fieldnames(id)));
    i = 0;
    k = k+1;
    conc = [];
    C1(:,k) = C(2+2*k,:);
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
    errorbar(conc,Dall,Dnegall,Dposall,'linestyle','none','Color',C1(:,k),'MarkerFaceColor',C1(:,k),'Marker',syms(k),'MarkerSize',10)
end
lgd{4} = 'F87 Cheng';
lgd{5} = 'F127 Cheng';
lgd{6} = 'P123 Cheng';

plot(F87c,F87d,'linestyle','none','Color',C1(:,1),'Marker',syms(1),'MarkerSize',10,'linewidth',2)
plot(F127c,F127d,'linestyle','none','Color',C1(:,2),'Marker',syms(2),'MarkerSize',10,'linewidth',2)
plot(P123c,P123d,'linestyle','none','Color',C1(:,3),'Marker',syms(3),'MarkerSize',10,'linewidth',2)

axis([5 45 1e-4 1e2])
set(gca,'YScale','log');
xlabel(['wt%'])
ylabel('D/D_0 [\mum^2s^{-1}]')
title('Diffusivity of BSA in Pluronic at 25\circC')
legend(lgd)

fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 6];             % define location to save the images 
plot_name = ['Cheng_vs_CSV_DvC_25C'];
plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
