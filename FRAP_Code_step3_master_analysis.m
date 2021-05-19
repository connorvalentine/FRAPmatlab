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
    save_im = 'n';
%% cheng comparison data
F127c = [7;10;12;15;17.5;20;23;25;27];
F127d = [0.2;0.1;0.08;0.04;0.025;0.012;0.009;0.006;0.003];

F87c = [15;20;25;30;37;40;43];
F87d = [0.06;0.035;0.017;0.005;0.0017;0.001;0.0004];

P123c = [10;15;20;30;33;35];
P123d = [0.18;0.12;0.07;0.004;0.0004;0.00025];

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

%% Plotting all the data at each temperature in subplots
%normalizing the data by the bsa in solution D0 data
for topfield = fieldnames(av)'
    temperature_field = topfield{1}
    for field = fieldnames(id)'
        pluronic_field = field{1}
        BSA_fields = fieldnames(avBSA.(temperature_field));
        BSA_field = BSA_fields{1};
        
        D_temp = av.(temperature_field).(pluronic_field).D;
        Dneg_temp = av.(temperature_field).(pluronic_field).Dneg;
        Dpos_temp = av.(temperature_field).(pluronic_field).Dpos;
        
        D0 = avBSA.(temperature_field).(BSA_field).D
        Dneg0 = avBSA.(temperature_field).(BSA_field).Dneg;
        Dpos0 = avBSA.(temperature_field).(BSA_field).Dpos;
    
        D_normalized = D_temp/D0;
        D_normalized_neg = D_normalized.*sqrt((Dneg_temp./D_temp).^2 + (Dneg0./D0).^2);
        D_normalized_pos = D_normalized.*sqrt((Dpos_temp./D_temp).^2 + (Dpos0./D0).^2);
        
        av.(temperature_field).(pluronic_field).D =  D_normalized;
        av.(temperature_field).(pluronic_field).Dneg = D_normalized_neg;
        av.(temperature_field).(pluronic_field).Dpos =  D_normalized_pos;
    end
end

%% Plotting all temp data on one pluronic plot
colors = parula(69);
close all
fig = figure(123);
set(123, 'WindowStyle', 'Docked');

for i = 1:length(fieldnames(av))
    temperatures = fieldnames(av);
    field = temperatures{i};
    temperature = str2num(field(4:5));
    color = colors(i*15,:);
    for k = 1:3
        pluronics = fieldnames(av.(field))';
        pluronic = pluronics{k};
        temp = av.(field).(pluronic);
        
        subplot(2,3,k)
            
            
            errorbar(temp.c,temp.D,temp.Dneg,temp.Dpos,'d','color',color,'markerfacecolor',color)
            hold on
            set(gca,'YScale','log');
            axis([0.9*min(temp.c) 1.1*max(temp.c) 1e-5 1])
            ylabel('D/D_0 [\mum^2s^{-1}]')
            yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

        subplot(2,3,k+3)  
            
            plot(temp.c,temp.fm0,'o','color',color,'markerfacecolor','none','linewidth',3)
            hold on
            axis([0.9*min(temp.c) 1.1*max(temp.c) 0 1 ])
            xlabel([pluronic,' wt%'])
            yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
            ylabel('f_m')

    end
    subplot(2,3,2)
    legend('D: 25C','D: 35C','D: 45C','D: 55C','location','north','orientation','horizontal','NumColumns', 2)
    subplot(2,3,5)
    legend('D: 25C','D: 35C','D: 45C','D: 55C','location','south','orientation','horizontal','NumColumns', 2)
    box on
end

if save_im == 'y'
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 14 6];             % define location to save the images 
    plot_name = ['All Pluronics','DvC'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end

%% Plotting all fm0 data on one plot (without fit
fig = figure(12345);
set(12345, 'WindowStyle', 'Docked');

for i = 1:length(fieldnames(av))
    temperatures = fieldnames(av); % temperature of the dataset
    field = temperatures{i};
    temperature = str2num(field(4:5));
    color = colors(i*15,:);
    k = 0;
    for pluronics = fieldnames(av.(field))'
        k = k+1;
        pluronic = pluronics{1};
        temp = av.(field).(pluronic);
        
        subplot(1,3,k)
            hold on
            plot(temp.c,temp.fm0,'d','color',color,'markerfacecolor',color)
            axis([0.9*min(temp.c) 1.1*max(temp.c) 0.1 1.2 ])
            xlabel([pluronic,' wt%'])
            ylabel('fraction of mobile proteins')
            if i ==4
            legend('25C','35C','45C','55C','location','best','NumColumns',2)
            else
            end
    end
end

if save_im == 'y'
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 3];             % define location to save the images 
    plot_name = ['All Pluronics','Fm0vC'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end
%% Plotting D vs temp for each concentration
for i = 1:length(fieldnames(av))
    temperatures = fieldnames(av);
    field = temperatures{i};
    temperature = str2num(field(4:5));
    color = colors(i*15,:);
    k = 0;
    
    for pluronics = fieldnames(av.(field))'
        k = k+1;
        pluronic = pluronics{1};
        temp = av.(field).(pluronic);
        
        if i == 1
            fig = figure(k);
            set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float 
        else 
        end
        
        figure(k);
            for j = 1:6
                subplot(2,3,j)
                    hold on
                    set(gca,'YScale','log');
                    errorbar(temperature,temp.D(j),temp.Dneg(j),temp.Dpos(j),'d','color',color,'markerfacecolor',color)
                    axis([20 60 1e-5 1])
                    xlabel([pluronic,' ',num2str(temp.c(j)),' wt%'])
                    ylabel('D/D_0')
                    yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
            end
    end
    
end

if save_im == 'y'
    for k = 1:3
        pluronics = fieldnames(av.(field))';
        pluronic = pluronics{k};
        fig = figure(k);
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 10 10];             % define location to save the images 
        plot_name = [pluronic,'_DvT'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    end
else
end

%% Plotting fm vs temp for each concentration

for i = 1:length(fieldnames(av))
    temperatures = fieldnames(av);
    field = temperatures{i};
    temperature = str2num(field(4:5));
    color = colors(i*15,:);
    k = 0;
    
    for pluronics = fieldnames(av.(field))'
        k = k+1;
        pluronic = pluronics{1};
        temp = av.(field).(pluronic);
        
        if i == 1
            fig = figure(k);
            set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float 
        else 
        end
        
        figure(k*45);
        set(k*45, 'WindowStyle', 'Docked');
            for j = 1:6
                subplot(2,3,j)
                    hold on
                    set(gca,'YScale','linear');
                    plot(temperature,temp.fm0(j),'d','color',color,'markerfacecolor',color)
                    axis([20 60 0 1.2])
                    xlabel([pluronic,' ',num2str(temp.c(j)),' wt%'])
                    ylabel('f_m')
            end
    end
    
end

if save_im == 'y'
for k = 1:3
    pluronics = fieldnames(av.(field))';
    pluronic = pluronics{k};
    fig = figure(2*45);
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10 10];             % define location to save the images 
    plot_name = [pluronic,'_fmvT'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
end
else
end
%% old
% subplot(1,3,1)
%     set(gca,'YScale','log');
%     hold all
%     errorbar(pd.c,pd.D,pd.Dneg,pd.Dpos,'db')
%     axis([0.9*min(pd.c) 1.1*max(pd.c) 0.7*min(pd.D) 2*max(pd.D)])
%     xlabel([id.(position).plur,' wt%'])
%     ylabel('D/D_0 [\mum^2s^{-1}]')
% 
% subplot(1,3,2)
%     hold on
%     errorbar(pd.c,pd.r,pd.rlb,pd.rub,'or')
%     axis([0.9*min(pd.c) 1.1*max(pd.c) 0.8*min(pd.r) 1.2*max(pd.r)])
%     xlabel([id.(position).plur,' wt%'])
%     ylabel('laser radius [um]')
%     
% subplot(1,3,3)
%     hold on
% %     plot(pd.c,pd.fm,'om')
%     errorbar(pd.c,pd.fm,pd.fmlb,pd.fmub,'om')
%     xlabel([id.(position).plur,' wt%'])
%     ylabel('mobile fraction [fm]')
%     axis([0.9*min(pd.c) 1.1*max(pd.c) 0.5 1.2])
% 
% 
% 
% %% plotting    
% close all
% fig = figure('name','All','visible','on');
% set(fig, 'WindowStyle', 'Docked'); 
% 
% 
% % need for loop to unpack
% k = 0;
% lgd = [];
% syms = ['d','s','o'];
% for topfield = fieldnames(id)'
%     pluronic = topfield{1};
%     C = parula(3*length(fieldnames(id)));
%     i = 0;
%     k = k+1;
%     conc = [];
%     C1(:,k) = C(2+2*k,:);
%     for field = fieldnames(id.(pluronic))'
%         position = field{1};
%         fits2 = fits.(pluronic);
%         id2 = id.(pluronic);
%         i = i+1;
%         D = fits2.(position).D/55.1;
%         c = id2.(position).plwt;
%         pxl = fits2.(position).pixel_size;
%         r = fits2.(position).radius*pxl;
%         rall(i) =r;
%         conc(i) = c;
%         Dall(i) = D; % from cheng to normalize by free soln BSA D
%         Dneg = fits2.(position).errD(1)/55.1;
%         Dpos = fits2.(position).errD(2)/55.1;
%         Dnegall(i) = Dneg;
%         Dposall(i) = Dpos;
%     end
%     lgd{k} = pluronic;
%     hold on 
%     errorbar(conc,Dall,Dnegall,Dposall,'linestyle','none','Color',C1(:,k),'MarkerFaceColor',C1(:,k),'Marker',syms(k),'MarkerSize',10)
% end
% lgd{4} = 'F87 Cheng';
% lgd{5} = 'F127 Cheng';
% lgd{6} = 'P123 Cheng';
% 
% plot(F87c,F87d,'linestyle','none','Color',C1(:,1),'Marker',syms(1),'MarkerSize',10,'linewidth',2)
% plot(F127c,F127d,'linestyle','none','Color',C1(:,2),'Marker',syms(2),'MarkerSize',10,'linewidth',2)
% plot(P123c,P123d,'linestyle','none','Color',C1(:,3),'Marker',syms(3),'MarkerSize',10,'linewidth',2)
% 
% axis([5 45 1e-4 1e2])
% set(gca,'YScale','log');
% xlabel(['wt%'])
% ylabel('D/D_0 [\mum^2s^{-1}]')
% title('Diffusivity of BSA in Pluronic at 25\circC')
% legend(lgd)
% 
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 8 6];             % define location to save the images 
% plot_name = ['Cheng_vs_CSV_DvC_25C'];
% plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
% print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
