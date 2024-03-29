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
    set(0,'DefaultAxesFontSize',16)
% Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial','DefaultTextFontSize',16)
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
    
    global boxfolder outputfolder plotfolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,'z_outputs','main_analysis');
    
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
colors3(1,:) = all_colors(40,:);
colors3(2,:) = all_colors(20,:);
colors3(3,:) = all_colors(1,:);
colors3(4,:) = all_colors(40,:);
colors3(5,:) = all_colors(20,:);
colors3(6,:) = all_colors(1,:);



save_im = 'y';

[main, mainT, mainSAXS, DLS] = data_loader();
%% plot defaults 
MarkerSize = 12;
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
%             if strcmp(pluronic,'F87') 
%                 lgd_names.(protein_name).(pluronic) = {'25 wt%';'30 wt%';'35 wt%';'37.5 wt%';'40 wt%';'42.5 wt%'};
%             elseif strcmp(pluronic,'F127') 
%                 lgd_names.(protein_name).(pluronic) = {'17.5 wt%';'20 wt%';'22.5 wt%';'25 wt%';'27.5 wt%';'30 wt%'};
%             elseif strcmp(pluronic,'P123') 
%                 lgd_names.(protein_name).(pluronic) = {'20 wt%';'25 wt%';'27.5 wt%';'30 wt%';'32.5 wt%';'35 wt%'};
%             else
%             end
            if strcmp(pluronic,'F87') 
                lgd_names.(protein_name).(pluronic) = {'37.5 wt%';'40 wt%';'42.5 wt%'};
            elseif strcmp(pluronic,'F127') 
                lgd_names.(protein_name).(pluronic) = {'25 wt%';'27.5 wt%';'30 wt%'};
            elseif strcmp(pluronic,'P123') 
                lgd_names.(protein_name).(pluronic) = {'30 wt%';'32.5 wt%';'35 wt%'};
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

%% Plotting DLS data (protein Rh vs Temperature 

fig = figure();

protein_count = 0;
temperatures = [25 35 45 55]';
for proteins = fieldnames(DLS)'
    protein_count = protein_count +1;
    protein_name = proteins{1};
    color = colors4(protein_count,:);
    
    data = DLS.(protein_name);
    if strcmp(protein_name, 'LYS')
        i = 1;
    elseif strcmp(protein_name, 'HSA')
        i = 2;
    elseif strcmp(protein_name, 'CHA')
        i = 3;
    elseif strcmp(protein_name, 'BSA')
        i = 4;
    end
    
    subplot(2,4,i)
        p1 = errorbar(temperatures,data.rh,data.rh_err/2,data.rh_err/2);
        hold on
        p1.Marker = 's';
        p1.Color = color;
        p1.MarkerFaceColor = color;
        p1.MarkerSize = MarkerSize;
        p1.LineStyle = 'none';

        ax = gca;
        ax.YScale = 'linear';

        ax.YLim = [1.5 5.5];
        ax.XLim = [20 60];
        ax.XLabel.String = ['Temperature [ \circC ]'];
        ax.YLabel.String = 'R_{Hydrodynamic} [ \mum ]';
        yticks([2,3,4,5]);
        lgd = legend(protein_name);
        lgd.Orientation = 'vertical';
        lgd.NumColumns = 2;
        lgd.Location = 'best';
        lgd.LineWidth = 0.5; 
        
    subplot(2,4,i+4)
        p1 = errorbar(temperatures,data.D,data.D_err,data.D_err);
        hold on
        p1.Marker = 's';
        p1.Color = color;
        p1.MarkerFaceColor = color;
        p1.LineStyle = 'none';
        p1.MarkerSize = MarkerSize;

        ax = gca;
        ax.YScale = 'log';
        ax.YLim = [1e1 5e2];
        ax.XLim = [20 60];
        ax.XLabel.String = ['Temperature [ \circC ]'];
        ax.YLabel.String = 'D [\mum^2s^{-1}]';
        yticks([10,100]);
        lgd = legend(protein_name);
        lgd.Orientation = 'vertical';
        lgd.NumColumns = 2;
        lgd.Location = 'best';
        lgd.LineWidth = 0.5; 
       
end

if save_im == 'y'
    fig.PaperPosition = [0 0 14 8];  
    fig.PaperUnits = 'inches';
    plot_name = ['DLS'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end

%% Plotting the SAXS results from APS
pluronics = {'F87';'F127'};
count = 0;
Temperatures = [25;35;45;55];
for p = 1:2
    pluronic = pluronics{p};
    count = count+1;
    for concs = fieldnames(mainSAXS.LYS.(pluronic))'
       
        conc = concs{1};
        data = mainSAXS.LYS.(pluronic).(conc);
        
        subplot(2,2,count)
            for k = 1:4
%                 a_calc = unit_cell_extrapolation(pluronic,Temperatures(k),data.c(k));
                p1 = plot(data.c(k),data.a(k));
                hold on
%                 p2 = plot(data.c(k),a_calc );
%                 p2.Marker = 'd';
                p1.Marker = 's';
                p1.Color = colors4(k,:);
                p1.MarkerFaceColor = colors4(k,:);
                p1.MarkerSize = MarkerSize;
                p1.LineStyle = 'none';
                
                markers(k) = p1;
            end
            
            ax = gca;
            ax.YScale = 'linear';
            if strcmp(pluronic,'F87')
                ax.XLim = [36 44];
                ax.YLim = [14 18];
            else
                ax.XLim = [24 31];
                ax.YLim = [20 24];
            end
            ax.XLabel.String = ['Concentration [ wt% ]'];
            ax.YLabel.String = 'Unit Cell Size [ nm ]';
            lgd = legend(markers, {'25 \circC','35 \circC','45 \circC','55 \circC'});
            lgd.Title.String = [pluronic];
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'north';
            lgd.LineWidth = 0.5; 

        subplot(2,2,count+2)
             for k = 1:4
                p1 = plot(data.c(k),data.rm(k));
                hold on
                p1.Marker = 's';
                p1.Color = colors4(k,:);
                p1.MarkerFaceColor = colors4(k,:);
                p1.MarkerSize = MarkerSize;
                p1.LineStyle = 'none';
                markers(k) = p1;
            end
            hold on
            p1.Marker = 's';
            p1.Color = color;
            p1.MarkerFaceColor = color;
            p1.MarkerSize = MarkerSize;
            p1.LineStyle = 'none';
            
            ax = gca;
            ax.YScale = 'linear';
            if strcmp(pluronic,'F87')
                ax.XLim = [36 44];
                ax.YLim = [6 8];
            else
                ax.XLim = [24 31];
                ax.YLim = [8 10];
            end
            
%             ax.XLim = [20 60];
            ax.XLabel.String = ['Concentration [ wt% ]'];
            ax.YLabel.String = 'Micelle radius [ nm ]';
            lgd = legend(markers,{'25 \circC','35 \circC','45 \circC','55 \circC'});
            lgd.Title.String = [pluronic];
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'north';
            lgd.LineWidth = 0.5; 
    end
end

if save_im == 'y'
    fig.PaperPosition = [0 0 12 10];  
    fig.PaperUnits = 'inches';
    plot_name = ['SAXS'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
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
            if k == 3
                continue 
            end
            subplot(1,2,k);
                p1 = errorbar(data.c,data.D,data.Dneg,data.Dpos);
                hold on
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.MarkerSize = MarkerSize;
                p1.LineStyle = 'none';
                p1.LineWidth = 2;

                ax = gca;
                ax.YScale = 'log';
                if  strcmp(pluronic, 'F127')
                    ax.XLim = [24 31];
                elseif strcmp(pluronic, 'F87')
                    ax.XLim = [36 44];
                elseif strcmp(pluronic, 'P123')
                     ax.XLim = [19 36];
                else
                end
                
                ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = [pluronic ,' wt%'];
                ax.YLabel.String = 'D/D_0';
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

            if save_im == 'y' && temperature == 55
                fig.PaperPosition = [0 0 12 5];  
                fig.PaperUnits = 'inches';
                plot_name = ['d v c ', protein_name];
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
            if k == 3
                continue 
            end
            subplot(1,2,k);
                p1 = errorbar(data.c,data.fm0,data.fmlb,data.fmub);
                hold on
                p1.Marker = 'o';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.MarkerSize = MarkerSize;
                p1.LineStyle = 'none';
                p1.LineWidth = 2;

                ax = gca;
                ax.YScale = 'linear';
                if  strcmp(pluronic, 'F127')
                    ax.XLim = [24 31];
                elseif strcmp(pluronic, 'F87')
                    ax.XLim = [36 44];
                elseif strcmp(pluronic, 'P123')
                    ax.XLim = [19 36];
                else
                end
                
                ax.YLim = [0 1]; 
                ax.XLabel.String = [pluronic ,' wt%'];
                ax.YLabel.String = 'Mobile Protein Fraction';

            if temperature == 55
                lgd = legend({'25 \circC','35 \circC','45 \circC','55 \circC'});
                lgd.Orientation = 'vertical';
                lgd.NumColumns = 2;
                lgd.Location = 'best';
                lgd.Title.String = [protein_name,' at T = '];
                lgd.LineWidth = 0.5; 
            else
            end        

                fig.PaperPosition = [0 0 12 5];  

            if save_im == 'y' && temperature == 55
                fig.PaperUnits = 'inches';
                plot_name = ['fm v c ', protein_name];
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
        if i >2
            continue
        end
        pluronic = pluronics{1};
        
        subplot(1,2,i);
        length(fieldnames(mainT.(protein_name)));
        c_count = 0;
            for concentrations = fieldnames(mainT.(protein_name).(pluronic))'

                c_count = c_count+1;
                
                if strcmp(protein_name, 'BSA') && c_count < 4 
                    continue
                end
                
                c = concentrations{1};
                data = mainT.(protein_name).(pluronic).(c);

                p1 = errorbar(data.T,data.D,data.Dneg,data.Dpos);
                hold on
                p1.Marker = 'd';
                p1.Color = colors3(c_count,:);
                p1.MarkerFaceColor = colors3(c_count,:);
                p1.LineStyle = 'none';
                p1.MarkerSize = MarkerSize;
                p1.LineWidth = 2;
            end
            
                ax = gca;
                ax.XLim = [20 60]; 
                ax.YScale = 'log';
                ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = 'Temperature [\circC]';
                ax.YLabel.String = 'D/D_0';
                
%                 if strcmp(protein_name,'BSA')
%                     lgd = legend(lgd_names.(protein_name).(pluronic));
%                     lgd.Orientation = 'vertical';
%                     lgd.NumColumns = 2;
%                     lgd.Location = 'best';
%                     lgd.Title.String = [protein_name,' in ',pluronic];
%                     lgd.LineWidth = 0.5; 
%                 else
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
%                 end
    end

        fig.PaperPosition = [0 0 12 5];  
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = ['d v t ', protein_name];
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
        if i >2
            continue
        end
        subplot(1,2,i);
        length(fieldnames(mainT.(protein_name)));
        
            c_count = 0;
            for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
                c_count = c_count+1;
                
                if strcmp(protein_name, 'BSA') && c_count < 4 
                    continue
                end
                c = concentrations{1};
                
                data = mainT.(protein_name).(pluronic).(c);

                p1 = errorbar(data.T,data.fm0,data.fmlb,data.fmub);
                hold on
                p1.Marker = 'o';
                p1.Color = colors3(c_count,:);
                p1.MarkerFaceColor = colors3(c_count,:);
                p1.LineStyle = 'none';
                p1.MarkerSize = MarkerSize;
                p1.LineWidth = 2;
            end

                
                ax = gca;
                ax.XLim = [20 60]; 
                ax.YScale = 'linear';
                ax.YLim = [0 1]; 
                ax.XLabel.String = 'Temperature [\circC]';
                ax.YLabel.String = 'Mobile Protein Fraction';
                
%                 if strcmp(protein_name,'BSA')
%                     lgd = legend(lgd_names.(protein_name).(pluronic));
%                     lgd.Orientation = 'vertical';
%                     lgd.NumColumns = 1;
%                     lgd.Location = 'northwest';
%                     lgd.Title.String = [protein_name,' in ',pluronic];
%                     lgd.LineWidth = 0.5; 
%                     lgd.Location = 'southwest';
%                 else
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name,' in ',pluronic];
                    lgd.LineWidth = 0.5; 
                    lgd.Location = 'southwest';
%                 end
    end
%     if strcmp(protein_name,'BSA')
%         fig.PaperPosition = [0 0 16 5]; 
%         
%     else
        fig.PaperPosition = [0 0 12 5];  
%     end
    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = ['f v t ', protein_name];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end


%% Plotting D vs rh from DLS (figure == all temperatutrs),(subplots == pluronic at all temperatures)
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
        if i >2
            continue
        end
        subplot(1,2,i);
        c_count = 0;
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            c = concentrations{1};
            data = mainT.(protein_name).(pluronic).(c);

            p1 = errorbar(data.rh_prot,data.D,data.Dneg,data.Dpos);
            hold on
            p1.Marker = 'd';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [1 6]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1e0]; 
        ax.XLabel.String = 'R_h Protein';
        ax.YLabel.String = 'D/D_0';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'northeast';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end

end
    fig.PaperPosition = [0 0 12 5];  

    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = ['d v rh ALL Temperatures'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
    
%% Plotting D vs rh from DLS (figure == 8 subplots),(subplots == pluronic at one Temperature)
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
        if i >2
            continue
        end
        
        c_count = 0;
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            c = concentrations{1};
            data = mainT.(protein_name).(pluronic).(c);
            
            for k = 1:4
                n = (4*(i-1)) + k;
                subplot(2,4,n);
                p1 = errorbar(data.rh_prot(k),data.D(k),data.Dneg(k),data.Dpos(k));
                hold on
                p1.Marker = 'd';
                p1.Color = colors3(c_count,:);
                p1.MarkerFaceColor = colors3(c_count,:);
                p1.LineStyle = 'none';
                p1.MarkerSize = MarkerSize-2;
                p1.LineWidth = 2;
                ax = gca;
                ax.XLim = [1 6]; 
                ax.YScale = 'log';
                ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = 'R_h Protein';
                ax.YLabel.String = 'D/D_0';
                
        
                if strcmp(protein_name,'BSA') && i ==2
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 2;
                    lgd.Location = 'southwest';
                    lgd.Title.String = [pluronic, ': T = ', num2str(data.T(k)), '\circC'];
                    lgd.LineWidth = 0.5; 
                    lgd.AutoUpdate = 'off';
                elseif strcmp(protein_name,'BSA') && i ==1
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 2;
                    lgd.Location = 'northeast';
                    lgd.Title.String = [pluronic, ': T = ', num2str(data.T(k)), '\circC'];
                    lgd.LineWidth = 0.5; 
                    lgd.AutoUpdate = 'off';
                end
            end
        end




    end

end
    fig.PaperPosition = [0 0 12*1.5 5*1.6];  

    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = ['d v rh 2x4'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
    
%% Plotting fm vs rh from DLS (figure == all temperatutrs),(subplots == pluronic at all temperatures)
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
        if i >2
            continue
        end
        subplot(1,2,i);
        c_count = 0;
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            c = concentrations{1};
            data = mainT.(protein_name).(pluronic).(c);

            p1 = errorbar(data.rh_prot,data.fm0,data.fmlb,data.fmub);
            hold on
            p1.Marker = 'o';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [1 6]; 
        ax.YScale = 'linear';
        ax.YLim = [0 1]; 
        ax.XLabel.String = 'R_{h,Protein}';
        ax.YLabel.String = 'Immobile Protein Fraction';
        yticks([0 0.2 0.4 0.6 0.8 1])

        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'southwest';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end


    end

end
    fig.PaperPosition = [0 0 12 5];  

    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = ['fm v rh ALL Temperatures'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
    
%% Plotting fm vs rh from DLS (figure == 8 subplots),(subplots == pluronic at one Temperature)
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
        if i >2
            continue
        end
        
        c_count = 0;
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            c = concentrations{1};
            data = mainT.(protein_name).(pluronic).(c);
            
            for k = 1:4
                n = (4*(i-1)) + k;
                subplot(2,4,n);
                p1 = errorbar(data.rh_prot(k),data.fm0(k),data.fmlb(k),data.fmub(k));
                hold on
                p1.Marker = 'o';
                p1.Color = colors3(c_count,:);
                p1.MarkerFaceColor = colors3(c_count,:);
                p1.LineStyle = 'none';
                p1.MarkerSize = MarkerSize-2;
                p1.LineWidth = 2;
                
                ax = gca;
                ax.XLim = [1 6]; 
                ax.YScale = 'linear';
                ax.YLim = [0 1]; 
                ax.XLabel.String = 'R_{h,Protein}';
                ax.YLabel.String = 'Mobile Protein Fraction';
                yticks([0 0.2 0.4 0.6 0.8 1])
                
                if strcmp(protein_name,'BSA')
                    lgd = legend(lgd_names.(protein_name).(pluronic));
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 1;
                    lgd.Location = 'southwest';
                    lgd.Title.String = [pluronic, ': T = ', num2str(data.T(k)), '\circC'];
                    lgd.LineWidth = 0.5; 
                    lgd.AutoUpdate = 'off';
                end
            end
        end
    end
end
    fig.PaperPosition = [0 0 12*1.5 5*1.6];  

    if save_im == 'y'
        fig.PaperUnits = 'inches';
        plot_name = ['fm v rh 2x4'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
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

        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c = concentrations{1};
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.rm;
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            p1.Marker = 'd';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [0.1 1.2]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1e0]; 
        ax.XLabel.String = 'R_{h,Protein} / R_{micelle}';
        ax.YLabel.String = 'D/D_0';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'northeast';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end
    end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 5];
    fig.PaperUnits = 'inches';
    plot_name = ['d v rh_rm ALL Temperatures'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end
%% Plotting fm vs rh/rm from DLS (figure == protein),(subplots == pluronic)
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

        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c = concentrations{1};
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.rm;
            p1 = errorbar(xdata,data.fm0,data.fmlb,data.fmub);
            hold on
            p1.Marker = 'o';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [0.1 1.2]; 
        ax.YScale = 'linear';
        ax.YLim = [0 1]; 
        ax.XLabel.String = 'R_{h,Protein} / R_{micelle}';
        ax.YLabel.String = 'Mobile Protein Fraction';

        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'northeast';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end
    end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 5];
    fig.PaperUnits = 'inches';
    plot_name = ['fm v rh_rm ALL Temperatures'];
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
        if i ==3
            break
        else
        end
        pluronic = pluronics{1};
        
        subplot(1,2,i);
            c_count = 0;

        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c = concentrations{1};
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.ro;
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            p1.Marker = 'd';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [0.1 1.2]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1e0]; 
        ax.XLabel.String = 'R_{h,Protein} / R_{octahedral}';
        ax.YLabel.String = 'D/D_0';
        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'northeast';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end
    end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 5];
    fig.PaperUnits = 'inches';
    plot_name = ['d v rh_ro ALL Temperatures'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end
%% Plotting fm vs rh/ro from DLS (figure == protein),(subplots == pluronic)
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

        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c = concentrations{1};
            c_count = c_count+1;
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.ro;
            p1 = errorbar(xdata,data.fm0,data.fmlb,data.fmub);
            hold on
            p1.Marker = 'o';
            p1.Color = colors3(c_count,:);
            p1.MarkerFaceColor = colors3(c_count,:);
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [0.1 1.2]; 
        ax.YScale = 'linear';
        ax.YLim = [0 1]; 
        ax.XLabel.String = 'R_{h,Protein} / R_{octahedral}';
        ax.YLabel.String = 'Mobile Protein Fraction';

        
        if strcmp(protein_name,'BSA')
            lgd = legend(lgd_names.(protein_name).(pluronic));
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 1;
            lgd.Location = 'southwest';
            lgd.Title.String = [pluronic];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end
    end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 5];
    fig.PaperUnits = 'inches';
    plot_name = ['fm v rh_ro ALL Temperatures'];
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
            c_count = 0;
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c = concentrations{1};
            c_count = c_count+1;
            
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);
            xdata = data.rh_prot./SAXSdata.rm;
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            if strcmp(pluronic,'F87')
                p1.Marker = 's';
            else
                p1.Marker = '^';
            end 
            p1.Color = colors3(c_count,:);
            if SAXSdata.extrapolated(1) == "y"
                p1.MarkerFaceColor = 'none';
            else
                p1.MarkerFaceColor = colors3(c_count,:);
            end
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [0.1 1.2]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1]; 
        ax.XLabel.String = 'R_{h, Protein} / R_{micelle}';
        ax.YLabel.String = 'D/D_0';
        
        if strcmp(protein_name,'HSA')
            lgd = legend(lgd_names.all);
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = ['wt% , pluronic'];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end
    end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 5];
    fig.PaperUnits = 'inches';
    plot_name = ['d vs rh_rm COMBO'];
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
        c_count = 0;
        for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
            c = concentrations{1};
            c_count = c_count+1;
            
            if strcmp(protein_name, 'BSA') && c_count < 4 
                continue
            end
            data = mainT.(protein_name).(pluronic).(c);
            SAXSdata = mainSAXS.(protein_name).(pluronic).(c);

            xdata = data.rh_prot./SAXSdata.ro;
            
            p1 = errorbar(xdata,data.D,data.Dneg,data.Dpos);
            hold on
            if strcmp(pluronic,'F87')
                p1.Marker = 's';
                p1.MarkerFaceColor = colors3(c_count,:);
            else
                p1.Marker = '^';
                p1.MarkerFaceColor = 'none';
            end   
            p1.Color = colors3(c_count,:);
%             if SAXSdata.extrapolated(1) == "y"
%                 p1.MarkerFaceColor = 'none';
%             else
%                 p1.MarkerFaceColor = colors3(c_count,:);
%             end

            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize;
            p1.LineWidth = 2;
        end
            
        ax = gca;
        ax.XLim = [0.1 1.2]; 
        ax.YScale = 'log';
        ax.YLim = [1e-4 1]; 
        ax.XLabel.String = 'R_{h, Protein} / R_{octahedral}';
        ax.YLabel.String = 'D/D_0';
        
        if strcmp(protein_name,'HSA')
            lgd = legend(lgd_names.all);
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'northeast';
            lgd.Title.String = ['wt% , pluronic'];
            lgd.LineWidth = 0.5; 
            lgd.AutoUpdate = 'off';
        else
        end
    end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 7 6];
    fig.PaperUnits = 'inches';
    plot_name = ['d vs rh_ro COMBO'];
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
        if i ==3
            break
        else
        end
        subplot(1,2,i);
        length(fieldnames(mainT.(protein_name)));
            c_count = 0;
            for concentrations = fieldnames(mainT.(protein_name).(pluronic))'
                c = concentrations{1};
                c_count = c_count+1;

                if strcmp(protein_name, 'BSA') && c_count < 4 
                    continue
                end
                data = mainT.(protein_name).(pluronic).(c);
                ydata = log(data.D_raw);
                xdata = 1./(data.T+273);

                p1 = plot(xdata,ydata);
                hold on
                p1.Marker = 'o';
                p1.Color = colors3(c_count,:);
                p1.MarkerFaceColor = colors3(c_count,:);
                p1.LineStyle = 'none';
            end
            
                ax = gca;
%                 ax.XLim = [20 60]; 
%                 ax.YScale = 'log';
%                 ax.YLim = [1e-4 1e0]; 
                ax.XLabel.String = '1/Temperature [1/K]';
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

    if save_im == 'y'
        fig.PaperPosition = [0 0 12 5];  
        fig.PaperUnits = 'inches';
        plot_name = ['ln[d] v T ', protein_name];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end
end

%% Plotting all D_raw v C (figure == protein),(subplots == pluronic)
close all
protein_count = 0;
fits = struct();
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
            if k == 3
                break
            else
            end
            data = main.(protein_name).(temperature_name).(pluronic);
            
            
            if strcmp(protein_name,'BSA')
                x = (data.c./100)';
                y = log(data.D_raw');
                x = x(4:6);
                y = y(4:6);

                yneg = 55.1*data.Dneg(4:6)';
                ypos = 55.1*data.Dpos(4:6)';
            else
                x = (data.c./100)';
                y = log(data.D_raw');
                yneg = 55.1*data.Dneg';
                ypos = 55.1*data.Dpos';
            end

            
            % fit the data 
            ft = fittype('interstitial_hopping_fit_ln(x,B,G)');
                 options = fitoptions(ft);
                  options.StartPoint = [0.1,0.3];
                  options.Lower = [0.01, 0];
                  options.Upper = [1000, 2];
                  options.TolFun = 1e-10;
                 options.Robust = 'LAR';  
            % perform the fit to the data
            [f,gof,output] = fit(x,y,ft,options);
            x_fit = linspace(0.1,0.5,10)'; 
            D_fit = exp(f(x_fit)); 

            fits.(protein_name).(temperature_name).(pluronic) = f;
            
            subplot(1,2,k);
                
                p1 = errorbar(data.c/100,data.D_raw,data.Dneg,data.Dpos);
                hold on
                p2 = plot(x_fit,D_fit);
                p2.Color = color;
                p1.Marker = 'd';
                p1.Color = color;
                p1.MarkerFaceColor = color;
                p1.LineStyle = 'none';
                p1.MarkerSize = MarkerSize-2;
                p1.LineWidth = 2;
                
                legend_markers(t) = p1;
                ax = gca;
                ax.YScale = 'log';
                
                if  temperature == 55
                    lgd = legend(legend_markers, {'25 \circC','35 \circC','45 \circC','55 \circC'});
                    lgd.Orientation = 'vertical';
                    lgd.NumColumns = 2;
                    lgd.Location = 'best';
                    lgd.Title.String = [protein_name];
                    lgd.LineWidth = 0.5; 
                else
                end  

                if  strcmp(pluronic, 'F127')
%                     ax.XLim = [17 31];
                    ax.XLim = [0.24 0.31];
                elseif strcmp(pluronic, 'F87')
%                     ax.XLim = [24 44];
                    ax.XLim = [0.36 0.44];
                elseif strcmp(pluronic, 'P123')
%                      ax.XLim = [19 36];
                     ax.XLim = [0.19 0.36];
                else
                end
                
                ax.YLim = [1e-2 5e1]; 
                ax.XLabel.String = [pluronic ,' wt fraction'];
                ax.YLabel.String = 'D [\mum^2s^{-1}]';
                yticks([1e-3 1e-2 1e-1 1e0 1e1 1e2])
        end
   
    end
    
        if save_im == 'y' && temperature == 55
            fig.PaperPosition = [0 0 12 5];  
            fig.PaperUnits = 'inches';
            plot_name = ['dRAW v c fit ' , protein_name];
            plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
            print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
        else
        end
end


%% Plotting all D_raw v C for F87 with SAXS params and DLS params
crystal_model_table_F87 = table('size',[16,4],'VariableTypes',{'double','double','double','string'},'variableNames',{'B','G','T','prot'});

close all
protein_count = 0;
n = 1;
global d a 
fig = figure(69);
set(fig,'Position',[2600 300 1500 500]);
fits = struct();
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
    

    t = 0;
    for temperatures = fieldnames(main.(protein_name))' %all25C, all35C etc
        t = t+1;
        temperature_name = temperatures{1};
        temperature = str2num(temperature_name(4:5));
        color = colors4(t,:);
        
        pluronic = 'F87';
        global pluronic_name T
        pluronic_name = pluronic; 
        T = temperature;
        data = main.(protein_name).(temperature_name).(pluronic);
        
        w = 1;
        a = zeros(3,1);
        for concs = fieldnames(mainSAXS.(protein_name).F87)'
            conc = concs{1};

            aw = mainSAXS.(protein_name).F87.(conc).a(t)/1000; % nm to micron convesion
            a(w) = aw;
            w = w +1;
        end
%         a = a(end); % matlab isnt playin with coefficients as vectors
        
        if strcmp(protein_name,'BSA')
            x = (data.c./100)';
            y = data.D_raw';
            x = x(4:6);
            y = y(4:6);
            
            yneg = 55.1*data.Dneg(4:6)';
            ypos = 55.1*data.Dpos(4:6)';
            
        else
            x = (data.c./100)';
            y = data.D_raw';
            yneg = 55.1*data.Dneg';
            ypos = 55.1*data.Dpos';
        end
        d = 2*main.(protein_name).(temperature_name).(pluronic).rh_prot(1);
        
%         for concs = fieldnames(mainSAXS.BSA.F87)
            

        % fit the data 
        
        % fit the data 
        x0 = [10 1];
        xlb = [0 0];
        xub = [1000000 1e2];
        xdata = x;
        ydata = log(y);
        
        output_coefficients = lsqcurvefit(@interstitial_hopping_SAXS_DLS_lsq,x0,xdata,ydata,xlb,xub);
        B = output_coefficients(1);
        G = output_coefficients(2);
        
%         ft = fittype('interstitial_hopping_SAXS_DLS(x,B,G)');
%              options = fitoptions(ft);
%               options.StartPoint = [0.1,0.3];
%               options.Lower = [0.01, 0];
%               options.Upper = [1000, 2];
%               options.TolFun = 1e-10;
%              options.Robust = 'LAR';  
%         % perform the fit to the data
%         y_fit = log(y);
%         [f,gof,output] = fit(x,y_fit,ft,options);
%         x_fit = linspace(0.35,0.5,10)'; 
%         D_fit = exp(f(x_fit)); 
%         
        x_fit = linspace(0.3,0.5,10)'; 
        D_fit = exp(interstitial_hopping_SAXS_DLS_lsq([B;G],x_fit)); 

        subplot(2,2,protein_count);

            p1 = errorbar(x,y,yneg,ypos);
            hold on
            p2 = plot(x_fit,D_fit);
            p2.Color = color;
            p2.LineWidth = 2;
            
            
            p1.Marker = 'd';
            p1.Color = color;
            p1.MarkerFaceColor = color;
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
            legend_markers(t) = p1;

            ax = gca;
            ax.YScale = 'log';
            ax.XLim = [0.36 0.44];
            ax.YLim = [1e-2 5e1]; 
            ax.XLabel.String = [pluronic ,' wt fraction'];
            ax.YLabel.String = 'D [\mum^2s^{-1}]';
            yticks([1e-3 1e-2 1e-1 1e0 1e1 1e2])
            
        if temperature == 55
            lgd = legend(legend_markers,{'25 \circC','35 \circC','45 \circC','55 \circC'});
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = [protein_name,' in ',pluronic];
            lgd.LineWidth = 0.5; 
        else
        end        
        temp =struct('fits',f);
        fits.('fits')(n) = temp;
        crystal_model_table_F87.B(n) = exp(B);
        crystal_model_table_F87.G(n) = G;
        crystal_model_table_F87.T(n) = temperature;
        crystal_model_table_F87.prot(n) = protein_name;
        n = n+1;
      end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 10];  
    fig.PaperUnits = 'inches';
    plot_name = ['dRAW v c SAXS_fit F87'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end

%% Plotting all D_raw v C for F127 with SAXS params and DLS params
crystal_model_table_F127 = table('size',[16,4],'VariableTypes',{'double','double','double','string'},'variableNames',{'B','G','T','prot'});

protein_count = 0;
n = 1;
global d a 
fig = figure(75);
set(fig,'Position',[1700 300 1200 500]);
fits = struct();
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
    

    t = 0;
    for temperatures = fieldnames(main.(protein_name))' %all25C, all35C etc
        t = t+1;
        temperature_name = temperatures{1};
        temperature = str2num(temperature_name(4:5));
        color = colors4(t,:);
        
        pluronic = 'F127';
        global pluronic_name T
        pluronic_name = pluronic; 
        T = temperature;
        data = main.(protein_name).(temperature_name).(pluronic);
        
        w = 1;
        a = zeros(3,1);
        for concs = fieldnames(mainSAXS.LYS.F127)'
            conc = concs{1};

            aw = mainSAXS.(protein_name).F127.(conc).a(t)/1000; % nm to micron convesion
            a(w) = aw;
            w = w +1;
        end
%         a = a(end); % matlab isnt playin with coefficients as vectors
        
        if strcmp(protein_name,'BSA')
            x = (data.c./100)';
            y = data.D_raw';
            x = x(4:6);
            y = y(4:6);
            
            yneg = 55.1*data.Dneg(4:6)';
            ypos = 55.1*data.Dpos(4:6)';
            
        else
            x = (data.c./100)';
            y = data.D_raw';
            yneg = 55.1*data.Dneg';
            ypos = 55.1*data.Dpos';
        end
        d = 2*main.(protein_name).(temperature_name).(pluronic).rh_prot(1);
        
        % fit the data 
        x0 = [10 1];
        xlb = [0 0];
        xub = [1000 1e2];
        xdata = x;
        ydata = log(y);
        
        output_coefficients = lsqcurvefit(@interstitial_hopping_SAXS_DLS_lsq,x0,xdata,ydata,xlb,xub);
        B = output_coefficients(1);
        G = output_coefficients(2);
        
%         ft = fittype('interstitial_hopping_SAXS_DLS(x,B,G)');
%              options = fitoptions(ft);
%               options.StartPoint = [0.1,0.3];
%               options.Lower = [0.01, 0];
%               options.Upper = [1000, 2];
%               options.TolFun = 1e-10;
%              options.Robust = 'LAR';  
        % perform the fit to the data

        x_fit = linspace(0.2,0.4,10)'; 
        D_fit = exp(interstitial_hopping_SAXS_DLS_lsq([B;G],x_fit)); 

        subplot(2,2,protein_count);

            p1 = errorbar(x,y,yneg,ypos);
            hold on
            p2 = plot(x_fit,D_fit);
            p2.Color = color;
            p2.LineWidth = 2;
            
            p1.Marker = 'd';
            p1.Color = color;
            p1.MarkerFaceColor = color;
            p1.LineStyle = 'none';
            p1.MarkerSize = MarkerSize-2;
            p1.LineWidth = 2;
            legend_markers(t) = p1;

            ax = gca;
            ax.YScale = 'log';
            ax.XLim = [0.24 0.31];
            ax.YLim = [1e-2 5e1]; 
            ax.XLabel.String = [pluronic ,' wt fraction'];
            ax.YLabel.String = 'D [\mum^2s^{-1}]';
            yticks([1e-3 1e-2 1e-1 1e0 1e1 1e2])
            
        if temperature == 55
            lgd = legend(legend_markers,{'25 \circC','35 \circC','45 \circC','55 \circC'});
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;
            lgd.Location = 'best';
            lgd.Title.String = [protein_name,' in ',pluronic];
            lgd.LineWidth = 0.5; 
        else
        end        
        temp =struct('fits',f);
        fits.('fits')(n) = temp;
        crystal_model_table_F127.B(n) = exp(B);
        crystal_model_table_F127.G(n) = G;
        crystal_model_table_F127.T(n) = temperature;
        crystal_model_table_F127.prot(n) = protein_name;
        n = n+1;
      end
end
if save_im == 'y'
    fig.PaperPosition = [0 0 12 10];  
    fig.PaperUnits = 'inches';
    plot_name = ['dRAW v c SAXS_fit F127'];
    plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
    print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else
end


%% now plotting the B and G for F87
fig = figure(70);
set(fig,'Position',[1700 300 1200 500]);

subplot(1,2,1)
    hold on
    p1 = plot(crystal_model_table_F87.T(1:4),crystal_model_table_F87.B(1:4));
    p1.Marker = 'd';
    p1.LineStyle = 'none';
    p1.Color = colors4(4,:);
    p1.MarkerFaceColor = colors4(4,:);
    p1.MarkerSize = MarkerSize-2;
    p1.LineWidth = 2;
    
    p2 = plot(crystal_model_table_F87.T(5:8),crystal_model_table_F87.B(5:8));
    p2.Marker = 's';
    p2.LineStyle = 'none';   
    p2.Color = colors4(3,:);
    p2.MarkerFaceColor = colors4(3,:);
    p2.MarkerSize = MarkerSize-2;
    p2.LineWidth = 2;
    
    p3 = plot(crystal_model_table_F87.T(9:12),crystal_model_table_F87.B(9:12));
    p3.Marker = 'o';
    p3.LineStyle = 'none';  
    p3.Color = colors4(2,:);
    p3.MarkerFaceColor = colors4(2,:);
    p3.MarkerSize = MarkerSize-2;
    p3.LineWidth = 2;
    
    p4 = plot(crystal_model_table_F87.T(13:16),crystal_model_table_F87.B(13:16));
    p4.Marker = '*';
    p4.LineStyle = 'none';  
    p4.Color = colors4(1,:);
    p4.MarkerFaceColor = colors4(1,:);
    p4.MarkerSize = MarkerSize-2;
    p4.LineWidth = 2;
    
    ax = gca;
    ax.YScale = 'log';
    ax.XLim = [20 60];
    ax.XLabel.String = 'Temperature';
    ax.YLabel.String = 'Beta';
    ax.YLim = [1e2 1e7]; 

    lgd = legend([p1,p2,p3,p4],'BSA','LYS','CHA','HSA');
    lgd.Location = 'southwest';
    lgd.Title.String = ['Protein'];
    lgd.LineWidth = 0.5; 
    

subplot(1,2,2)
    hold on
    p1 = plot(crystal_model_table_F87.T(1:4),crystal_model_table_F87.G(1:4));
    p1.Marker = 'd';
    p1.LineStyle = 'none';
    p1.Color = colors4(4,:);
    p1.MarkerFaceColor = colors4(4,:);
    p1.MarkerSize = MarkerSize-2;
    p1.LineWidth = 2;
    
    p2 = plot(crystal_model_table_F87.T(5:8),crystal_model_table_F87.G(5:8));
    p2.Marker = 's';
    p2.LineStyle = 'none'; 
    p2.Color = colors4(3,:);
    p2.MarkerFaceColor = colors4(3,:);
    p2.MarkerSize = MarkerSize-2;
    p2.LineWidth = 2;
    
    p3 = plot(crystal_model_table_F87.T(9:12),crystal_model_table_F87.G(9:12));
    p3.Marker = 'o';
    p3.LineStyle = 'none';  
    p3.Color = colors4(2,:);
    p3.MarkerFaceColor = colors4(2,:);
    p3.MarkerSize = MarkerSize-2;
    p3.LineWidth = 2;
    
    p4 = plot(crystal_model_table_F87.T(13:16),crystal_model_table_F87.G(13:16));
    p4.Marker = '*';
    p4.LineStyle = 'none';  
    p4.Color = colors4(1,:);
    p4.MarkerFaceColor = colors4(1,:);
    p4.MarkerSize = MarkerSize-2;
    p4.LineWidth = 2;
    
    ax = gca;
    ax.XLim = [20 60];
    ax.YScale = 'Linear';
    ax.XLabel.String = 'Temperature';
    ax.YLabel.String = 'Gamma';
     ax.YLim = [0 30]; 

    lgd = legend([p1,p2,p3,p4],'BSA','LYS','CHA','HSA');
    lgd.Location = 'northwest';
    lgd.Title.String = ['Protein'];
    lgd.LineWidth = 0.5; 
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;

    if save_im == 'y' 
fig.PaperPosition = [0 0 12 5];  
        fig.PaperUnits = 'inches';
        plot_name = ['informed_hope_model F87'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end


% now plotting the B and G for F127
fig = figure(71);
set(fig,'Position',[1700 300 1200 500]);

subplot(1,2,1)
    hold on
    p1 = plot(crystal_model_table_F127.T(1:4),crystal_model_table_F127.B(1:4));
    p1.Marker = 'd';
    p1.LineStyle = 'none';
    p1.Color = colors4(4,:);
    p1.MarkerFaceColor = colors4(4,:);
    p1.MarkerSize = MarkerSize-2;
    p1.LineWidth = 2;
    
    p2 = plot(crystal_model_table_F127.T(5:8),crystal_model_table_F127.B(5:8));
    p2.Marker = 's';
    p2.LineStyle = 'none';   
    p2.Color = colors4(3,:);
    p2.MarkerFaceColor = colors4(3,:);
    p2.MarkerSize = MarkerSize-2;
    p2.LineWidth = 2;
    
    p3 = plot(crystal_model_table_F127.T(9:12),crystal_model_table_F127.B(9:12));
    p3.Marker = 'o';
    p3.LineStyle = 'none';  
    p3.Color = colors4(2,:);
    p3.MarkerFaceColor = colors4(2,:);
    p3.MarkerSize = MarkerSize-2;
    p3.LineWidth = 2;
    
    p4 = plot(crystal_model_table_F127.T(13:16),crystal_model_table_F127.B(13:16));
    p4.Marker = '*';
    p4.LineStyle = 'none';  
    p4.Color = colors4(1,:);
    p4.MarkerFaceColor = colors4(1,:);
    p4.MarkerSize = MarkerSize-2;
    p4.LineWidth = 2;
    
    ax = gca;
    ax.YScale = 'log';
    ax.XLim = [20 60];
    ax.XLabel.String = 'Temperature';
    ax.YLabel.String = 'Beta';
    ax.YLim = [1e2 1e7]; 

    lgd = legend([p1,p2,p3,p4],'BSA','LYS','CHA','HSA');
    lgd.Location = 'northwest';
    lgd.Title.String = ['Protein'];
    lgd.LineWidth = 0.5; 

subplot(1,2,2)
    hold on
    p1 = plot(crystal_model_table_F127.T(1:4),crystal_model_table_F127.G(1:4));
    p1.Marker = 'd';
    p1.LineStyle = 'none';
    p1.Color = colors4(4,:);
    p1.MarkerFaceColor = colors4(4,:);
    p1.MarkerSize = MarkerSize-2;
    p1.LineWidth = 2;
    
    p2 = plot(crystal_model_table_F127.T(5:8),crystal_model_table_F127.G(5:8));
    p2.Marker = 's';
    p2.LineStyle = 'none'; 
    p2.Color = colors4(3,:);
    p2.MarkerFaceColor = colors4(3,:);
    p2.MarkerSize = MarkerSize-2;
    p2.LineWidth = 2;
    
    p3 = plot(crystal_model_table_F127.T(9:12),crystal_model_table_F127.G(9:12));
    p3.Marker = 'o';
    p3.LineStyle = 'none';  
    p3.Color = colors4(2,:);
    p3.MarkerFaceColor = colors4(2,:);
    p3.MarkerSize = MarkerSize-2;
    p3.LineWidth = 2;
    
    p4 = plot(crystal_model_table_F127.T(13:16),crystal_model_table_F127.G(13:16));
    p4.Marker = '*';
    p4.LineStyle = 'none';  
    p4.Color = colors4(1,:);
    p4.MarkerFaceColor = colors4(1,:);
    p4.MarkerSize = MarkerSize-2;
    p4.LineWidth = 2;
    
    ax = gca;
    ax.XLim = [20 60];
    ax.YScale = 'Linear';
    ax.XLabel.String = 'Temperature';
    ax.YLabel.String = 'Gamma';
    ax.YLim = [0 30]; 

    lgd = legend([p1,p2,p3,p4],'BSA','LYS','CHA','HSA');
    lgd.Location = 'northwest';
    lgd.Title.String = ['Protein'];
    lgd.LineWidth = 0.5; 
            lgd.Orientation = 'vertical';
            lgd.NumColumns = 2;

  
    if save_im == 'y' 
        fig.PaperPosition = [0 0 12 5];
        fig.PaperUnits = 'inches';
        plot_name = ['informed_hope_model F127'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
    else
    end

%% 
