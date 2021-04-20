% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - move the plot folder to the outputs folder where the matlab data is
% stored
% - move experiments off of this hard-drive 
% - revisit the movie making to be able to see drift, etc
% - look up what the laser focus knob is actually doing
% - work on the master step3 file so can compare data faster
% - error flag for initial bleaching depth 

%%%%%%%%%%% primary
% - is photobleaching a permanent phenomena? irreversible?
% - why is radius so different than cheng??
% - cap fm off at one for the fits?

%%%%%%%%%%% thoughts
% - uniformity metrics for the bleached spot and the reference region

%%%%%%%%%%% Experiments to run %%%%%%%%%%%%%%%%%%%%%
% - energy balance on laser to make sure it isnt heating up

%%%%%%%%%%% Diffusivity fitting %%%%%%%%%%%%%%%%%%%%%5
% - other models for FRAP and normalization?
% - apply Ionic hopping mechanism model to frap data to see if that is a
% good metric?

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
    folder1 = 'F127_BSA_45C';
    folder2 = 'trial_6';

% Part 2: Would you like to save the plots generated? 
% Note:         Select 'y' or 'n'
save_im = 'y'; 

% Part 3: Are you troubleshooting the fits? 
% Note:         Select 'y' or 'n'
trouble = 'n'; 

%% %%%%%%%%%%%%%%%% Parameters Section: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Data Folder, structures, and sample ID data
global boxfolder outputfolder plotfolder datafolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    datafolder = fullfile(boxfolder,folder1,folder2); % full path where outputs structures will be saved
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,folder1,folder2,'z_plots'); % plot folder within the data folder for the defined experiment
    
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

[id,fits,fig,stats] = fun_radius_finder(id,fits,save_im);

disp("bleached circles found in " + string(round(toc)) +" s");
%% Reading in the rest of the frap data 
tic 

[alldata] = fun_frap_frames(id,fits,alldata);

disp(["Data loaded in " + string(round(toc)) + " s."])
%% Performing fits of the normalized data Cheng style
tic 

[fits] = fun_intensity_fits(fits,alldata);

disp(["all data fit in " + string(round(toc)) + " s."])
%%
% adding variables to id that can be useful for plotting 
% need for loop to unpack
i = 0;
pd = struct(); %plot data structure
for field = fieldnames(id)'
    position = field{1};
    i = i+1;
    pd.('r')(i) = fits.(position).FWHM .* 0.5 .*fits.(position).pixel_size;
    pd.('rlb')(i) = (fits.(position).FWHM - fits.(position).dFWHM(1)).* 0.5 .*fits.(position).pixel_size;
    pd.('rub')(i) = (fits.(position).dFWHM(2) - fits.(position).FWHM).* 0.5 .*fits.(position).pixel_size;
    pd.('c')(i) = id.(position).plwt;
    pd.('D')(i) = fits.(position).D/55.1; % from cheng to normalize by free soln BSA D
    pd.('Dneg')(i) = fits.(position).errD(1)/55.1;
    pd.('Dpos')(i) = fits.(position).errD(2)/55.1;
    pd.('fm0')(i) = fits.(position).fm0;   
    pd.('fm')(i) = fits.(position).fit_info.f;
    pd.('fmlb')(i) = fits.(position).fit_info.f - fits.(position).ci(1,1); 
    pd.('fmub')(i) =  fits.(position).ci(2,1) - fits.(position).fit_info.f;
end

% plotting the analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array

    c_temp = [num2str(id.(position).plwt), ' wt%'];
    fig = figure('name',c_temp,'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    
    % plotting data from image intensities, without the fits 
    t = [alldata.(position).dt]';
    t = t(2:end);
    I = [alldata.(position).I_ti]';
    I = I(2:end);
    IN = [alldata.(position).IN_ti]';
    IN = IN(2:end);
    ref_ratio = [alldata.(position).ref_ratio]';
    ref_ratio = ref_ratio(2:end);
    Ifitparam = fits.(position).('Ifitparam');

    % plotting the data
    hold on
%     plot(t,I,'x','color',C(40,:))
    plot(t,IN,'d','color',C(40,:),'markerfacecolor',C(40,:))
    % if not troubleshooting, plot the fits as well
    if trouble == 'n'
    t_fit = [fits.(position).t_fit]';
    I_fit = [fits.(position).I_fit]';
    plot(t_fit,I_fit,'--')
    xc = fits.(position).xc;
    shapec = fits.(position).shapec;
    fill(xc',shapec,'r','FaceAlpha',0.1,'LineStyle','--','linewidth',1);
    plot([0,t_fit(end)],[Ifitparam,Ifitparam])
    else 
    end
    title(position);
    axis([0,max(t),0,1.2]);
    xlabel('Time [s]');
    ylabel('Normalized Intensity');

% save images if selected above
if save_im == 'y'
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 8 6];             % define location to save the images 
        a = fieldnames(id);
        pp = a{1};
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
        plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_intensity_fit'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png  
else    
end

end

% diffusivity data vs conc
%pd = plot data structure
if trouble == 'n'
fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 

subplot(1,3,1)
    set(gca,'YScale','log');
    hold all
    errorbar(pd.c,pd.D,pd.Dneg,pd.Dpos,'db')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.7*min(pd.D) 2*max(pd.D)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('D/D_0 [\mum^2s^{-1}]')

subplot(1,3,2)
    hold on
    errorbar(pd.c,pd.r,pd.rlb,pd.rub,'or')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.8*min(pd.r) 1.2*max(pd.r)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('laser radius [um]')
    
subplot(1,3,3)
    hold on
%     plot(pd.c,pd.fm,'om')
    errorbar(pd.c,pd.fm,pd.fmlb,pd.fmub,'om')
    xlabel([id.(position).plur,' wt%'])
    ylabel('mobile fraction [fm]')
    axis([0.9*min(pd.c) 1.1*max(pd.c) 0.5 1.2])
else
end

% save images if selected above
if save_im == 'y'
        a = fieldnames(id);
        posn = a{1};
        fig.PaperUnits = 'inches';
        fig.PaperPosition = [0 0 10 4];             % define location to save the images 
        struct_name = [id.(posn).plur,'_',id.(posn).prot,'_',id.(posn).temp,'_',folder2];
        plot_name = [struct_name,'_DvC'];
        plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
        print(fig,plot_path, '-painters', '-dpng', '-r600')    % saving the figure as a high quality png    
else    
end
%% Exporting the data structures for future use.
% can also save individual fields if we want to get crazy
struct_name = [folder1,'_',folder2]; % specify file name prefix
cd(outputfolder);% first cd to output folder (where we want to save the data)

alldata_name = [struct_name,'_','data','.mat'];% save alldata struct
save(alldata_name,'alldata'); 
fits_name = [struct_name,'_','fits','.mat']; % save fits struct
save(fits_name,'fits');
id_name = [struct_name,'_','id','.mat']; % save id struct
save(id_name,'id');
id_name = [struct_name,'_','pd','.mat']; % save pd struct (plot data)
save(id_name,'pd');

cd(mainfolder);% cd back to main folder
%% complete
disp('Analysis Complete')

