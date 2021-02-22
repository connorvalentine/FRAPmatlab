% FRAP analysis V1 
% Connor Valentine
%% to do list
%%%%%%%%%%% urgent
% - recomment/clean code
% - save thge D and c vectors for plotting into a new structure to export
% - add mobile fraction vs C to output
% - make a sheet showing which experiments are good/favorites
% - add background flattening, and "inverted peak finding"

%%%%%%%%%%% primary
% - is photobleaching a permanent phenomena? irreversible?
% - concentric ring reference region?
% - why is P123 data so much higher?
% - why is radius so different than cheng??
% - multiple reference regions to help ignore irregularities
% - why are Ifit values so low?
% - error analysis of r_i vs tau
% - Ri should be the same more or less for every one- flag for this error
% - error propagation of std_dev of radius?

%%%%%%%%%%% thoughts
% - uniformity metrics for the bleached spot and the reference region
% - reference sample in each array? HDFL gel?
% - why is ifit taking so long

%%%%%%%%%%% Experiments to run %%%%%%%%%%%%%%%%%%%%%
% - Qunatify I0 vs bleaching time 
% - energy balance on laser to make sure it isnt heating up

%%%%%%%%%%% Diffusivity fitting %%%%%%%%%%%%%%%%%%%%%5
% - make fitting go faster
% - other models for FRAP and normalization?
% - confirm that the diffusivity is calculated correctly
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
    folder1 = 'F87_BSA_35C';
    folder2 = 'trial_1';

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

%% Reading in First Frap image to find bleached circle
tic; % timing
[id,fits,fig] = fun_first_bleached_frame(id,fits,save_im);

disp("bleached circles found in " + string(round(toc)) +" s");

%% Reading in the pre-bleach images next
[fits] = fun_prebleach_frame(id,fits);

%% Reading in the rest of the frap data 
   tic % timing
for field = fieldnames(id)'
    position = field{1};
 
    % Make a list of each image name in our position folder
    folder3 = 'frap';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    
    % define metadata filename for the current position
    global metadata_filename 
    txt = dir(fullfile(data(1).folder,'*.txt'));% lists all files with .txt ending
    metadata_filename = fullfile(txt.folder,txt.name);
    
    % pull out the prebleach image info
    refI_t0 = fits.(position).refI_t0;
    I_t0 = fits.(position).I_t0;
    idx_ref = fits.(position).idx_ref;
    idx_bleach = fits.(position).idx_bleach;

    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time_V4(data,t); 
        fits.(position).('pixel_size') = id.(position).pixel_size_manual;
        
%         if data(1).pixel_size == 0
%             fits.(position).('pixel_size') = 0.323; %for 20x in case it was not found 
%             fits.(position).('pixel_size') = 0.645; % 10x case 
%         else
%             fits.(position).('pixel_size') = data(1).pixel_size;
%         end

        % add pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        im = double(imread(f1));
        I_ti = mean(im(idx_bleach))./I_t0; % normalized by initial Intensity
        refI_ti = mean(im(idx_ref));
        if refI_ti ==0
            refI_ti = 0.001;
        else
        end
        normalized_i = (I_ti).*(refI_t0/refI_ti);
       % normalized_i = (I_ti./I_t0);

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('refI_ti') = refI_ti;
        data(t).('ref_ratio') = (refI_t0/refI_ti);
        data(t).('IN_ti') = normalized_i;
        alldata.(position) = data;
    end
    
    % add end frame to plots 
    radius = fits.(position).radius;
    center = fits.(position).center;
    imbds = fits.(position).imbds;
    disp([position,' loaded'])
end
a = round(toc);
disp(['Data loaded in ',num2str(a),' s.'])


%% Adding some additional things to the data structure
% 1. elapsed time from first "initial" frame
% 2. fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam

for field = fieldnames(alldata)'
    tic % timing start
    position = field{1}; % use{} bcuz field is a cell array
    % elapsed time from received time of first frap frame.
    time1 = alldata.(position)(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(alldata.(position))
        t = alldata.(position)(i).r_time;
        t = datevec(t);
        dt = etime(t,time1); % calculate elapsed time from the first image in seconds
        alldata.(position)(i).('dt') = dt; %save to data structure
    end
    
    if trouble == 'n' % ignores the fit process if trouble == 'y'
    % fitting the data to the equation in Cheng paper for FRAP
    % use x for the variable instead of t
    x = [alldata.(position).dt]';
    x = x(t1:end); % could have to adapt this later
    t_fit = linspace(0,x(end),250)';
    norm_i = [alldata.(position).IN_ti]';
    norm_i = norm_i(t1:end);
    
    % adding fit parameters 
    global Ifit
        I_t1 = alldata.(position)(2).I_ti; % skip laser frame
        R_t1 = alldata.(position)(2).refI_ti; %% check
        I_t0 = fits.(position).I_t0;
        R_t0 = fits.(position).refI_t0;
        Ifit = (I_t1./I_t0).*(R_t0./R_t1);
        disp(["Ifit= " + string(Ifit)])
    
    % now setting fit options and fitting
    % bounds on f are set by + or - 5% of the last intensity value now
    fm0 = norm_i(end);
    fm0_lb = 0.95*fm0
    fm0_ub = 1.05*fm0
    if fm0_ub > 1
        fm0_ub = 1
    else
    end
    ft = fittype('FRAPfit(x,f,k,tau)');
        options = fitoptions(ft); 
        options.StartPoint = [fm0 1 1000];
        options.Lower = [fm0_lb,0.1,1];
        options.Upper = [fm0_ub,10,20000];
        options.Robust = 'Bisquare';  
    [f,gof] = fit(x,norm_i, ft,options);
    I_fit = f(t_fit);
    
    % confidence intervales from the fit, f
    ci = confint(f,0.95);
    x_ci = t_fit;
    y_ci = predint(f,x_ci,0.95);
    %shapec and xc can be used to plot the shaded confidence region
    xc = [x_ci;flip(x_ci)];
    curve1 = y_ci(:,1);
    curve2 = y_ci(:,2);
    shapec = [curve1; flip(curve2)];
    
    % calculated diffusivity value
    ri = fits.(position).radius.*fits.(position).pixel_size;
    tau = f.tau;
    dtau = ci(:,3);
    D = (ri.^2)./(4*tau);
    Dlow = ((0.95*ri).^2)./(4*max(dtau)); % adds 5% variation to ri
    Dhigh = ((1.05*ri).^2)./(4*min(dtau)); % adds 5% variation to ri
    Dneg = D-Dlow;
    Dpos = Dhigh-D;
    errD = [Dneg Dpos];
    % creating new data structure for fit data
    fits.(position).('fit_info') = f;
    fits.(position).('ci') = ci;
    fits.(position).('D') = D;
    fits.(position).('errD') = errD;
    fits.(position).('gof') = gof;
    fits.(position).('I_fit') = I_fit;
    fits.(position).('t_fit') = t_fit;
    fits.(position).('xc') = xc;
    fits.(position).('shapec') = shapec;
    
    % check fm value
    I_inf = alldata.(position)(end).IN_ti;
    I_1 = alldata.(position)(1).IN_ti;
    fm_temp = (I_inf - I_1)./(I_t0 - I_1);
    fits.(position).('fm_calc') = fm_temp;
    % output values from fit 
    b = round(toc); % timing
    disp([position,' data fit in ',num2str(b),' s.'])
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(round(f.tau))]);
    else
    end
    % timing 
    b = round(toc);
end
disp('time added')

%% Exporting the data structures for future use.
% can also save individual fields if we want to get crazy
% first specify output folder (where we want to save the data
% maybe would be good to do conc instead of pos# here? 

% fix this section to just use folder1 and folder2 for naming
% then can also get temperature to be a number 
cd(outputfolder);
struct_name = [folder1,'_',folder2]; % specify file name prefix
alldata_name = [struct_name,'_','data','.mat'];
fits_name = [struct_name,'_','fits','.mat'];
id_name = [struct_name,'_','id','.mat']; % fwd id struct along with it.
save(alldata_name,'alldata'); 
save(fits_name,'fits');
save(id_name,'id');
cd(mainfolder);
%% plotting the analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array
    c_temp = id.(position).plwt;
    c_temp = [num2str(c_temp), ' wt%'];
    fig = figure('name',c_temp,'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    
    % plotting data from image intensities, without the fits 
    t = [alldata.(position).dt]';
    t = t(t1:end);
    I = [alldata.(position).I_ti]';
    I = I(t1:end);
    IN = [alldata.(position).IN_ti]';
    IN = IN(t1:end);
    ref_ratio = [alldata.(position).ref_ratio]';
    ref_ratio = ref_ratio(t1:end);

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



%% diffusivity data vs conc
if trouble == 'n'
fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 

% need for loop to unpack
i = 0;
fm = [];
fmlb = [];
fmub = [];
for field = fieldnames(alldata)'
    position = field{1};
    i = i+1;
    D = fits.(position).D/55.1;
    c = id.(position).plwt;
    pxl = fits.(position).pixel_size;
    r = fits.(position).radius*pxl;
    rall(i) =r;
    conc(i) = c;
    Dall(i) = D; % from cheng to normalize by free soln BSA D
    Dneg = fits.(position).errD(1)/55.1;
    Dpos = fits.(position).errD(2)/55.1;
    Dnegall(i) = Dneg;
    Dposall(i) = Dpos;
    fm(i) = fits.(position).fit_info.f;
    fmlb(i) = fm(i) - fits.(position).ci(1,1); 
    fmub(i) =  fits.(position).ci(2,1) -fm(i);
end

subplot(1,3,1)
    set(gca,'YScale','log');
    hold all
    errorbar(conc,Dall,Dnegall,Dposall,'db')
    axis([0.9*min(conc) 1.1*max(conc) 0.7*min(Dall) 2*max(Dall)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('D/D_0 [\mum^2s^{-1}]')

subplot(1,3,2)
    hold on
    plot(conc,rall,'or')
    axis([0.9*min(conc) 1.1*max(conc) 0.8*min(rall) 1.2*max(rall)])
    xlabel([id.(position).plur,' wt%'])
    ylabel('laser radius [um]')
    
subplot(1,3,3)
    hold on
    errorbar(conc,fm,fmlb,fmub,'om')
    axis([0.9*min(conc) 1.1*max(conc) 0.8 1.1])
    xlabel([id.(position).plur,' wt%'])
    ylabel('mobile fraction [fm]')
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

%% complete
disp('Analysis Complete')