% FRAP analysis V1 
% Connor Valentine
%% to do list
% - have a more rigorous way of choosing a smaller radius for pixel analysis
% - remove speckles?
% - save plots into a results folder with names
% - make text file, or matlab file, with the starting inputs to run the
% script in each subfolder. for instance the file name, concentration in
% order, etc.
%%%%%%%%%%%%%%%%%future future future goals %%%%%%%%%%%%%%%%%%%%%%%%%
% - split this into two matlab scripts, one to load the pics and one to
% analyze 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - first script should output for each position: 
% - a data structure saved with:
% 1. Identifying information: what day this data was taken
% 2. temperature, concentration of pluronic, conc of BSA
% 2.1 necessary intensity data, etc for each frame
% 3. save them all to a master folder? could make fitting the data easier
% later, perhaps with filenames: Plur_conc_temp_protein?

% - second script for fitting the data
% 1. takes in the strctures from 1) and fits the data

% - make fitting go faster
% - confirm that the diffusivity is calculated correctly
% - add concentration too
% - apply Ionic hopping mechanism model to frap data to see if that is a
% good metric?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5%
%% Part 1: Initialize some basic parameters and defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change default axes fonts.
    set(0,'DefaultAxesFontName', 'Arial')
    set(0,'DefaultAxesFontSize',9.5)
% Change default text fonts.
    set(0,'DefaultTextFontname', 'Arial','DefaultTextFontSize',9.5)
% default axis settings for better plots
    set(groot,'defaultaxeslinewidth',1)
    set(groot,'DefaultLineLineWidth',1)
    set(groot, 'DefaultAxesBox', 'on')
    set(groot, 'DefaultAxesYGrid', 'off')
    set(groot, 'DefaultAxesXGrid', 'off')
% Reset the environment
    clear all;
    clear global all;
    clc
    close all;
% Establish the main directory, 
    global mainfolder 
    mainfolder = cd;
    
%% choose the folders
% prebleach folder must be named 'prebleach'
% frap folder must be named 'frap'

% % F127 dataset
%     pluronic = 'F127';
%     folder1 = 'F127_BSA_25C';
%     folder2 ='trial_1';
%     conc = [17.5;20;22.5;25;27.5;30]; 
%     temperature = '25C'; % later could have this pickup from the filename?
%     dThresh = 0.7;
%     dy = 300;   
%     dx = 300;

% % F87 dataset
%     pluronic = 'F87';
%     folder1 = 'F87_BSA_25C';
%     folder2 ='trial_1';
%     conc = [25;30;35;37.5;40;42.5]; 
%     temperature = '25C'; % later could have this pickup from the filename?
%     dThresh = 0.9; %0.7 didnt work
%     dy = 350;   
%     dx = 350;
    
% P123 dataset trial 1
    pluronic = 'P123';
    folder1 = 'P123_BSA_25C';
    folder2 ='trial_1';
    conc = [20;25;27.5;30;32.5;35]; % later could have this pickup from a text file?
    temperature = '25C'; % later could have this pickup from the filename?
    dThresh = 1; %0.7 didnt work
    dy = 350;   
    dx = 350;
    
% frame to analyze as first bleached frame (frame 1)
    t1 = 2;
   
% folder where data is sorted data
    global boxfolder outputfolder plotfolder
    boxfolder = 'C:\Users\user\Box\Sorted Data FRAP\Sorted Data';
    outputfolder = fullfile(boxfolder,'z_outputs');
    plotfolder = fullfile(boxfolder,folder1,folder2,'z_plots');
    
%initialize the structure to store all data in 
alldata = struct();
fits = struct();
id = struct();
images = struct();

% get the subfolder names from the folders instead of naming by number
    folder3 = 'frap';
    d = dir(fullfile(boxfolder,folder1,folder2,folder3));
% remove all files (isdir property is 0)
    dfolders = d([d(:).isdir]);
% remove '.' and '..'from the list
	dfolders = dfolders(~ismember({dfolders(:).name},{'.','..'}));
    pos_struct = rmfield(dfolders,{'bytes','date','isdir','datenum','folder'});  % cleaning up the data structure
    npos = length(pos_struct); % number of positions
    
% re-order positions so they're arranged from lowest number to highest
% instead of alphabetical
positions = zeros(npos,1);
for p = 1:npos
    A = regexp(pos_struct(p).name,'\d*','Match');
    positions(p) = str2double(A{1});
end 
    %resort the position order
    [sortedpos,sortorder] = sort(positions);

% form position data structure with ID info 
for p = 1:npos
    % First we make the folder name for the position we want
    % add the positions in numerical order:
    idx = sortorder(p);
    position = pos_struct(idx).name;
    id.(position).('plur') = pluronic;
    id.(position).('plwt') = conc(p);
    id.(position).('prot') = 'BSA';
    id.(position).('protc') = '4p2';
    id.(position).('temp') = temperature;
end

%% Reading in data V2
close all
clc

% analyze frap data first.
    % timing
    tic;
    count = 1;
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
    
     % Analyze first bleached frame
    f1 = fullfile(data(t1).folder,data(t1).name);
    im = imread(f1);

    im1=double(im);
    % for plotting images
    imlb = min(im(:));
    imub = mean2(im);
    %threshold = (max(im1(:)) - min(im1(:)));
    T = graythresh(im)*dThresh; % works better with 0.7 for some

    % First we do some pre-analysis on the images to find the bleached circle
    bw1 = imbinarize(im,T); % make image black and white
    bwf = imcomplement(bw1);% flip black and white    
    bw2 = bwareaopen(bwf,200); % remove specks and dots smaller than lim2 pixels in area
    bw3 = imfill(bw2,8,'holes'); % fill in the holes so we can find the area
    % filling in the boundaries in black and white image to make solid shapes
    [B,L] = bwboundaries(bw3,'noholes');

    % These objects are now analyzed by regionprops(). 
    % The output is a table, called stats, that has the information for each object
    stats = regionprops('table', L, 'Centroid', 'EquivDiameter','Circularity','Solidity','Area','PixelList'); 
    main_idx = find(stats.Solidity == max(stats.Solidity));

    circularity = stats.Circularity(main_idx);
    idx = stats.PixelList(main_idx);
    idx = cell2mat(idx);
    idx_bleach = sub2ind(size(im1),idx(:,2),idx(:,1));
    center = stats.Centroid(main_idx,:); 
    radius = (stats.EquivDiameter(main_idx)/2);

    % reference region dx and dy defined above 
    
    % convert to pixel # list instead of coords
    idx_ref = sub2ind(size(im1),idx(:,2)+dy,idx(:,1)+dx);  
        
       % %% show the first image after photobleaching
        % modify bw3 to show reference region 
        bw4 = bw1;
        bw4(idx_ref) = 0;
        
        %figure(f1);
        fig = figure('name',position,'visible','on');
        set(fig, 'WindowStyle', 'Docked');  %figure will dock instead of free float
        
        subplot(1,2,1);
            hold on 
            imshow(im,[imlb,imub],'Border','tight','InitialMagnification', 'fit');
            viscircles(center,radius,'linewidth',0.2,'color','g');
            % Increase the size and position of the image 
            pos1 = get(gca,'Position'); pos1(1) = 0; pos1(2) = 0.2; pos1(3) = 0.5; 
            set(gca,'Position',pos1); 
            hold off 
        subplot(1,2,2)
            hold on
            % first plot the black and white image
            imshow(bw4,'Border','tight','InitialMagnification', 'fit')
            h= viscircles(center,radius,'linewidth',0.2,'color','g');
            % Increase the size and position of the image 
            pos2 = get(gca,'Position'); pos2(1)= 0.5; pos2(2) = 0.2; pos2(3) = 0.5;
            set(gca,'Position',pos2);
            hold off 
        
            % saving the figure as a high quality png 
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8 6];             % define location to save the images 
            a = fieldnames(id);
            pp = a{1};
            struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
            plot_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','_frames'];
            plot_path = fullfile(plotfolder,[plot_name,'.png']); % can change saved name here
            print(fig,plot_path, '-painters', '-dpng', '-r600')

            % save the circle and reference regions into fits 
            fits.(position).('radius') = radius;
            fits.(position).('center') = center;
            fits.(position).('idx_ref') = idx_ref;
            fits.(position).('idx_bleach') = idx_bleach;
            fits.(position).('imbds') = [imlb,imub]; 
            fits.(position).('images') = images;
end  
disp("bleached circles found in " + string(round(toc)) +" s");
% now repeating for pre-bleached images


%% 
for field = fieldnames(id)'
    position = field{1};
    % timing
    tic
    % Make a list of each image name in our position folder
    folder3 = 'prebleach';
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    % Analyze the pre-bleach image
    f1 = fullfile(data(1).folder,data(1).name);
    im0=double(imread(f1));    
    idx_ref = fits.(position).idx_ref;
    idx_bleach = fits.(position).idx_bleach;
    radius = fits.(position).radius;
    center = fits.(position).center; 
    
    % intensity of the reference region before bleaching
    refI_t0 = mean(im0(idx_ref));
    % intensity of the spot before bleaching 
    I_t0 = mean(im0(idx_bleach));
        
%         %show the first image after photobleaching
%         fig2 = figure('visible','on');
%         set(fig2, 'WindowStyle', 'Docked');  %figure will dock instead of free float
%             hold on 
%             imshow(im0,[min(im0(:)),mean2(im0)],'Border','tight','InitialMagnification', 'fit');
%             viscircles(center,radius,'linewidth',0.2,'color','g');
%             % Increase the size and position of the image 
%             hold off 
   
        % save the circle and reference regions into fits 
        % I0, ref0, what else?
            fits.(position).('refI_t0') = refI_t0;
            fits.(position).('I_t0') = I_t0;
   
end  
%%
for field = fieldnames(id)'
    position = field{1};
    tic % timing
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

        if data(1).pixel_size == 0
            fits.(position).('pixel_size') = 0.323; %for 20x in case it was not found 
        else
            fits.(position).('pixel_size') = data(1).pixel_size;
        end
        % add pixel intensity information
        f1 = fullfile(data(t).folder,data(t).name);
        im = double(imread(f1));

        I_ti = mean(im(idx_bleach));
        refI_ti = mean(im(idx_ref));

        normalized_i = (I_ti./I_t0).*refI_t0/(refI_ti);
        normalized_i = (I_ti./I_t0);

        % adding to data structure
        data(t).('I_ti') = I_ti;
        data(t).('refI_ti') = refI_ti;
        data(t).('IN_ti') = normalized_i;
        alldata.(position) = data;
    end
    
    % add end frame to plots 
    radius = fits.(position).radius;
    center = fits.(position).center;
    imbds = fits.(position).imbds;
    figure(p);
        subplot(1,2,2);
        hold on 
        imshow(im,imbds,'Border','tight','InitialMagnification', 'fit');
        viscircles(center,radius,'linewidth',0.2,'color','g');
        % Increase the size and position of the image 
        pos2 = get(gca,'Position'); pos2(1)= 0.5; pos2(2) = 0.2; pos2(3) = 0.5;
        set(gca,'Position',pos2);
        hold off 
    a = round(toc);
    disp([position,' loaded in ',num2str(a),' s.'])
end

disp("Data Loaded and intensity measured")

%% Adding some additional things to the data structure
% 1. elapsed time from first "initial" frame
% 2. fitting the data to equation outlined by Cheng 
% fit notes: Ifit at end of fit equation ((1-fm)*Ifit) is basically the
% normalized intensity of first bleached spot after laser beam

for field = fieldnames(alldata)'
    tic
    position = field{1}; % use{} bcuz field is a cell array

    % 1. adding the elapsed time for each frame, based on the first recieved
    % time in the data structure. In this case, the first received time is the
    % time of the first frame in "initial" folder
    time1 = alldata.(position)(1).r_time;
    time1 = datevec(time1);
    for i = 1:length(alldata.(position))
        t = alldata.(position)(i).r_time;
        t = datevec(t);
        dt = etime(t,time1); % calculate elapsed time from the first image in minutes
        % add to the data structure
        alldata.(position)(i).('dt') = dt;
    end
    
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
    ft = fittype('FRAPfit(x,f,k,tau)');
        options = fitoptions(ft); 
        options.StartPoint = [0.75 1 10];
        options.Lower = [0.15,0,0];
        options.Upper = [1,10,5000];
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
    Dl = (ri.^2)./(4*tau + max(dtau));
    Dr = (ri.^2)./(4*tau - min(dtau));
    Dneg = D-Dl;
    Dpos = Dr-D;
    errD = [Dneg Dpos];
    % creating new data structure for fit data
    fits.(position).('fit_info') = f;
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
    % timing 
    b = round(toc);
    disp([position,' data fit in ',num2str(b),' s.'])
    disp(["f= " + string(f.f) + ",k= "+string(f.k)+",tau= "+string(round(f.tau))]);
end

disp('time added and fits performed')
%% Exporting the data structures for future use.
% can also save individual fields if we want to get crazy
% first specify output folder (where we want to save the data
% maybe would be good to do conc instead of pos# here? 
cd(outputfolder);

% save the data structures for later plotting/analysis
a = fieldnames(id);
pp = a{1};

struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
alldata_name = [struct_name,'_','data','.mat'];
fits_name = [struct_name,'_','fits','.mat'];
id_name = [struct_name,'_','id','.mat'];
save(alldata_name,'alldata');
save(fits_name,'fits');
save(id_name,'id');
cd(mainfolder);
%% plotting the analyzed data 
C = jet;
for field = fieldnames(alldata)'
    position = field{1}; % use{} bcuz field is a cell array

    fig = figure('name',position,'visible','on');
    set(fig, 'WindowStyle', 'Docked');
    
%     subplot(2,1,1)
    title(position)
    % choosing data to plot
    t = [alldata.(position).dt]';
    t = t(t1:end);
    I = [alldata.(position).IN_ti]';
    I = I(t1:end);
    % fit info
    t_fit = [fits.(position).t_fit]';
    I_fit = [fits.(position).I_fit]';
    xc = fits.(position).xc;
    shapec = fits.(position).shapec;
    
    % adding photobleaching effect
    refI_ti = [alldata.(position).refI_ti]'./[fits.(position).refI_t0]';
    refI_ti = refI_ti(t1:end);
    
    % plotting and marker choice, linestyle, etc
    fill(xc',shapec,'r','FaceAlpha',0.1,'LineStyle','--','linewidth',1);
    hold on
    plot(t_fit,I_fit,'--')
    plot(t,I,'x','color',C(40,:))
    

%     plot(t,refI_ti,'d')
    axis([0,max(t),0,1])
    xlabel('time');
%     subplot(2,1,2);
%         center = fits.(position).center;
%         radius = fits.(position).radius;
%         hold on 
%         image = alldata.(position)(end).image;
%         imshow(image,[min(image(:)),mean2(image)],'Border','tight','InitialMagnification', 'fit');
%         viscircles(center,radius,'linewidth',0.2,'color','g');
%         hold off 
end
    
%% diffusivity data vs conc

fig = figure('name','All','visible','on');
set(fig, 'WindowStyle', 'Docked'); 

% need for loop to unpack
hold on
i = 0;
for field = fieldnames(alldata)'
    position = field{1};
    i = i+1;
    D = fits.(position).D/55.1;
    c = id.(position).plwt;
    Dall(i) = D; % from cheng to normalize by free soln BSA D
    Dneg = fits.(position).errD(1)/55.1;
    Dpos = fits.(position).errD(2)/55.1;
    plot(c,D,'db')
    errorbar(c,D,Dneg,Dpos,'ob')
end
    axis([0.8*min(conc) 1.2*max(conc) 0.8*min(Dall) 1.2*max(Dall)])
    set(gca,'YScale','log');
    xlabel('F127 wt%')
    ylabel('D/D_0 [\mum^2s^{-1}]')
