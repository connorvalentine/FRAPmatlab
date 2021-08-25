

function [alldata] = fun_frap_frames_drift_tracking(id,fits,alldata)
global boxfolder folder1 folder2 makemovies moviefolder
folder3 = 'frap';
counter = 0;

for field = fieldnames(id)' % iterate through the position list in id structure
    tic
    position = field{1};
%     position = 'Pos9';
    folder3 = 'frap';

    % make one figure that will be updated at every 10th timepoint just to
    % watch it run
    fig2 = figure('name',position,'visible','on');
    set(fig2, 'WindowStyle', 'Docked');  %figure will dock instead of free float
    set(gcf, 'Color','white')
    set(gca, 'nextplot','replacechildren', 'Visible','off');

    % Make a list of each image name in our position folder
    list = dir(fullfile(boxfolder,folder1,folder2,folder3,position,'*.tif')); % lists all files with .tif ending in the position folder
    data = rmfield(list,{'bytes','date','isdir','datenum'});  % cleaning up the data structure
    n_images = length(list);
    
    if makemovies == 'y'
            if counter == 3
                counter = 0;
            else
            end
            counter = counter +1;
        a = fieldnames(id);
        pp = a{1};    
        struct_name = [id.(pp).plur,'_',id.(pp).prot,'_',id.(pp).temp,'_',folder2];
       	movie_name = [struct_name,'_',num2str(round(id.(position).plwt)),'wtp','movie',num2str(counter)];
        movie_path = fullfile(moviefolder,movie_name); % can change saved name here
        vidObj = VideoWriter([movie_path, '.avi']);
        vidObj.Quality = 100;
        vidObj.FrameRate = 10;
        open(vidObj);
    else 
    end
        
    %# create movie
    % iterate through the images in data structure
    for t = 1:n_images 
        % add time information to the data structure
        [data] = fun_time(data,t); 
       % add elapsed time from received time of first frap frame to data struct.
        time1 = data(1).r_time;
        time1 = datevec(time1);
        time_i = data(t).r_time;
        time_i  = datevec(time_i);
        dt = etime(time_i,time1); % calculate elapsed time from the first image in seconds
        data(t).('dt') = dt; %save to data structure
        try 
            % calculate pixel intensity information
            f1 = fullfile(data(t).folder,data(t).name);
            imi = double(imread(f1));
            if t == 1
                previous_center = fits.(position).center;
                previous_radius = fits.(position).radius;
                ref_0_i = fits.(position).ref_0;
                ref_i = fits.(position).ref_1;
                I_ti = fits.(position).I_t1;
            else 
                [center,ref_i,I_ti,ref_0_i,drift_flag] = fun_radius_finder_i(position,imi,fits,previous_center,previous_norm_ratio,t,dt);
                if makemovies == 'y'
                    writeVideo(vidObj, getframe(gca));
                else
                end
                previous_center = center;
                
                if drift_flag == 1
                    disp('center has drifted too close to edge of frame')
                break
                % center has drifted too far, stop analyzing this frame
                else 
                end
            end
            
            if ref_i == 0
                ref_i = 0.001; % avoids a divide by zero error
                disp(position)
                disp('reference region mistake at frame #' +string(t))
            else
                norm_ratio =  ref_0_i/ref_i;
                previous_norm_ratio = norm_ratio; % for feeding into next loop
            end
            
            % normalized by prebleach Intensity(I_t0)
            % double normalized by the intensity ratio in the reference region
            % (refI_t0) divdided by the ref region ROI at time ti refI_ti
            normalized_i = (I_ti./fits.(position).I_t0).*(norm_ratio); 

            % adding to data structure

            data(t).('I_ti') = I_ti;
            data(t).('ref_i') = ref_i;
            data(t).('ref_ratio') = norm_ratio;
            data(t).('IN_ti') = normalized_i;
            
%             if makemovies == 'y'
%                 close(gcf)
%                 %# save as AVI file, and open it using system video player
%                 close(vidObj);
%             else 
%             end

        catch e
            warning('weird thing happened')
            fprintf(1,'There was an error! The message was:\n%s',e.message);
            data(t).('dt') = dt; %save to data structure
            data(t).('I_ti') = data(t-1).I_ti;
            data(t).('ref_i') = data(t-1).ref_i;
            data(t).('ref_ratio') = data(t-1).ref_ratio;
            data(t).('IN_ti') = data(t-1).IN_ti;
        end     
    end

    % add the completed data structure to alldata
    alldata.(position) = data;
    disp([position, 'loaded in ' + string(round(toc)) + " s."])
end