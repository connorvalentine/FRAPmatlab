%% main analysis code to average/clean the data
% right now this function takes all the parameters in plot data and
% averages them.

% it only propagates error for Diffusion right now but that could be
% expanded

function [av,all] = fun_data_average(data,id,av)
for topfield = fieldnames(data)'
    pluronic = topfield{1};
    
    pd = data.(pluronic);
    disp(pd.D)
    pd_av = pd;
    
    % find the temperature 
    names = fieldnames(id.(pluronic))';
    temperature = id.(pluronic).(names{1}).temp;
    temperature = ['all' temperature]; % field names cannot start with a number
    c = pd.c;
    [C,ia,ic] = unique(c);
    % ic is the indices of the unique concentrations now. so if two rows
    % have the same number in ic, they are the same conc.
    
    % find the mean of each number here
    for fields = fieldnames(pd)'
        
        meand = [];
        stdd = [];
        field = fields{1};
        disp(field)
        data_to_average = pd.(field);
        for i = 1:length(C)
            % average all the values for the concentration, but use find to
            % ignore entries == 0
            temp = data_to_average(find(ic ==i));
            meand(i) = mean(temp(find(temp)));
            stdd(i) = std(temp(find(temp)));
            nums = temp;
        end
        
        % propagate uncertainty from the averaging of D as well.
        if strcmp(field,'D') 
            D_std = stdd; % standard dev from the mean in D
            disp(meand)
        else 
        end
        
        if strcmp(field,'Dneg') ||strcmp(field,'Dpos') 
            meand = sqrt((D_std.^2) + (meand.^2));            
        else
        end
        
        pd_av.(field) = meand;
    end 
    av.(char(temperature)).(pluronic) = pd_av;
    all.(char(temperature)).(pluronic) = pd;
end

av = av;
all = all;
end