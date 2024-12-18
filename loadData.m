function data_out = loadData(filename, CAR_bool)
    % Load the EEG data, preprocess with CAR if needed
    %
    % Input:
    %   filename (string) - filename of the data 
    %   CAR_bool (boolean) - true: run rereferencing using CAR; false: do 
    %   not run CAR
    %
    % Output:
    %   data_out (struct) - cleaned data in the structure necessary for
    %   FieldTrip
    % 

    % load data 
    data_load = load(filename);
    data_before_car = data_load.data_clean;

    % ignore 'resp' channel if it exists
    if ismember("Resp", data_before_car.label)
        num_clean_channels = length(data_before_car.label) - 1; % ignore the "Resp" channel
        data_before_car.label = data_before_car.label(1:num_clean_channels); 
        data_before_car.trial = cellfun(@(m)m(1:num_clean_channels, :), data_before_car.trial, 'UniformOutput',false);
    end
    
    if CAR_bool == true
        disp('running CAR...');
        % implement CAR 
        cfg = [];
        cfg.channel = 'all'; % this is the default
        cfg.reref = 'yes';
        cfg.refmethod = 'avg';
        cfg.refchannel = 'all';
        data_out = ft_preprocessing(cfg, data_before_car);
    else
        disp('skipping CAR...');
        data_out = data_before_car;
    end
end
