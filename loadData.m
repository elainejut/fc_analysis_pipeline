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
