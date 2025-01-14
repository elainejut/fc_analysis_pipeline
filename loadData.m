function data_out = loadData(filename, reref_mode, downsample_bool)
    % Load the EEG data, preprocess with CAR if needed
    %
    % Input:
    %   filename (string) - filename of the data 
    %   reref_mode (string) - 'none': no rereferencing; 'car': common
    %   average referencing; 'linked': linked mastoid (use T7 and T8 for
    %   our data)
    %   downsample_bool (boolean) - true: downsample from 500 to 250 Hz;
    %   false: keep Fs=500 Hz
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
   
    if strcmp(reref_mode, 'car')
        disp('running CAR...');
        % implement CAR 
        cfg = [];
        cfg.channel = 'all'; % this is the default
        cfg.reref = 'yes';
        cfg.refmethod = 'avg';
        cfg.implicitref = 'FCz';
        cfg.refchannel = 'all';
        data_eeg = ft_preprocessing(cfg, data_before_car);
    elseif strcmp(reref_mode, 'linked')
        disp('running linked-ear rereferencing...')
        cfg = [];
        cfg.channel = 'all'; % this is the default
        cfg.reref = 'yes';
        cfg.refmethod = 'avg';
        cfg.implicitref = 'FCz'; % see below
        cfg.refchannel = {'T7', 'T8'}; % T7 and T8 as virtual ear lobes
        data_eeg = ft_preprocessing(cfg, data_before_car);
    else
        disp('skipping CAR and linked-ear rereferencing...');
        % cfg = [];
        % cfg.channel = 'all'; % this is the default
        % cfg.reref = 'no';
        % cfg.implicitref = 'FCz';
        % data_eeg = ft_preprocessing(cfg, data_before_car);
        data_eeg = data_before_car;
    end

    if downsample_bool == true
        disp('downsampling to 250 Hz...');
        cfg = [];
        cfg.resamplefs = 250;
        cfg.detrend = 'no';
        data_out = ft_resampledata(cfg, data_eeg);
    else
        disp('skipping downsampling (Fs=500 Hz)...');
        data_out = data_eeg;
    end

end
