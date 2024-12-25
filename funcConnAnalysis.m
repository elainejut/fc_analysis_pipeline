function connectivity_analysis_result = funcConnAnalysis(data, fc_method, freq_band, all_electrodes, electrode_organizations)
    % Perform connectivity analysis based on the method specified
    %
    % Input:
    %   data (struct) - data for analysis with FieldTrip analysis
    %   fc_method (string) - specify the type of FC analysis to do
    %       options: ['icoh', 'corr', 'granger']
    %   freq_band (string) - specify the frequency band to use
    %       options: ['theta', 'alpha', 'low_beta', 'high_beta', 'gamma']
    %
    % Output:
    %   connectivity_analysis_result (struct) - structure containing
    %   resulting connectivity matrices
    %       1. connectivity - initial connectivity output from FieldTrip 
    %       toolbox i
    %       2. conn_mat_orig - averaged conn_mat as output by FieldTrip
    %       (will not include the channels that are missing) 
    %       3. conn_mat_all_chan
    %       4. conn_mat_organized - conn_mat with electrode channels (ie.
    %       order of x-axis and y-axis)organized based on brain region 
    %       (organization should be consistent across all FC analyses)
    % 

    % create dictionary of frequency bands and their thresholds
    FREQ_DICT = containers.Map({'delta', 'theta', 'alpha', 'low_beta', 'high_beta', 'gamma', 'broadband'}, {[1 4], [4 8], [8 13], [13 20], [20 30], [30 45], [1 45]});
    disp(FREQ_DICT(freq_band));

    if strcmp(fc_method, 'icoh')
        disp('running functional connectivity analysis with imaginary part of coherence...');
        % frequency analysis --> maybe can also make this into its own function
        cfg           = [];
        cfg.method    = 'mtmfft';
        cfg.taper     = 'dpss';
        cfg.output    = 'fourier';
        cfg.foi    = FREQ_DICT(freq_band); % cfg.foi or cfg.foilim
        cfg.tapsmofrq = 2;
        cfg.channel   = 'all';
        cfg.keeptrials = 'yes';
        freq            = ft_freqanalysis(cfg, data);

        % disp("CHECKING FREQ BANDS: ");
        % disp(data);
        % disp(freq);
        % disp(freq.freq(1));
        % disp(freq.freq(end));

        % implement imaginary part of coherence functional connectivity 
        % analysis with FieldTrip
        cfg         = [];
        cfg.method  ='coh';
        cfg.complex = 'absimag';
        connectivity   = ft_connectivityanalysis(cfg, freq);
        
        % find the average connectivity matrix (averaged across each
        % frequency bin (dimension 3))
        conn_mat = mean(connectivity.cohspctrm, 3);
        conn_mat_diag = conn_mat - diag(diag(conn_mat)); % make sure diagonal is 0 % ERROR: TAKING THIS OUT BECAUSE DIAGONAL SHOULD BE 1 AND I SHOULDN'T NEED TO SET THIS


    elseif strcmp(fc_method, 'granger')
        disp('running functional connectivity analysis with granger causality...');
        % mvar analysis
        cfg = [];
        cfg.order = 8; % Define model order
        % cfg.keeptrial = 'yes';
        cfg.toolbox = 'biosig'; % Use the BSMART toolbox (FieldTrip default)
        mdata = ft_mvaranalysis(cfg, data);

        % freq analysis on mvar data
        cfg        = [];
        cfg.method = 'mvar';
        cfg.foi = FREQ_DICT(freq_band);
        mfreq      = ft_freqanalysis(cfg, mdata);

        % find granger causality
        cfg           = [];
        cfg.method    = 'granger';
        % cfg.granger.sfmethod = 'bivariate';  % Use bivariate Granger causality
        connectivity       = ft_connectivityanalysis(cfg, mfreq);

        % find the average connectivity matrix (averaged across each
        % frequency bin (dimension 3))
        conn_mat = mean(connectivity.grangerspctrm, 3);
        conn_mat_diag = conn_mat - diag(diag(conn_mat)); % make sure diagonal is 0
    
    elseif strcmp(fc_method, 'dtf')
        disp('running functional connectivity analysis with directed transfer function...');
        % mvar analysis
        cfg = [];
        cfg.order = 8; % Define model order
        % cfg.keeptrial = 'yes';
        cfg.toolbox = 'biosig'; % Use the BSMART toolbox (FieldTrip default)
        mdata = ft_mvaranalysis(cfg, data);

        % freq analysis on mvar data
        cfg        = [];
        cfg.method = 'mvar';
        cfg.foi = FREQ_DICT(freq_band);
        mfreq      = ft_freqanalysis(cfg, mdata);

        % find granger causality
        cfg           = [];
        cfg.method    = 'dtf';
        % cfg.granger.sfmethod = 'bivariate';  % Use bivariate Granger causality
        connectivity       = ft_connectivityanalysis(cfg, mfreq);

        % find the average connectivity matrix (averaged across each
        % frequency bin (dimension 3))
        conn_mat = mean(connectivity.dtfspctrm, 3);
        conn_mat_diag = conn_mat - diag(diag(conn_mat)); % make sure diagonal is 0
    
    elseif strcmp(fc_method, 'pdc')
        disp('running functional connectivity analysis with partial directed coherence...');
        % mvar analysis
        cfg = [];
        cfg.order = 8; % Define model order
        % cfg.keeptrial = 'yes';
        cfg.toolbox = 'biosig'; % Use the BSMART toolbox (FieldTrip default)
        mdata = ft_mvaranalysis(cfg, data);

        % freq analysis on mvar data
        cfg        = [];
        cfg.method = 'mvar';
        cfg.foi = FREQ_DICT(freq_band);
        mfreq      = ft_freqanalysis(cfg, mdata);

        % find granger causality
        cfg           = [];
        cfg.method    = 'pdc';
        % cfg.granger.sfmethod = 'bivariate';  % Use bivariate Granger causality
        connectivity       = ft_connectivityanalysis(cfg, mfreq);

        % find the average connectivity matrix (averaged across each
        % frequency bin (dimension 3))
        conn_mat = mean(connectivity.pdcspctrm, 3);
        conn_mat_diag = conn_mat - diag(diag(conn_mat)); % make sure diagonal is 0


    elseif strcmp(fc_method, 'mi')
        disp('running functional connectivity analysis with mutual information...');
        % analysis with FieldTrip
        cfg         = [];
        cfg.method  = 'mi';
        connectivity   = ft_connectivityanalysis(cfg, data);
        
        % find the connectivity matrix 
        conn_mat = connectivity.mi;
        conn_mat_diag = conn_mat - diag(diag(conn_mat)); % make sure diagonal is 0

    elseif strcmp(fc_method, 'corr')
        disp('running functional connectivity analysis with Pearson correlation...');
        % analysis with FieldTrip
        cfg         = [];
        cfg.method  ='corr';
        connectivity   = ft_connectivityanalysis(cfg, data);
        
        % find the connectivity matrix 
        conn_mat = connectivity.corr;
        conn_mat_diag = conn_mat - diag(diag(conn_mat)); % make sure diagonal is 0

    end
   
    % % extend from conn_mat 
    % labels_mask = ~ismember(connectivity.label, all_chan_labels);
    % % Apply the mask to rows
    % conn_mat(labels_mask, :) = NaN;
    % % Apply the mask to columns
    % conn_mat(:, labels_mask) = NaN;
    
    % disp("[funcConnAnalysis] debugging");
    % disp(connectivity);
    % disp(length(conn_mat));
    % disp(size(connectivity.cohspctrm));
    % disp(connectivity.label);

    conn_mat_all_chan = insertMissingElectrodes(all_electrodes, connectivity.label, conn_mat_diag);

    % disp(size(conn_mat_all_chan));
    % disp('conn_mat_all_chan:');
    % figure;imagesc(conn_mat_all_chan);
    % xlabel('To Node');
    % ylabel('From Node');
    % % xticks(data_load.data_interp_avgref.label);
    % cbh = colorbar();

    % conn_mat_by_letter = conn_mat_all_chan(electrode_organizations.by_letter.idx, electrode_organizations.by_letter.idx);
    % conn_mat_by_region = conn_mat_all_chan(electrode_organizations.by_region.idx, electrode_organizations.by_region.idx);

    connectivity_analysis_result = struct();
    connectivity_analysis_result.connectivity = connectivity;
    connectivity_analysis_result.conn_mat_orig = conn_mat;
    connectivity_analysis_result.conn_mat_all_chan = conn_mat_all_chan;
    % connectivity_analysis_result.conn_mat_by_letter = conn_mat_by_letter;
    % connectivity_analysis_result.conn_mat_by_region = conn_mat_by_region;


end
