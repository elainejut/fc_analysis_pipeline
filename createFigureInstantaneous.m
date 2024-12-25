function createFigureInstantaneous(fc_pre_group_all_chan, fc_post_group_all_chan, mode, p_idx, ELECTRODE_ORGANIZATIONS, layout, freq_band, save_dir)
    % [description]
    %
    % Input:
    %   p_idx = value representing which participant; or 'avg' to sum
    %   across all participants? (not implemented yet)
    %
    % Output:
    %   
    % 

    % find max value of both pre and post matrices to create standard
    % z-range
    

    % CONNECTIVITY MATRIX
    % Original electrode order
    pre_matrix = fc_pre_group_all_chan(:, :, p_idx);
    post_matrix = fc_post_group_all_chan(:, :, p_idx);

    pre_max = max(pre_matrix, [], 'all');
    post_max = max(post_matrix, [], 'all');
    max_conn = max([pre_max, post_max]);

    % Reordered by electrode region
    reorderedMatrix = pre_matrix(ELECTRODE_ORGANIZATIONS.by_letter.idx, ELECTRODE_ORGANIZATIONS.by_letter.idx);
    reorderedLabels = ELECTRODE_ORGANIZATIONS.by_letter.label;
    f = figure('Visible','off');
    imagesc(reorderedMatrix, [0 max_conn]); % [0 1] for standardization between all 
    % title('Normalized W-statistic as effect size');
    xlabel('To Node');
    ylabel('From Node');
    axis square; % Makes the axes equal in size for a square plot
    % Add labels to rows and columns
    % set(gca,'fontsize', ) 
    xticks(1:length(reorderedLabels)); % Set x-axis ticks at integer positions
    yticks(1:length(reorderedLabels)); % Set y-axis ticks at integer positions
    xticklabels(reorderedLabels); % Set x-axis labels
    yticklabels(reorderedLabels); % Set y-axis labels
    xtickangle(90); % Rotates x-axis labels by 45 degrees for better readability
    ax = gca;
    ax.XAxis.FontSize = 6;
    ax.YAxis.FontSize = 6;
    ax.FontWeight = 'bold';
    N = 256; % number of colorsdd
    cmap = brewermap(N, '-RdBu');
    cbh = colormap(cmap);
    colorbar;
    saveas(f, sprintf("%s/fc_pre_group_all_chan_%03d.png", save_dir, p_idx));

    reorderedMatrix = post_matrix(ELECTRODE_ORGANIZATIONS.by_letter.idx, ELECTRODE_ORGANIZATIONS.by_letter.idx);
    reorderedLabels = ELECTRODE_ORGANIZATIONS.by_letter.label;
    f = figure('Visible','off');
    imagesc(reorderedMatrix, [0 max_conn]); % [0 1] for standardization between all 
    % title('Normalized W-statistic as effect size');
    xlabel('To Node');
    ylabel('From Node');
    axis square; % Makes the axes equal in size for a square plot
    % Add labels to rows and columns
    % set(gca,'fontsize', ) 
    xticks(1:length(reorderedLabels)); % Set x-axis ticks at integer positions
    yticks(1:length(reorderedLabels)); % Set y-axis ticks at integer positions
    xticklabels(reorderedLabels); % Set x-axis labels
    yticklabels(reorderedLabels); % Set y-axis labels
    xtickangle(90); % Rotates x-axis labels by 45 degrees for better readability
    ax = gca;
    ax.XAxis.FontSize = 6;
    ax.YAxis.FontSize = 6;
    ax.FontWeight = 'bold';
    N = 256; % number of colorsdd
    cmap = brewermap(N, '-RdBu');
    cbh = colormap(cmap);
    colorbar;
    saveas(f, sprintf("%s/fc_post_group_all_chan_%03d.png", save_dir, p_idx));

    % % CONNECTIVITY CIRCLE
    % % Reorganize electrode order for circle graph
    % newOrder = ELECTRODE_ORGANIZATIONS.by_region.idx;
    % reorderedMatrix = fc_wilcoxon_effect(newOrder, newOrder);
    % reorderedLabels = ELECTRODE_ORGANIZATIONS.by_region.label;
    % % Separate into "increase" and "decrease"
    % mask_increase = (reorderedMatrix > 0);
    % increase = zeros(size(reorderedMatrix));
    % increase(mask_increase) = reorderedMatrix(mask_increase);
    % mask_decrease = (reorderedMatrix < 0);
    % decrease = zeros(size(reorderedMatrix));
    % decrease(mask_decrease) = reorderedMatrix(mask_decrease);
    % % Normalize for visualization (increase matrix)
    % norm_increase = abs(increase) / max(abs(increase(:))); % Normalize matrix values
    % norm_increase_thresh = norm_increase>0.9;
    % N = length(norm_increase_thresh);
    % cmap = brewermap(N, 'Reds');
    % f = figure('Visible','off');
    % circularGraph(norm_increase_thresh, 'Colormap', cmap, 'Label', reorderedLabels);
    % saveas(f, sprintf("%s/conn_circle_increase.png", save_dir)); 
    % 
    % % Normalize for visualization (decrease matrix)
    % norm_decrease = abs(decrease) / max(abs(decrease(:))); % Normalize matrix values
    % norm_decrease_thresh = norm_decrease>0.9;
    % % Visualize with CircularGraph toolbox
    % N = length(norm_decrease_thresh);
    % cmap = brewermap(N, 'Blues');
    % f = figure('Visible','off');
    % circularGraph(norm_decrease_thresh, 'Colormap', cmap, 'Label', reorderedLabels);
    % saveas(f, sprintf("%s/conn_circle_decrease.png", save_dir)); 
    % 
    % % TOPOGRAPHIC BRAIN
    % if mode == "undirected"
    %     assert(all(abs(sum(pre_matrix, 2, "omitnan")-sum(pre_matrix, 1, "omitnan").') < 0.00001)); % checking that matrix is symmetric so sum across rows and columns should be equal
    %     per_node_pre = sum(pre_matrix, 2, "omitnan"); % Mean of rows (averaging across all pairs for each electrode)    
    %     per_node_post = sum(post_matrix, 2, "omitnan"); % Mean of rows (averaging across all pairs for each electrode)
    %     pre_max = max(per_node_pre);
    %     post_max = max(per_node_post);
    %     max_conn = max([pre_max, post_max]);
    % 
    %     % pre_matrix
    %     topo_data = [];
    %     topo_data.label = layout.label(1:64,:);           % Electrode labels
    %     topo_data.avg = per_node_pre';       % Connectivity values (1x64)
    %     topo_data.time = 1;                       % Dummy time point
    %     topo_data.dimord = 'chan_time';           % Specify data dimension
    %     f = figure('Visible','off');
    %     cfg = [];
    %     cfg.layout = layout;                 % Use loaded layout
    %     cfg.parameter = 'avg';               % Specify data parameter
    %     cfg.zlim = [0 max_conn];
    %     cfg.colorbar = 'yes';
    %     cfg.comment = sprintf("freq=%s", freq_band);
    %     cfg.figure = gca;
    %     ft_topoplotTFR(cfg, topo_data);
    %     N = 256; % number of colorsdd
    %     cmap = brewermap(N, '-RdBu');
    %     cbh = colormap(cmap);
    %     colorbar;
    %     saveas(f, sprintf("%s/fc_pre_topoplot_%03d.png", save_dir, p_idx)); 
    % 
    %     % post_matrix
    %     topo_data = [];
    %     topo_data.label = layout.label(1:64,:);           % Electrode labels
    %     topo_data.avg = per_node_post';       % Connectivity values (1x64)
    %     topo_data.time = 1;                       % Dummy time point
    %     topo_data.dimord = 'chan_time';           % Specify data dimension
    %     f = figure('Visible','off');
    %     cfg = [];
    %     cfg.layout = layout;                 % Use loaded layout
    %     cfg.parameter = 'avg';               % Specify data parameter
    %     cfg.zlim = [0 max_conn];
    %     cfg.colorbar = 'yes';
    %     cfg.comment = sprintf("freq=%s", freq_band);
    %     cfg.figure = gca;
    %     ft_topoplotTFR(cfg, topo_data);
    %     N = 256; % number of colorsdd
    %     cmap = brewermap(N, '-RdBu');
    %     cbh = colormap(cmap);
    %     colorbar;
    %     saveas(f, sprintf("%s/fc_post_topoplot_%03d.png", save_dir, p_idx)); 

    % elseif mode == "directed"
    %     per_node_pre_out = sum(pre_matrix, 2, "omitnan"); % Mean of rows (averaging across all pairs for each electrode)    
    %     per_node_pre_in = sum(pre_matrix, 1, "omitnan"); % Mean of rows (averaging across all pairs for each electrode)    
    %     per_node_post_out = sum(post_matrix, 2, "omitnan"); % Mean of rows (averaging across all pairs for each electrode)    
    %     per_node_post_in = sum(post_matrix, 1, "omitnan"); % Mean of rows (averaging across all pairs for each electrode)    
    %     pre_out_max = max(per_node_pre_out);
    %     pre_in_max = max(per_node_pre_in);
    %     post_out_max = max(per_node_post_out);        
    %     post_in_max = max(per_node_post_in);
    %     max_conn = max([pre_out_max, pre_in_max, post_out_max, post_in_max]);
    % 
    %     % pre_matrix 
    %     % outbound
    %     topo_data = [];
    %     topo_data.label = layout.label(1:64,:);           % Electrode labels
    %     topo_data.avg = per_node_pre_out';       % Connectivity values (1x64)
    %     topo_data.time = 1;                       % Dummy time point
    %     topo_data.dimord = 'chan_time';           % Specify data dimension
    %     f = figure('Visible','off');
    %     cfg = [];
    %     cfg.layout = layout;                 % Use loaded layout
    %     cfg.parameter = 'avg';               % Specify data parameter
    %     cfg.zlim = [0 max_conn];
    %     cfg.colorbar = 'yes';
    %     cfg.comment = sprintf("freq=%s", freq_band);
    %     cfg.figure = gca;
    %     ft_topoplotTFR(cfg, topo_data);
    %     N = 256; % number of colorsdd
    %     cmap = brewermap(N, '-RdBu');
    %     cbh = colormap(cmap);
    %     colorbar;
    %     saveas(f, sprintf("%s/fc_pre_topoplot_outbound_%03d.png", save_dir, p_idx)); 
    %     % inbound
    %     topo_data = [];
    %     topo_data.label = layout.label(1:64,:);           % Electrode labels
    %     topo_data.avg = per_node_pre_in';       % Connectivity values (1x64)
    %     topo_data.time = 1;                       % Dummy time point
    %     topo_data.dimord = 'chan_time';           % Specify data dimension
    %     f = figure('Visible','off');
    %     cfg = [];
    %     cfg.layout = layout;                 % Use loaded layout
    %     cfg.parameter = 'avg';               % Specify data parameter
    %     cfg.zlim = [0 max_conn];
    %     cfg.colorbar = 'yes';
    %     cfg.comment = sprintf("freq=%s", freq_band);
    %     cfg.figure = gca;
    %     ft_topoplotTFR(cfg, topo_data);
    %     N = 256; % number of colorsdd
    %     cmap = brewermap(N, '-RdBu');
    %     cbh = colormap(cmap);
    %     colorbar;
    %     saveas(f, sprintf("%s/fc_pre_topoplot_inbound_%03d.png", save_dir, p_idx)); 
    % 
    %     % post_matrix
    %     % outbound
    %     topo_data = [];
    %     topo_data.label = layout.label(1:64,:);           % Electrode labels
    %     topo_data.avg = per_node_post_out';       % Connectivity values (1x64)
    %     topo_data.time = 1;                       % Dummy time point
    %     topo_data.dimord = 'chan_time';           % Specify data dimension
    %     f = figure('Visible','off');
    %     cfg = [];
    %     cfg.layout = layout;                 % Use loaded layout
    %     cfg.parameter = 'avg';               % Specify data parameter
    %     cfg.zlim = [0 max_conn];
    %     cfg.colorbar = 'yes';
    %     cfg.comment = sprintf("freq=%s", freq_band);
    %     cfg.figure = gca;
    %     ft_topoplotTFR(cfg, topo_data);
    %     N = 256; % number of colorsdd
    %     cmap = brewermap(N, '-RdBu');
    %     cbh = colormap(cmap);
    %     colorbar;
    %     saveas(f, sprintf("%s/fc_post_topoplot_outbound_%03d.png", save_dir, p_idx)); 
    %     % inbound
    %     topo_data = [];
    %     topo_data.label = layout.label(1:64,:);           % Electrode labels
    %     topo_data.avg = per_node_post_in';       % Connectivity values (1x64)
    %     topo_data.time = 1;                       % Dummy time point
    %     topo_data.dimord = 'chan_time';           % Specify data dimension
    %     f = figure('Visible','off');
    %     cfg = [];
    %     cfg.layout = layout;                 % Use loaded layout
    %     cfg.parameter = 'avg';               % Specify data parameter
    %     cfg.zlim = [0 max_conn];
    %     cfg.colorbar = 'yes';
    %     cfg.comment = sprintf("freq=%s", freq_band);
    %     cfg.figure = gca;
    %     ft_topoplotTFR(cfg, topo_data);
    %     N = 256; % number of colorsdd
    %     cmap = brewermap(N, '-RdBu');
    %     cbh = colormap(cmap);
    %     colorbar;
    %     saveas(f, sprintf("%s/fc_post_topoplot_inbound_%03d.png", save_dir, p_idx));
    % 


    % end
    % 
    % 
    % % NETWORK CONNECTIVITY ON TOPOGRAPHIC BRAIN ??
    % 

end
