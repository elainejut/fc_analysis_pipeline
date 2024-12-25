function createFigureDirected(fc_wilcoxon_effect, ELECTRODE_ORGANIZATIONS, layout, freq_band, save_dir)
    % [description]
    %
    % Input:
    %   
    %
    % Output:
    %   
    % 

    % CONNECTIVITY MATRIX
    % Original electrode order
    f = figure('Visible','off');
    imagesc(fc_wilcoxon_effect, [-1 1]);
    % title('Normalized W-statistic as effect size');
    xlabel('To Node');
    ylabel('From Node');
    axis square;
    ax = gca;
    ax.XAxis.FontSize = 6;
    ax.YAxis.FontSize = 6;
    ax.FontWeight = 'bold';
    N = 256; % number of colorsdd
    cmap = brewermap(N, '-RdBu');
    cbh = colormap(cmap);
    colorbar;
    saveas(f, sprintf("%s/orig_conn_mat.png", save_dir));

    % Reordered by electrode region
    reorderedMatrix = fc_wilcoxon_effect(ELECTRODE_ORGANIZATIONS.by_letter.idx, ELECTRODE_ORGANIZATIONS.by_letter.idx);
    reorderedLabels = ELECTRODE_ORGANIZATIONS.by_letter.label;
    f = figure('Visible','off');
    imagesc(reorderedMatrix, [-1 1]);
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
    saveas(f, sprintf("%s/conn_mat_by_letter.png", save_dir));

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
    % connectivity_mean = mean(fc_wilcoxon_effect, 2); % Mean of rows (averaging across all pairs for each electrode)    
    % topo_data = [];
    % topo_data.label = layout.label(1:64,:);           % Electrode labels
    % topo_data.avg = connectivity_mean';       % Connectivity values (1x64)
    % topo_data.time = 1;                       % Dummy time point
    % topo_data.dimord = 'chan_time';           % Specify data dimension
    % f = figure('Visible','off');
    % cfg = [];
    % cfg.layout = layout;                 % Use loaded layout
    % cfg.parameter = 'avg';               % Specify data parameter
    % cfg.zlim = [-1 1];
    % cfg.colorbar = 'yes';
    % cfg.comment = sprintf("freq=%s", freq_band);
    % cfg.figure = gca;
    % ft_topoplotTFR(cfg, topo_data);
    % N = 256; % number of colorsdd
    % cmap = brewermap(N, '-RdBu');
    % cbh = colormap(cmap);
    % colorbar;
    % saveas(f, sprintf("%s/topoplot.png", save_dir)); 
    % 
    % 
    % % NETWORK CONNECTIVITY ON TOPOGRAPHIC BRAIN ??
    % 

end
