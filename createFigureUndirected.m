function createFigureUndirected(wilcoxon_results, wilcoxon_results_out, wilcoxon_results_in, ELECTRODE_ORGANIZATIONS, layout, freq_band, save_dir)
    % [description]
    %
    % Input:
    %   
    %
    % Output:
    %   
    % 

    fc_wilcoxon_effect = wilcoxon_results.both.w_normalized;
    % fc_wilcoxon_effect_decrease = wilcoxon_results.right.w_normalized;
    % fc_wilcoxon_orig_sig_decrease = wilcoxon_results.right.orig_significant_pairs;
    % fc_wilcoxon_sig_decrease = wilcoxon_results.right.significant_pairs;
    % fc_wilcoxon_effect_increase = wilcoxon_results.left.w_normalized; % should be the same as fc_wilcoxon_effect_decrease
    % fc_wilcoxon_orig_sig_increase = wilcoxon_results.left.orig_significant_pairs;
    % fc_wilcoxon_sig_increase = wilcoxon_results.left.significant_pairs;
    fc_wilcoxon_orig_sig_both_01 = wilcoxon_results.both.orig_significant_pairs_01;
    fc_wilcoxon_orig_sig_both_001 = wilcoxon_results.both.orig_significant_pairs_001;
    fc_wilcoxon_sig_change = wilcoxon_results.both.significant_pairs;

    wilcoxon_out_effect = wilcoxon_results_out.right.w_normalized;
    wilcoxon_in_effect = wilcoxon_results_in.right.w_normalized;
    wilcoxon_sig_change = wilcoxon_results_out.both.significant_pairs; % note that this is in original electrode order
    % 
    % % CONNECTIVITY MATRIX
    % % Original electrode order
    % f = figure('Visible','off');
    % imagesc(fc_wilcoxon_effect_decrease, [-1 1]);
    % % title('Normalized W-statistic as effect size');
    % xlabel('To Node');
    % ylabel('From Node');
    % % Find significant connections and overlay dots
    % % rows = fc_wilcoxon_sig_decrease(:, 1); 
    % % cols = fc_wilcoxon_sig_decrease(:, 2); % Indices of significant connections
    % % scatter(rows, cols, 5, 'w', 'filled');   % Overlay white dots (size 50)
    % axis square;
    % ax = gca;
    % ax.XAxis.FontSize = 6;
    % ax.YAxis.FontSize = 6;
    % ax.FontWeight = 'bold';
    % N = 256; % number of colorsdd
    % cmap = brewermap(N, '-RdBu');
    % cbh = colormap(cmap);
    % colorbar;
    % saveas(f, sprintf("%s/orig_conn_mat.png", save_dir));

    % Reordered by electrode region 
    newOrder = ELECTRODE_ORGANIZATIONS.by_letter.idx;
    reorderedMatrix = fc_wilcoxon_effect(newOrder, newOrder);
    reorderedLabels = ELECTRODE_ORGANIZATIONS.by_letter.label;
    L = tril(reorderedMatrix); %To get the lower triangular part.
    % M = L.*(L > 0.5); % To introduce some artificial zeros into the lower triangular part.
    f = figure('Visible','off');
    imagesc(L, [-1 1]);
    % imagesc([1:size(L,2)]-0.5, [1:size(L,1)]-0.5, L, [-1 1]); % img is your image
    line(repmat([1.5:1:65.5],2,1),repmat([0;67], 1, 65), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.2); % vertical
    line(repmat([0;67], 1, 65), repmat([1.5:1:65.5],2,1), 'Color', [0.8 0.8 0.8], 'LineWidth', 0.2); % horizontal
    % grid();
    % xlabel('To Node');
    % ylabel('From Node');
    % xticks(1:length(reorderedLabels)); % Set x-axis ticks at integer positions
    % yticks(1:length(reorderedLabels)); % Set y-axis ticks at integer positions
    % xticklabels(reorderedLabels); % Set x-axis labels
    % yticklabels(reorderedLabels); % Set y-axis labels
    % xtickangle(90); % Rotates x-axis labels by 45 degrees for better readability
    % ax = gca;
    % ax.XAxis.FontSize = 6;
    % ax.YAxis.FontSize = 6;
    % ax.FontWeight = 'bold';
    N = 256; % number of colorsdd
    cmap = brewermap(N, '-RdBu');
    cbh = colormap(cmap);
    % colorbar;
    hold on;
    % line([0,65.5], [2.5,2.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [7.5, 7.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [16.5, 16.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [23.5, 23.5], 'Color', 'black', 'LineWidth', 3);
    line([0,65.5], [27.5, 27.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [35.5, 35.5], 'Color', 'black', 'LineWidth', 3);
    line([0,65.5], [41.5, 41.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [43.5, 43.5], 'Color', 'black', 'LineWidth', 3);
    line([0,65.5], [47.5, 47.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [43.5, 43.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [56.5, 56.5], 'Color', 'black', 'LineWidth', 3);
    line([0,65.5], [61.5, 61.5], 'Color', 'black', 'LineWidth', 3);
    % line([0,65.5], [64.5, 64.5], 'Color', 'black', 'LineWidth', 3);
    % line([7.5, 7.5], [0,65.5],'Color', 'black', 'LineWidth', 3);
    % line([2.5, 2.5], [0,65.5],'Color', 'black', 'LineWidth', 3);
    % line([16.5, 16.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % line([23.5, 23.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    line([27.5, 27.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % line([35.5, 35.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    line([41.5, 41.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % line([43.5, 43.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    line([47.5, 47.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % line([43.5, 43.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % line([56.5, 56.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    line([61.5, 61.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % line([64.5, 64.5], [0,65.5], 'Color', 'black', 'LineWidth', 3);
    % Create a mapping from original indices to new order positions
    indexMap = zeros(1, length(newOrder)); % Initialize mapping array
    indexMap(newOrder) = 1:length(newOrder); % Map original indices to new order positions
    % Apply the mapping to the sig_pairs matrix
    reordered_sig_pairs_01 = indexMap(fc_wilcoxon_orig_sig_both_01); % Replace indices with their new positions
    if ~isempty(reordered_sig_pairs_01)
        x_vals = reordered_sig_pairs_01(:, 2); 
        y_vals = reordered_sig_pairs_01(:, 1);
        scatter(x_vals, y_vals, 7, 'filled', 'MarkerFaceColor', [1, 0.65, 0]);
    end
    reordered_sig_pairs_001 = indexMap(fc_wilcoxon_orig_sig_both_001); % Replace indices with their new positions
    if ~isempty(reordered_sig_pairs_001)
        x_vals = reordered_sig_pairs_001(:, 2); 
        y_vals = reordered_sig_pairs_001(:, 1);
        scatter(x_vals, y_vals, 7, 'filled', 'MarkerFaceColor', [0.55, 1, 0]);
    end
    axis square;
    saveas(f, sprintf("%s/conn_mat_by_letter.png", save_dir));

    % % From node
    % reorderedMatrix = wilcoxon_out_effect(ELECTRODE_ORGANIZATIONS.by_letter.idx).';
    % reorderedLabels = ELECTRODE_ORGANIZATIONS.by_letter.label;
    % f = figure('Visible','off');
    % imagesc(reorderedMatrix, [-1 1]);
    % % xlabel('To Node');
    % ylabel('From Node');
    % % Add labels to rows and columns
    % % set(gca,'fontsize', ) 
    % % xticks(1:length(reorderedLabels)); % Set x-axis ticks at integer positions
    % yticks(1:length(reorderedLabels)); % Set y-axis ticks at integer positions
    % % xticklabels(reorderedLabels); % Set x-axis labels
    % yticklabels(reorderedLabels); % Set y-axis labels
    % % xtickangle(90); % Rotates x-axis labels by 45 degrees for better readability
    % ax = gca;
    % % ax.XAxis.FontSize = 6;
    % ax.YAxis.FontSize = 6;
    % ax.FontWeight = 'bold';
    % N = 256; % number of colorsdd
    % cmap = brewermap(N, '-RdBu');
    % cbh = colormap(cmap);
    % colorbar;
    % saveas(f, sprintf("%s/conn_mat_out.png", save_dir));

    % % To node
    % reorderedMatrix = wilcoxon_in_effect(ELECTRODE_ORGANIZATIONS.by_letter.idx);
    % reorderedLabels = ELECTRODE_ORGANIZATIONS.by_letter.label;
    % f = figure('Visible','off');
    % imagesc(reorderedMatrix, [-1 1]);
    % xlabel('To Node');
    % % ylabel('From Node');
    % xticks(1:length(reorderedLabels)); % Set x-axis ticks at integer positions
    % % yticks(1:length(reorderedLabels)); % Set y-axis ticks at integer positions
    % xticklabels(reorderedLabels); % Set x-axis labels
    % % yticklabels(reorderedLabels); % Set y-axis labels
    % xtickangle(90); % Rotates x-axis labels by 45 degrees for better readability
    % ax = gca;
    % ax.XAxis.FontSize = 6;
    % % ax.YAxis.FontSize = 6;
    % ax.FontWeight = 'bold';
    % N = 256; % number of colorsdd
    % cmap = brewermap(N, '-RdBu');
    % cbh = colormap(cmap);
    % colorbar;
    % saveas(f, sprintf("%s/conn_mat_in.png", save_dir));

    % CONNECTIVITY CIRCLE
    % Reorganize electrode order for circle graph
    newOrder = ELECTRODE_ORGANIZATIONS.by_region.idx;
    reorderedMatrix = fc_wilcoxon_effect; % (newOrder, newOrder);
    reorderedLabels = ELECTRODE_ORGANIZATIONS.by_region.label;
    % get significant connections (same as connectivity matrix above
    if ~isempty(fc_wilcoxon_orig_sig_both_01)
        x_vals = fc_wilcoxon_orig_sig_both_01(:, 1); 
        y_vals = fc_wilcoxon_orig_sig_both_01(:, 2);
    else
        x_vals = []; 
        y_vals = [];
    end
    mask_sig_01 = zeros(size(reorderedMatrix));
    linear_indices = sub2ind(size(reorderedMatrix), x_vals, y_vals);
    mask_sig_01(linear_indices) = 1;
    if ~isempty(fc_wilcoxon_orig_sig_both_001)
        x_vals = fc_wilcoxon_orig_sig_both_001(:, 1); 
        y_vals = fc_wilcoxon_orig_sig_both_001(:, 2);
    else
        x_vals = []; 
        y_vals = [];
    end
    mask_sig_001 = zeros(size(reorderedMatrix));
    linear_indices = sub2ind(size(reorderedMatrix), x_vals, y_vals);
    mask_sig_001(linear_indices) = 1;
    % Separate into "increase" and "decrease"
    mask_increase = (reorderedMatrix > 0);
    mask_increase_sig_01 = mask_sig_01 & mask_increase;
    mask_increase_sig_001 = mask_sig_001 & mask_increase;
    increase_01 = zeros(size(reorderedMatrix));
    increase_01(mask_increase_sig_01) = reorderedMatrix(mask_increase_sig_01);
    increase_001 = zeros(size(reorderedMatrix));
    increase_001(mask_increase_sig_001) = reorderedMatrix(mask_increase_sig_001);
    mask_decrease = (reorderedMatrix < 0);
    mask_decrease_sig_01 = mask_sig_01 & mask_decrease;
    mask_decrease_sig_001 = mask_sig_001 & mask_decrease;
    decrease_01 = zeros(size(reorderedMatrix));
    decrease_01(mask_decrease_sig_01) = reorderedMatrix(mask_decrease_sig_01);
    decrease_001 = zeros(size(reorderedMatrix));
    decrease_001(mask_decrease_sig_001) = reorderedMatrix(mask_decrease_sig_001);
    % Normalize for visualization (increase matrix)
    if max(abs(increase_01(:))) == 0
        norm_increase_01 = abs(increase_01);
    else
        norm_increase_01 = abs(increase_01) / max(abs(increase_01(:))); % Normalize matrix values
    end
    if max(abs(increase_001(:))) == 0
        norm_increase_001 = abs(increase_001);
    else
        norm_increase_001 = abs(increase_001) / max(abs(increase_001(:))); % Normalize matrix values    % norm_increase_thresh = norm_increase>0.9;
    end
    N_01 = length(norm_increase_01);
    % N_001 = length(norm_increase_001);
    cmap = brewermap(N_01, 'Reds');
    f = figure('Visible','off');
    % circularGraph(norm_increase_001, 'Colormap', cmap, 'Label', reorderedLabels);
    modifiedCircularGraph(norm_increase_01(newOrder, newOrder), 'SecondAdjacencyMatrix', norm_increase_001(newOrder, newOrder), 'Colormap', cmap, 'Label', reorderedLabels);
    saveas(f, sprintf("%s/conn_circle_increase.png", save_dir)); 
    % Normalize for visualization (decrease matrix)
    if max(abs(increase_01(:))) == 0
        norm_decrease_01 = abs(decrease_01);
    else
        norm_decrease_01 = abs(decrease_01) / max(abs(decrease_01(:))); % Normalize matrix values
    end
    if max(abs(increase_001(:))) == 0
        norm_decrease_001 = abs(decrease_001);
    else
        norm_decrease_001 = abs(decrease_001) / max(abs(decrease_001(:))); % Normalize matrix values    % norm_increase_thresh = norm_increase>0.9;
    end
    % norm_decrease_thresh = norm_decrease>0.9;
    % Visualize with CircularGraph toolbox
    N_01 = length(norm_decrease_01);
    % N_001 = length(norm_decrease_001);
    cmap = brewermap(N_01, 'Blues');
    f = figure('Visible','off');
    % circularGraph(norm_decrease_01, 'Colormap', cmap, 'Label', reorderedLabels);
    check1 = sum(norm_decrease_01, 'all');
    check2 = sum(norm_decrease_001, 'all');
    
    modifiedCircularGraph(norm_decrease_01(newOrder, newOrder), 'SecondAdjacencyMatrix', norm_decrease_001(newOrder, newOrder), 'Colormap', cmap, 'Label', reorderedLabels);
    saveas(f, sprintf("%s/conn_circle_decrease.png", save_dir)); 
    % 
    % TOPOGRAPHIC BRAIN
    connectivity_sum = wilcoxon_out_effect;    % same as wilcoxon_in_effect for undirected
    topo_data = [];
    topo_data.label = layout.label(1:65,:);           % Electrode labels
    topo_data.avg = connectivity_sum';       % Connectivity values (1x65)
    topo_data.time = 1;                       % Dummy time point
    topo_data.dimord = 'chan_time';           % Specify data dimension
    f = figure('Visible','off');
    cfg = [];
    cfg.layout = layout;                 % Use loaded layout
    cfg.parameter = 'avg';               % Specify data parameter
    cfg.zlim = [-1 1];
    % cfg.colorbar = 'yes';
    cfg.comment = 'no'; % sprintf("freq=%s", freq_band);
    cfg.maskparameter = 'mask';
    cfg.highlight = 'on';
    cfg.highlightchannel = wilcoxon_sig_change;
    cfg.highlightsymbol = '.';
    cfg.highlightsize = 28;
    cfg.highlightcolor = [1 1 0.3];
    cfg.marker = 'on';
    cfg.figure = gca;
    ft_topoplotTFR(cfg, topo_data);
    N = 256; % number of colorsdd
    cmap = brewermap(N, '-RdBu');
    cbh = colormap(cmap);
    % colorbar;
    saveas(f, sprintf("%s/topoplot.png", save_dir)); 


    % NETWORK CONNECTIVITY ON TOPOGRAPHIC BRAIN ??


end
