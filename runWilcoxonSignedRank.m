function wilcoxon_results = runWilcoxonSignedRank(pre, post, alpha, tail_dir)
    % [description]
    %
    % Input:
    %   pre  - (64, 64, n) connectivity matrix for "pre" condition
    %   post - (64, 64, n) connectivity matrix for "post" condition
    %   alpha - significance level (default 0.05)
    %   tail_dir - direction of tail for signrank test ('tail', 'right' is
    %   checking whether pre-post > 0; ==> if yes, then there was a decrease)
    %
    % Output:
    % 

    [n_channels, ~, n_participants] = size(pre);

    % Initialize storage
    p_values = []; % Store p-values for all (i, j) pairs
    w_stats = [];

    % Check whether matrices are symmetric --> FOR NOW DO EVERYTHING AS IF
    % NOT SYMMETRIC; MAYBE INTRODUCE THIS LATER
    % if issymmetric(pre) && issymmetric(post)
    
    % Iterate through upper triangle of matrix (including diagonal)
    for i = 1:n_channels
        for j = 1:n_channels
            % Extract paired data across participants
            pre_values = squeeze(pre(i, j, :));  % Vector of size (n_participants, 1)
            post_values = squeeze(post(i, j, :)); % Vector of size (n_participants, 1)
            % Perform Wilcoxon signed-rank test
            try 
                [p,h,stats] = signrank(pre_values, post_values, 'tail', tail_dir); % Two-tailed test --> 'tail', 'right' is checking whether pre-post > 0; ==> if yes, then there was a decrease
                % disp(stats);
                p_values = [p_values; p]; % Append p-value to the list
                w_stats = [w_stats; stats.signedrank];
            catch ME
                if ME.message == "No data remaining after removal of NaNs."
                    p = NaN;
                    p_values = [p_values; p];
                    w_stats = [w_stats; p];
                end
            end 
        end
    end

    % Apply FDR correction
    % [fdr_corrected_pvals, rejected] = mafdr(p_values, 'BHFDR', true);
    fdr_corrected_pvals = mafdr(p_values, 'BHFDR', true);
    rejected = fdr_corrected_pvals < alpha; % Logical array for significant tests
    
    % Identify significant pairs
    significant_indices = find(rejected); % Indices of significant pairs
    
    % Map indices back to (i, j) electrode pairs
    significant_pairs = [];
    w_stat_vals = zeros(n_channels, n_channels);
    count = 0;
    for i = 1:n_channels
        for j = 1:n_channels
            count = count + 1;
            if i == j
                w_stat_vals(i, j) = 0;
            else
                w_stat_vals(i, j) = w_stats(count);
            end
            if ismember(count, significant_indices)
                significant_pairs = [significant_pairs; i, j]; % Append (i, j) pair
            end
        end
    end

    % normalize W statistic to be between [-1 1] --> interpret this as
    % "effect size" (Barnett et al., 2020)
    r_total = n_participants * (n_participants + 1) / 2;
    norm_factor = r_total / 2;
    w_normalized = (w_stat_vals - norm_factor) / norm_factor;

    % set up structure to return
    wilcoxon_results = struct();
    wilcoxon_results.w_stat_vals = w_stat_vals;
    wilcoxon_results.w_normalized = w_normalized;
    wilcoxon_results.significant_pairs = significant_pairs; 
    wilcoxon_results.orig_h = h; 
    wilcoxon_results.orig_p_values = p_values;
    wilcoxon_results.corrected_p_values = fdr_corrected_pvals; 

    % % Display result
    % disp('Normalized Wilcoxon W-statistics: (effect size - from Barnett et al., 2020)');
    % figure;imagesc(W_normalized);
    % xlabel('To Node');
    % ylabel('From Node');
    % % xticks(data_load.data_interp_avgref.label);
    % % cbh = colorbar();
    % N = 256; % number of colorsdd
    % % cmap = [linspace(0,1,N).' linspace(0,1,N).' ones(N,1)]; % decreasing R and G; B = 1
    % cmap = brewermap(N, '-RdBu');
    % cbh = colormap(cmap);
    % colorbar;
 
end
