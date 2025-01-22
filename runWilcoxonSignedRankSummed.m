function wilcoxon_results = runWilcoxonSignedRankSummed(pre, post, alpha, tail_dir)
    % run wilcoxon signed rank test on electrode-level connectivity
    %
    % Input:
    %   pre  - (64, n) summed (either across rows or columns) connectivity matrix for "pre" condition 
    %   post - (64, n) summed (either across rows or columns) connectivity matrix for "post" condition
    %   alpha - significance level (default 0.05)
    %   tail_dir - direction of tail for signrank test (right now stick
    %   with "both")
    %
    % Output:
    %   wilcoxon_results (struct)
    % 

    [n_channels, n_participants] = size(pre);

    % Initialize storage
    p_values = []; % Store p-values for all (i, j) pairs
    w_stats = [];

    for i = 1:n_channels
        % Extract paired data across participants
        pre_values = pre(i, :).';  % Vector of size (n_participants, 1)
        post_values = post(i, :).'; % Vector of size (n_participants, 1)
        % Perform Wilcoxon signed-rank test
        try 
            [p,h,stats] = signrank(pre_values, post_values, 'tail', tail_dir); 
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

    % Apply FDR correction
    fdr_corrected_pvals = mafdr(p_values, 'BHFDR', true);
    rejected = fdr_corrected_pvals < alpha; % Logical array for significant tests
    orig_rejected_05 = p_values < 0.05; % alpha;
    orig_rejected_01 = p_values < 0.01; % alpha; % stricter threshold instead of alpha
    orig_rejected_001 = p_values < 0.001; % alpha; % stricter threshold instead of alpha

    % Identify significant pairs
    significant_indices = find(rejected); % Indices of significant pairs
    orig_significant_indices_05 = find(orig_rejected_05);
    orig_significant_indices_01 = find(orig_rejected_01);
    orig_significant_indices_001 = find(orig_rejected_001);

    % Map indices back to (i, j) electrode pairs
    significant_pairs = [];
    orig_significant_pairs_05 = [];
    orig_significant_pairs_01 = [];
    orig_significant_pairs_001 = [];
    w_stat_vals = []; % zeros(n_channels);
    count = 0;
    for i = 1:n_channels
        count = count + 1;
        w_stat_vals = [w_stat_vals, w_stats(count)];
        if ismember(count, significant_indices)
            significant_pairs = [significant_pairs; i]; % Append (i) electrode
        end
        if ismember(count, orig_significant_indices_05)
            orig_significant_pairs_05 = [orig_significant_pairs_05; i]; % Append (i) electrode
        end
        if ismember(count, orig_significant_indices_01)
            orig_significant_pairs_01 = [orig_significant_pairs_01; i]; % Append (i) electrode
        end
        if ismember(count, orig_significant_indices_001)
            orig_significant_pairs_01 = [orig_significant_pairs_001; i]; % Append (i) electrode
        end
    end

    % normalize W statistic to be between [-1 1] --> interpret this as
    % "effect size" (Barnett et al., 2020)
    total_rank_sum = n_participants * (n_participants + 1) / 2;
    % norm_factor = r_total / 2;
    % w_normalized_temp = (w_stat_vals - norm_factor) / norm_factor;
    w_normalized = 2 * (w_stat_vals / total_rank_sum) - 1;

    % set up structure to return
    wilcoxon_results = struct();
    wilcoxon_results.w_stat_vals = w_stat_vals;
    wilcoxon_results.w_normalized = w_normalized; 
    wilcoxon_results.orig_h = h; 
    wilcoxon_results.orig_p_values = p_values;
    wilcoxon_results.orig_significant_pairs_05 = orig_significant_pairs_05;
    wilcoxon_results.orig_significant_pairs_01 = orig_significant_pairs_01;
    wilcoxon_results.orig_significant_pairs_001 = orig_significant_pairs_001;
    wilcoxon_results.corrected_p_values = fdr_corrected_pvals; 
    wilcoxon_results.significant_pairs = significant_pairs; 
 
end
