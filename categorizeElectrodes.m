function electrode_organizations = categorizeElectrodes(all_electrodes)
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

    % BY REGION: (useful for connectivity circle graph visualization) 

    % Define electrode regions
    remaining_electrodes = all_electrodes;
    Fp_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'Fp', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, Fp_electrodes);
    AF_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'AF', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, AF_electrodes);
    FC_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'FC', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, FC_electrodes);
    FT_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'FT', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, FT_electrodes);
    CP_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'CP', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, CP_electrodes);
    TP_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'TP', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, TP_electrodes);
    PO_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'PO', IgnoreCase=false), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, PO_electrodes);
    F_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'F', IgnoreCase=true), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, F_electrodes);
    C_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'C', IgnoreCase=true), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, C_electrodes);
    T_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'T', IgnoreCase=true), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, T_electrodes);
    P_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'P', IgnoreCase=true), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, P_electrodes);
    O_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'O', IgnoreCase=true), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, O_electrodes);
    I_electrodes = sort(remaining_electrodes(cell2mat(arrayfun(@(str)startsWith(str, 'I', IgnoreCase=true), remaining_electrodes, 'UniformOutput', false))));
    remaining_electrodes = setdiff(remaining_electrodes, I_electrodes); % should be empty 


    frontal_electrodes = {'Fp1', 'Fp2', 'Fz', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', ...
                          'AF3', 'AF4', 'AF7', 'AF8', 'AFz'};
    temporal_electrodes = {'T7', 'T8', 'FT7', 'FT8', 'FT9', 'FT10', 'TP7', 'TP8', 'TP9', 'TP10'};
    parietal_electrodes = {'Pz', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'POz', 'PO3', ...
                           'PO4', 'PO7', 'PO8', 'CP1', 'CP2', 'CP3', 'CP4', 'CP5', 'CP6', 'CPz'};
    occipital_electrodes = {'O1', 'O2', 'Oz', 'Iz'};
    central_electrodes = {'FC5', 'FC1', 'C3', 'Cz', 'C4', 'FC6', 'FC2', 'FC3', 'C1', 'C5', 'C6', 'C2', 'FC4'};
    
    % Define center, left, and right electrode names
    center_electrodes = {'Fz', 'AFz', 'Cz', 'Pz', 'POz', 'Iz', 'CPz', 'Oz'};
    left_electrodes = {'Fp1', 'F3', 'F7', 'FT9', 'FC5', 'FC1', 'C3', 'T7', 'TP9', ...
                      'CP5', 'CP1', 'P3', 'P7', 'O1', 'F1', 'FT7', 'FC3', 'C1', ...
                      'C5', 'TP7', 'CP3', 'P1', 'P5', 'PO7', 'PO3', 'AF3', 'AF7', 'F5'};
    right_electrodes = {'Fp2', 'F4', 'F8', 'FT10', 'FC6', 'FC2', 'C4', 'T8', 'TP10', ...
                       'CP6', 'CP2', 'P4', 'P8', 'O2', 'F2', 'FT8', 'FC4', 'C2', ... 
                       'C6', 'TP8', 'CP4', 'P2', 'P6', 'PO8', 'PO4', 'AF4', 'AF8', 'F6'};
    
    regions = struct();
    
    % find intersect to get all values
    regions.left_temporal = intersect(left_electrodes, temporal_electrodes);
    
    regions.left_parietal = intersect(left_electrodes, parietal_electrodes);
    regions.center_parietal = intersect(center_electrodes, parietal_electrodes);
 
    regions.left_occipital = intersect(left_electrodes, occipital_electrodes);
    regions.center_occipital = intersect(center_electrodes, occipital_electrodes);
    regions.right_occipital = intersect(right_electrodes, occipital_electrodes);

    regions.right_parietal = intersect(right_electrodes, parietal_electrodes);

    regions.right_temporal = intersect(right_electrodes, temporal_electrodes);
    regions.center_temporal = intersect(center_electrodes, temporal_electrodes);

    regions.right_central = intersect(right_electrodes, central_electrodes);

    regions.right_frontal = intersect(right_electrodes, frontal_electrodes);
    regions.center_frontal = intersect(center_electrodes, frontal_electrodes);
    regions.left_frontal = intersect(left_electrodes, frontal_electrodes);
    
    regions.center_central = intersect(center_electrodes, central_electrodes);
    regions.left_central = intersect(left_electrodes, central_electrodes);    

    region_idx = struct();

    % find intersect to get all values
    [~,region_idx.left_temporal]=ismember(regions.left_temporal, all_electrodes);
   
    [~,region_idx.left_parietal]=ismember(regions.left_parietal, all_electrodes);
    [~,region_idx.center_parietal]=ismember(regions.center_parietal, all_electrodes);

    [~,region_idx.left_occipital]=ismember(regions.left_occipital, all_electrodes);
    [~,region_idx.center_occipital]=ismember(regions.center_occipital, all_electrodes);
    [~,region_idx.right_occipital]=ismember(regions.right_occipital, all_electrodes);

    [~,region_idx.right_parietal]=ismember(regions.right_parietal, all_electrodes);

    [~,region_idx.right_temporal]=ismember(regions.right_temporal, all_electrodes);
    [~,region_idx.center_temporal]=ismember(regions.center_temporal, all_electrodes);

    [~,region_idx.right_central]=ismember(regions.right_central, all_electrodes);

    [~,region_idx.right_frontal]=ismember(regions.right_frontal, all_electrodes);
    [~,region_idx.center_frontal]=ismember(regions.center_frontal, all_electrodes);
    [~,region_idx.left_frontal]=ismember(regions.left_frontal, all_electrodes);

    [~,region_idx.center_central]=ismember(regions.center_central, all_electrodes);
    [~,region_idx.left_central]=ismember(regions.left_central, all_electrodes);

    % Concatenate all fields into one vector
    region_fields = struct2cell(regions); % Convert struct fields to a cell array
    regions_concat = cat(2, region_fields{:}); % Concatenate the contents horizontally
    region_idx_fields = struct2cell(region_idx); % Convert struct fields to a cell array
    region_concat_idx = cat(2, region_idx_fields{:}); % Concatenate the contents horizontally
    
    % BY LETTER: (useful for connectivity matrix visualization)

    by_letter_labels = cat(1, Fp_electrodes, AF_electrodes, F_electrodes, FC_electrodes, FT_electrodes, C_electrodes, CP_electrodes, T_electrodes, TP_electrodes, P_electrodes, PO_electrodes, O_electrodes, I_electrodes);
    [~, by_letter_idx] = ismember(by_letter_labels, all_electrodes);  
    % [by_letter_labels, by_letter_idx] = sort(all_electrodes);

    % define output structure 
    electrode_organizations = struct();
    electrode_organizations.by_letter = struct();
    electrode_organizations.by_letter.idx = by_letter_idx.';
    electrode_organizations.by_letter.label = by_letter_labels.';
    electrode_organizations.by_region = struct();
    electrode_organizations.by_region.struct = struct();
    electrode_organizations.by_region.struct.idx = region_idx.';
    electrode_organizations.by_region.struct.label = regions.';
    electrode_organizations.by_region.idx = region_concat_idx.';
    electrode_organizations.by_region.label = regions_concat.';
end
