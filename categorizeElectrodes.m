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

    % BY REGION: 

    % Define electrode regions
    frontal_electrodes = {'Fp1', 'Fp2', 'Fz', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', ...
                          'AF3', 'AF4', 'AF7', 'AF8', 'AFz'};
    temporal_electrodes = {'T7', 'T8', 'FT7', 'FT8', 'FT9', 'FT10', 'TP7', 'TP8', 'TP9', 'TP10'};
    parietal_electrodes = {'Pz', 'P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7', 'P8', 'POz', 'PO3', ...
                           'PO4', 'PO7', 'PO8', 'CP1', 'CP2', 'CP3', 'CP4', 'CP5', 'CP6', 'CPz'};
    occipital_electrodes = {'O1', 'O2', 'Oz', 'Iz'};
    central_electrodes = {'FC5', 'FC1', 'C3', 'Cz', 'C4', 'FC6', 'FC2', 'FC3', 'C1', 'C5', 'C6', 'C2', 'FC4', 'FCz'};
    
    % Define center, left, and right electrode names
    center_electrodes = {'Fz', 'AFz', 'Cz', 'Pz', 'POz', 'Iz', 'CPz', 'FCz', 'Oz'};
    left_electrodes = {'Fp1', 'F3', 'F7', 'FT9', 'FC5', 'FC1', 'C3', 'T7', 'TP9', ...
                      'CP5', 'CP1', 'P3', 'P7', 'O1', 'F1', 'FT7', 'FC3', 'C1', ...
                      'C5', 'TP7', 'CP3', 'P1', 'P5', 'PO7', 'PO3', 'AF3', 'AF7', 'F5'};
    right_electrodes = {'Fp2', 'F4', 'F8', 'FT10', 'FC6', 'FC2', 'C4', 'T8', 'TP10', ...
                       'CP6', 'CP2', 'P4', 'P8', 'O2', 'F2', 'FT8', 'FC4', 'C2', ... 
                       'C6', 'TP8', 'CP4', 'P2', 'P6', 'PO8', 'PO4', 'AF4', 'AF8', 'F6'};
    
    regions = struct();
    
    % find intersect to get all values
    regions.left_frontal = intersect(left_electrodes, frontal_electrodes);
    regions.center_frontal = intersect(center_electrodes, frontal_electrodes);
    regions.right_frontal = intersect(right_electrodes, frontal_electrodes);
    
    regions.left_central = intersect(left_electrodes, central_electrodes);
    regions.center_central = intersect(center_electrodes, central_electrodes);
    regions.right_central = intersect(right_electrodes, central_electrodes);
    
    regions.left_temporal = intersect(left_electrodes, temporal_electrodes);
    regions.center_temporal = intersect(center_electrodes, temporal_electrodes);
    regions.right_temporal = intersect(right_electrodes, temporal_electrodes);
    
    regions.left_parietal = intersect(left_electrodes, parietal_electrodes);
    regions.center_parietal = intersect(center_electrodes, parietal_electrodes);
    regions.right_parietal = intersect(right_electrodes, parietal_electrodes);
    
    regions.left_occipital = intersect(left_electrodes, occipital_electrodes);
    regions.center_occipital = intersect(center_electrodes, occipital_electrodes);
    regions.right_occipital = intersect(right_electrodes, occipital_electrodes);

    region_idx = struct();

    % find intersect to get all values
    [~,region_idx.left_frontal]=ismember(regions.left_frontal, all_electrodes);
    % region_idx.left_frontal = get_indices(regions.left_frontal, electrodes);
    % region_idx.center_frontal = get_indices(regions.center_frontal, electrodes);
    % region_idx.right_frontal = get_indices(regions.right_frontal, electrodes);
    % 
    % region_idx.left_central = get_indices(regions.left_central, electrodes);
    % region_idx.center_central = get_indices(regions.center_central, electrodes);
    % region_idx.right_central = get_indices(regions.right_central, electrodes);
    % 
    % region_idx.left_temporal = get_indices(regions.left_temporal, electrodes);
    % region_idx.center_temporal = get_indices(regions.center_temporal, electrodes);
    % region_idx.right_temporal = get_indices(regions.right_temporal, electrodes);
    % 
    % region_idx.left_parietal = get_indices(regions.left_parietal, electrodes);
    % region_idx.center_parietal = get_indices(regions.center_parietal, electrodes);
    % region_idx.right_parietal = get_indices(regions.right_parietal, electrodes);
    % 
    % region_idx.left_occipital = get_indices(regions.left_occipital, electrodes);
    % region_idx.center_occipital = get_indices(regions.center_occipital, electrodes);
    % region_idx.right_occipital = get_indices(regions.right_occipital, electrodes);

    
    % ALPHABETICAL:

    [alphabetical_labels, alphabetical_idx] = sort(all_electrodes);

    % define output structure 
    electrode_organizations = struct();
    electrode_organizations.alphabetical = struct();
    electrode_organizations.alphabetical.idx = alphabetical_idx.';
    electrode_organizations.alphabetical.label = alphabetical_labels.';
    electrode_organizations.by_region = struct();
    electrode_organizations.by_region.struct = struct();
    electrode_organizations.by_region.struct.idx = region_idx.';
    electrode_organizations.by_region.struct.label = regions.';
    % electrode_organizations.by_region.idx = region_concat_idx.';
    % electrode_organizations.by_region.label = regions_concat.';

end
