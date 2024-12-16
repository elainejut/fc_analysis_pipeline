function conn_mat_all_channels = insertMissingElectrodes(all_electrodes, good_electrodes, conn_mat)
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
    num_electrodes = length(all_electrodes);
    intersect = ismember(all_electrodes, good_electrodes);
    intersect_idx = find(intersect == 1);

    intermediate = NaN(num_electrodes, size(conn_mat, 2));
    intermediate(intersect_idx, :) = conn_mat;
    
    conn_mat_all_channels = NaN(size(intermediate, 1), num_electrodes);
    conn_mat_all_channels(:, intersect_idx) = intermediate;

    

    % 
    % missing_electrodes = ~ismember(all_electrodes, good_electrodes);
    % missing_electrodes_idx = find(missing_electrodes == 1);
    % % disp(missing_electrodes); % array with logical mask
    % 
    % % prev_idx = 1;
    % % for j = 1:length(missing_electrodes_idx)
    % %     desired_idx = missing_electrodes_idx(j);
    % %     B(prev_idx:2, prev_idx:3)
    % % 
    % %     prev_idx = desired_idx; 
    % % end
    % 
    % B = zeros(num_electrodes, size(A,2));
    % for i = 1:num_electrodes
    %     if missing_electrodes(i) == 1
    % 
    %     else
    % 
    %     end
    % end
    % 
    % 
    % B(A(:,1),:) = A;
    % B(B(:,1)==0,2:end) = NaN;
    % B(:,1) = 1:num_electrodes   % (Used 1:5 in my example below)
    % 
    % % conn_mat_all_channels = conn_mat(:);
    % % conn_mat_all_channels(end+1:9,end+1:9) = missing;
    % % 
    % % disp(conn_mat_all_channels);
    % 
end
