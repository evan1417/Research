% Final Tuning Curve Analysis Script
clear all
close all

%% Load in All Data from the holding folders
% Select Master Data Folder
% data_folder_dir = append(uigetdir(matlabroot,'Select main folder'), '\');
data_folder_dir = "C:\Users\ande1227\Desktop\mouse_analysis\";

% Loop through the data files present in the raw intan data folder and
% extract data
folder_info = dir(data_folder_dir);
for file_num = 3:5 % starts at 3 to miss the . and .. directories
    subfolder_path = data_folder_dir + folder_info(file_num).name + "\";
%     [data_folder_dir, folder_info(file_num).name, "\"];
    subfolder_info = dir(subfolder_path);
    
    % Will need to make a more robust check mechanism for loading in the
    % data files.  Currently hoping they are alwasys in the same order
    % Stim Time .mat file loading
    visual_stim_file = strcat(subfolder_path, subfolder_info(3).name);
    visual_stim_times = load(visual_stim_file);
    visual_stim_times = visual_stim_times.save_data;
    
    % Stim Data .csv file loading
    visual_ident_file = strcat(subfolder_path, subfolder_info(4).name);
    visual_ident_data = readtable(visual_ident_file);
    
    % Other Neuron Data
    cluster_grouping = readtable(strcat(subfolder_path, subfolder_info(5).name));
    spike_ids = readNPY(strcat(subfolder_path, subfolder_info(6).name));
    spike_times = readNPY(strcat(subfolder_path, subfolder_info(7).name));


    %% Begin Isolating Required Data
    % Isolate the stimulation times from the visual stimulation timing data.
    stim_time_diff = diff(visual_stim_times);
    stim_timing = find(stim_time_diff > 0.5) + 1;
    clear visual_stim_times stim_time_diff

    % Isolate the types of trials of trials.
    trial_types = table2array(unique(visual_ident_data(:, 2)));

    % Isolate the drifting grating data identifiers.
    trial_ident = table2array(visual_ident_data(:, 2));
    trial_dc = table2array(visual_ident_data(:, 10));
    trial_dd = table2array(visual_ident_data(:, 14));

    % Slice the data to only cover the drifting grating portion.
    dotfield_trials = strcmp(trial_ident, 'df');
    dotfield_dc = trial_dc(dotfield_trials == 1);
    dotfield_dd = trial_dd(dotfield_trials == 1);
    dotfield_stim_timing = stim_timing(dotfield_trials == 1);   
    num_trials = length(dotfield_stim_timing);

    % Remove extraneous data to save RAM.
    clear visual_ident_data trial_ident trial_dc trial_dd

    % Compute response vectors.  (Likely is the number of unique orientations
    % and spatial frequencies)
    unique_dd = unique(dotfield_dd);
    unique_dc = unique(dotfield_dc);
  
    num_dd = length(unique_dd);
    num_dc = length(unique_dc);

    %% Isolate Neurons

    % Determine the quality of each cluster
    good_clusters = strcmp(table2cell(cluster_grouping(:, 2)), 'good');
    good_cluster_ids = table2array(cluster_grouping(good_clusters == 1, 1));
    num_good_clusters = length(good_cluster_ids);

    % Identify the unique clusters present in the data
    clusters = unique(spike_ids);

    % Isolate the spike times for each cluster and the number of spikes per
    % cluster
    spikes = cell(length(num_good_clusters), 1);
    spike_count = NaN(length(num_good_clusters), 1);

    for k = 1:num_good_clusters
        spikes{k} = spike_times(spike_ids == good_cluster_ids(k));
        spike_count(k) = length(spikes{k});
    end

    % Remove neurons with low spiking rates.
    spike_threshold = 1;

    kept_indexes = spike_count > spike_threshold;
    num_cells = length(kept_indexes);
    kept_cells = spikes(kept_indexes == 1);

    %% Compute Tuning Curves

    % Calculate the tuning responses by summing the number of spikes in
    % response to a stimulus.
    tuning_responses = NaN(num_trials, num_cells);
    for k = 1:num_cells
        for g = 1:num_trials
            % This line needs to be updated to account for the exact end of a
            % stimulus, or a specified time after the exact end of a
            % stimulus. Would a buffer zone be warrented?
            tuning_responses(g, k) = sum(kept_cells{k} > dotfield_stim_timing(g) & kept_cells{k} < dotfield_stim_timing(g) + 0.5 * 20000);
            first_stim = kept_cells{k}(kept_cells{k} > dotfield_stim_timing(g) & kept_cells{k} <= dotfield_stim_timing(g) + 10000);
            if ~isempty(first_stim)
            
                stim_diff = first_stim(1) - dotfield_stim_timing(g);
                rast_save(g, k) = stim_diff;
            
            else
                rast_save(g, k) = NaN;
            end
            
%             stim_diff = first_stim - grating_stim_timing(g);
%             stim_save{sf, 
            
%             kept_cells{cell_num}(kept_cells{cell_num} > onset(q) & kept_cells{cell_num} <= onset(q) + 20000)
%                 stim_diff(q) = first_stim - onset(q);
%                 stim_save(k, q, cell_num) = stim_diff(q);
        end
    end
    
    
    % Save the raster data
%     final_rast_save{file_num

    % Calculate the actual tuning curves by orientation and spatial frequency.
    tuning_curve = NaN(num_dd, num_dc, num_cells);
    tuning_curve_z = NaN(num_dd, num_dc, num_cells);
    tuning_curve_std = NaN(num_dd, num_dc, num_cells);

%     cell_array = cell(num_sf, 1);

    for g = 1:num_dd
        test_data = NaN(num_dd, 20, num_cells);
    %     for k = 1:num_cells
    %         temp_array = NaN(num_ori, num_trials, num_cells);
        for j = 1:num_dc % switch from 1 or 5 to see the dot coherance thresh
            hello1 = tuning_responses(dotfield_dd == unique_dd(g) & dotfield_dc == unique_dc(j), :); % x: trial, y: cell
            hello2 = rast_save(dotfield_dd == unique_dd(g) & dotfield_dc == unique_dc(j), :);
            tuning_curve(g, j, :) = mean(tuning_responses(dotfield_dd == unique_dd(g) & dotfield_dc == unique_dc(j), :));
            for k = 1:num_cells
                test_data(j, :, k) = hello1(:, k)';
                rast_data(j, :, k) = hello2(:, k)';
            end
        end
        cell_array{g, file_num - 2} = test_data;
        final_rast_save{g, file_num - 2} = rast_data;
    end
    num_sqrt = ceil(sqrt(num_good_clusters));
    figure
    for g = 1:num_dc % switch from 1 or 5 to see the dot coherance thresh
        for k = 1:num_good_clusters
            subplot(num_sqrt, num_sqrt, k)
            hold on
            plot(tuning_curve(:, g, k))
            graph_title = int2str(k);
            title(k)
        end
    end
end

% Remove extraneous data
clear hello1 test_data

%% Population Level OSI Per Spatial Frequency
% The orientation data
orientations = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330];

% Seperate out the OSI values
for dc = 1:length(unique_dc)
    max_values = [];
    cell_num = 1;
    for mouse_num = 1:size(cell_array, 2)
        current_value = cell_array{dc, mouse_num};
        mean_value = mean(current_value, 2);
        mean_value = squeeze(mean_value);
        
        for neuron_num = 1:size(mean_value, 2)
            [max_values(cell_num, 1), index] = max(mean_value(:, neuron_num));
            
        
        
%             [max_values(cell_num, 1), index] = max(mean_value);
            max_values(cell_num, 2) = orientations(index);
            
            %Double check this logic, it is likely very wrong
            if index + 2 >= length(orientations)
                orth_ind = index + 4 - length(orientations);
            else
                orth_ind = index + 3;
            end
            
            max_values(cell_num, 3) = mean_value(orth_ind, neuron_num);
            max_values(cell_num, 4) = orth_ind;
            cell_num = cell_num + 1;
        end
    end
    orientation_array{dc} = max_values;
end

% Create the master syngap vs non-syngap list
syngap_mice = [1];
cell_num = 1;
mouse_ident = [];
for mouse_num = 1:size(cell_array, 2)
    num_neurons = size(cell_array{1, mouse_num}, 3);
    
    if find(syngap_mice == mouse_num)
        mouse_type = 1;
        mouse_add = ones(num_neurons, 1);
    else
        mouse_type = 0;
        mouse_add = zeros(num_neurons, 1);
    end
    
    
    mouse_ident = [mouse_ident; mouse_add];
    cell_num = cell_num + 1;
end

% Calculate the Occular selectivity index for each sf
for dc = 1:length(unique_dc)
    pulled_data = orientation_array{dc};
    OSI(:, dc) = (pulled_data(:, 1) - pulled_data(:, 3)) ./ (pulled_data(:, 1) + pulled_data(:, 3));
end

% OSI Thresholding Calculations (over threshold is considered tuned
% neurons)
OSI_threshold = 0.3;
over_threshold = OSI > 0.3;

%% Circular Variance Analysis
% Need to calculate tuning curve metrics
for dc = 1:length(unique_dc)
    cell_num = 1;
    total_mean_value = [];
    for mouse_num = 1:size(cell_array, 2)
        current_value = cell_array{dc, mouse_num};
        mean_value = mean(current_value, 2);
        mean_value = squeeze(mean_value);
        
        total_mean_value = [total_mean_value, mean_value];
%         for neuron = 1:size(mean_value, 2)
%             mean_values
%         mean_values(, sf) = 
    end
    dc_seperated_mean_values{dc} = total_mean_value;
end

% Now we need to calculate the data out on a per data basis
for dc = 1:length(unique_dc)
    value_array = dc_seperated_mean_values{dc};
    for cell_num = 1:size(value_array, 2)
        temp_array = value_array(:, cell_num);
        convt_orientations = orientations * (pi/180);
        sin_array = temp_array .* sin(orientations');
        cos_array = temp_array .* cos(orientations');
        A = sum(sin_array);
        B = sum(cos_array);
        
        theta = atan(B/A);
        R = (A^2 + B^2) ^ (1/2);
        
        per_dc_array(cell_num, 1) = R;
        per_dc_array(cell_num, 2) = theta;
    end
    circular_variance{dc} = per_dc_array;
end

%% Raster Plot and First Firing Analysis
for dc = 1:length(unique_dc)
    cell_number = 1;
    for mouse_num = 1:size(final_rast_save, 2)
        current_value = final_rast_save{dc, mouse_num};
        for neuron_num = 1:size(current_value, 3)
            if cell_number == 1
                overall_data = current_value(:, :, neuron_num);
                if mouse_ident(cell_number) == 1
                    rast_ident = ones(20, 1);
                else
                    rast_ident = zeros(20, 1);
                end
            else
                overall_data = [overall_data current_value(:, :, neuron_num)];
                if mouse_ident(cell_number) == 1
                    rast_ident = [rast_ident; ones(20, 1)];
                else
                    rast_ident = [rast_ident; zeros(20, 1)];
                end
            end
            cell_number = cell_number + 1;
        end
    end
end

% Split Raster Data
% find_data = find(rast_ident);
syngap_data = overall_data(:, find(rast_ident));
wt_data = overall_data(:, find(~rast_ident));

clear rast_ident overall_data rast_data rast_save final_rast_save

% Mean values for each orientation
mean_syngap = mean(nonzeros(syngap_data), 'all');
mean_wt = mean(nonzeros(wt_data), 'all');

% STD values for each orientation
std_syngap = std(double(nonzeros(syngap_data)), 0, 'all');
std_wt = std(double(nonzeros(wt_data)), 0, 'all');

% Graph the Firing Time Data
figure
hold on
errorbar(mean_syngap, std_syngap, 'o')
errorbar(mean_wt, std_wt, 'o')
title("Average Stimulus Response Time (Samples)")
xlabel("Image Number")
ylabel("# of Samples Past Stimulus Onset")
legend("Syngap", "WT")

% Mean converted data
mean_syngap_convt = mean(nonzeros(syngap_data) * 0.05, 'all');
mean_wt_convt = mean(nonzeros(wt_data) * 0.05, 'all');

% Std converted data
std_syngap_convt = std(double(nonzeros(syngap_data) * 0.05), 0, 'all');
std_wt_convt = std(double(nonzeros(wt_data) * 0.05), 0, 'all');

figure
hold on
errorbar(mean_syngap_convt, std_syngap_convt, 'o')
errorbar(mean_wt_convt, std_wt_convt, 'o')
title("Average Stimulus Response Time (ms)")
xlabel("Image Number")
ylabel("ms Past Stimulus Onset")
legend("Syngap", "WT")


%% Graphing of Tuning Curve Data
% OSI
figure
hold on
columns = 2;
rows = length(unique_dc);
plot_num = 1;
for dc = 1:length(unique_dc)
    subplot(rows, columns, plot_num)
    histogram(OSI(find(mouse_ident), dc), 10)
    subplot(rows, columns, plot_num + 1)
    histogram(OSI(find(~mouse_ident), dc), 10)
    plot_num = plot_num + 2
end

% Circular Variance
figure
hold on
plot_num = 1;
for dc = 1:length(unique_dc)
    subplot(rows, columns, plot_num)
    histogram(circular_variance{dc}(find(mouse_ident), 1), 10)
    subplot(rows, columns, plot_num + 1)
    histogram(circular_variance{dc}(find(~mouse_ident), 1), 10)
    plot_num = plot_num + 2;
end
    
    

% subplot(rows, columns)