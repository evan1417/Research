% extract receptive field data 
clear all 
close all 
% Get experiment parameters
% get stim times
open_mat = matfile('__name___CombinedVisualData.mat');
stimon = open_mat.save_data(1,:);
difference_stim = diff(stimon);
stimtime = find(difference_stim>0.5)+1;
clear stimon digin dstim % get rid of large dataset

% open files necessary for binning 
clustergroup = readtable('cluster_group.csv');
spiketime = readNPY('spike_times.npy');
spikecluster = readNPY('spike_clusters.npy');


% prune data to print good clusters ids
label = clustergroup(:,2);
cellVal = table2cell(label);
good_label = strcmp(cellVal,'good');
cluster_id = table2array(clustergroup(:,1));
good_id = cluster_id(good_label == 1);
number_ids = length(good_id);
good_id_cell = cell(number_ids,1);

% print the unique cluster numbers for spikes 
clusters = unique(spikecluster);
commonValues = intersect(good_id, clusters);
        

%cluster_group = clusterSort(spikecluster, spiketime);
ncells = length(commonValues);
spikes = cell(ncells,1);

% get trial types 
cam_mouse = readtable('cam_name.csv');
trial_types = table2array(unique(cam_mouse(:,2))); 
n_trial_types = length(trial_types);

% Just look at receptive field 
all_trial_types = table2array(cam_mouse(:,2));
all_rf_y = table2array(cam_mouse(:,7));
all_rf_x = table2array(cam_mouse(:,8));
nall_trial_types = length(all_trial_types);
receptive_field = strcmp(all_trial_types,'rf');
nrf = sum(receptive_field);
receptive_y = all_rf_y(receptive_field==1);
receptive_x = all_rf_x(receptive_field==1);
stimrf = stimtime(receptive_field==1);
ntrials = length(stimrf);
df = diff(stimrf);
binstimrf = stimrf;


binstimrf(1, 9861) = stimrf(1, 9860)+20000; % this should be changed with the updated python code


% Calculate out the last stim time plus the additional time
last_stim_time = stimrf(end) + 20000;


% group up the x,y coordinates with their respective stimtime 
col = stimrf.';
ncol = length(col);
z = [receptive_y receptive_x col];
lenz = length(z);


% combine good spike clusters with its spiketime
for k = 1:ncells
    spikes{k} = spiketime(spikecluster==commonValues(k));
end



z2 = z(:,1:2);

% change bad boys depending on the data set 
% input the y value first, then x for any abnormal data point
bad_boys = [140 100]; 
hh = ismember(z2, bad_boys, 'rows'); 
hh = double(hh); 


for k = 1:ncells 

    samples = spikes{k};
    [N,edges] = histcounts(samples,'BinEdges',binstimrf);

    for j = 1:ntrials
        M = zeros(232,136);
 
        if hh(j,1) == 1
            
        else
            % first value is the y and second value is x
            M(z(j, 1),z(j,2)) = 1;

            M = conv2(squeeze(M(:, :)),ones(8),'same');
        end        
   
        test = N(j) * M;
        if j == 1
            total_val = test;
        else
            total_val = total_val + test;
        end



    end
    
    graph1 = total_val * (1/20);
    
    
    figure(1)
    
    % set subplot to number of spikes 
    subplot(5,5,k);
    imagesc(graph1)
    colorbar
    
end 
