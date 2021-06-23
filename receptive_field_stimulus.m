% extract receptive field data 

% Get experiment parameters
% get stim times
open_mat = matfile('SynGAP1_Mouse1_CombinedVisualData.mat');
stimon = open_mat.save_data(1,:);
difference_stim = diff(stimon);
stimtime = find(difference_stim>0.5)+1;
clear stimon digin dstim % get rid of large dataset

% get trial types 
cam_mouse = readtable('cam-SynGAP1-Mouse1-PerStimData-2021_Jun_22_1230.csv');
trial_types = table2array(unique(cam_mouse(:,2))); 
n_trial_types = length(trial_types);

% Just look at receptive field 
all_trial_types = table2array(cam_mouse(:,2));
all_sf = table2array(cam_mouse(:,3));
all_ori = table2array(cam_mouse(:,4));
nall_trial_types = length(all_trial_types);
receptive_field = strcmp(all_trial_types,'rf');
nrf = sum(receptive_field);
receptive_sf = all_sf(receptive_field==1);
receptive_ori = all_ori(receptive_field==1);
stimrf = stimtime(receptive_field==1);
ntrials = length(stimrf);

% create bins of receptive field stimulus times 
figure(99)
for j = 1:ncells
    for i = spikes{j,1}   
        samples = i(:,1); 
        [N,edges] = histcounts(samples,'BinWidth',ntrials);
        subplot(4,6,j); 
        bar(N);
        title(j)
    end 
end


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

% combine good spike clusters with its spiketime
for k = 1:ncells
    spikes{k} = spiketime(spikecluster==commonValues(k));
end
