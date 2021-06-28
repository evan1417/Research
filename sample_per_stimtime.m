% extract receptive field data 

% Get experiment parameters
% get stim times
open_mat = matfile('SynGAP1_Mouse1_CombinedVisualData.mat');
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
cam_mouse = readtable('cam-SynGAP1-Mouse1-PerStimData-2021_Jun_22_1230.csv');
trial_types = table2array(unique(cam_mouse(:,2))); 
n_trial_types = length(trial_types);

% Just look at receptive field 
all_trial_types = table2array(cam_mouse(:,2));
all_rf_x = table2array(cam_mouse(:,7));
all_rf_y = table2array(cam_mouse(:,8));
nall_trial_types = length(all_trial_types);
receptive_field = strcmp(all_trial_types,'rf');
nrf = sum(receptive_field);
receptive_x = all_rf_x(receptive_field==1);
receptive_y = all_rf_y(receptive_field==1);
stimrf = stimtime(receptive_field==1);
ntrials = length(stimrf);
df = diff(stimrf);
stimrf(1,9861) = 84654830;


% group up the x,y coordinates with their respective stimtime 
col = stimrf.';
ncol = length(col);
%z = [receptive_x receptive_y col];
lenz = length(z);


% combine good spike clusters with its spiketime
for k = 1:ncells
    spikes{k} = spiketime(spikecluster==commonValues(k));
end

% create bins of receptive field stimulus times 
% figure(1)
for j = 1:ncells
    for i = spikes{j,1}   
        samples = i(:,1);
        [N,edges] = histcounts(samples,'BinEdges',stimrf);
        subplot(4,6,j);
        imagesc(N)
        colorbar
    end 
end

% create matrices for every x,y coordinate and plot them with a surrounding
% group of 1's
% for l = 1:493
%     M(l) = mat2cell(zeros(232,136),232,136);
%     M{l}(z(l,1),z(l,2)) = 1;
%     M{l} = conv2(M{l},ones(8),'same');
% end



% extract the number of stimuli in each bin per neuron 
