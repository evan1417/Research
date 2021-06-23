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

% create bins of spikeclusters with their times and plot as a figure
figure
for j = 1:ncells
    for i = spikes{j,1}   
        samples = i(:,1); 
        [N,edges] = histcounts(samples,'BinWidth',2000)
        subplot(4,8,j); 
        bar(N);
    end 
end











