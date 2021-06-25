open_mat = matfile('SynGAP1_Mouse1_CombinedVisualData.mat');
stimon = open_mat.save_data(1,:);
dstim = diff(stimon);
stimtime = find(dstim>0.5)+1;
clear stimon digin dstim % get rid of large dataset

cam = readtable('cam-SynGAP1-Mouse1-PerStimData-2021_Jun_22_1230.csv');

%parsing cam mouse for data related to receptive field
all_trial_types = table2array(cam(:,2));
nall_trial_types = length(all_trial_types);
allrfx = table2array(cam(:,7)); 
allrfy = table2array(cam(:,8));
receptive_field = strcmp(all_trial_types,'rf');
nrf = sum(receptive_field);

stimrf = stimtime(receptive_field==1);
ntrials = length(stimrf);

% x and ys that correspond to only the rf label
receptive_x = allrfx(receptive_field==1);
receptive_y = allrfy(receptive_field==1);

% group up the x,y coordinates with their respective stimtime 
col = stimrf.';
ncol = length(col);
z = [receptive_x receptive_y col];
lenz = length(z);

% create matrices for every x,y coordinate and plot them with a surrounding
% group of 1's
for l = 1:493
    M(l) = mat2cell(zeros(232,136),232,136);
    M{l}(z(l,1),z(l,2)) = 1;
    M{l} = conv2(M{l},ones(8),'same');
end


