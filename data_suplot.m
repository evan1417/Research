% extracts data points to put into a subplot for my poster 

%% set variables 
clear all 
close all 

%% load in data 
bad_syngap = load('bad_syngap.mat');
bad_syngap = cell2mat(struct2cell(bad_syngap));

good_syngap = load('good_syngap.mat');
good_syngap = cell2mat(struct2cell(good_syngap));

good_wt = load('good_wt.mat');
good_wt = cell2mat(struct2cell(good_wt));

bad_wt = load('bad_wt.mat');
bad_wt = cell2mat(struct2cell(bad_wt));


%% create the subplot of all four graphs 
sgtitle('Comparison of Syngap +/- and Wild Type Receptive Field Sizes', 'Fontsize',13)

subplot(2,2,1)
imagesc(good_syngap);
colorbar
title('Good Syngap +/- Receptive Field')

subplot(2,2,2)
imagesc(bad_syngap);
colorbar
title('Bad Syngap +/- Receptive Field')

subplot(2,2,3)
imagesc(good_wt);
colorbar
title('Good WT Receptive Field')

subplot(2,2,4)
imagesc(bad_wt);
colorbar
title('Bad WT Receptive Field')
