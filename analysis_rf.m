clear all 
close all 

%% open excel file and extract syngap and het data 
data = readtable('RF_data_area.xlsx');
syn = data(1:44,26);
syn = table2array(syn);

het = data(1:51,10);
het = table2array(het);


%% create statistical graphs 

figure(1)
errorbar(1,mean(syn),std(syn),'o-', 'Markersize', 10, 'Linewidth', 1)
hold on 
errorbar(2,mean(het),std(het),'o-', 'Markersize', 10, 'Linewidth', 1)
axis([0 3 -100 1000])
title('Syngap +/- Mutant vs Wild-Type Receptive Field Mean Area', 'Fontsize', 20)
legend({'Blue = Syngap+/- Mutant','Red = Wild Type'},'Location','northeast', 'FontSize',16)
xlabel('Syngap1 and WT', 'Fontsize', 12)
ylabel('Mean Area of the Receptive Field', 'Fontsize', 12)


figure(2)
histogram(syn,75)
hold on
histogram(het,75)
axis([0 2500 0 8.5])
title('Syngap +/- Mutant vs Wild-Type Receptive Field Area Data', 'Fontsize', 20)
legend({'Blue = Syngap+/- Mutant','Red = Wild Type'},'Location','northeast', 'FontSize',16)
xlabel('Area of the Receptive Field', 'Fontsize', 12)
ylabel('Number of Neurons', 'Fontsize', 12)


% calculate p value 
[h,p] = ttest2(syn,het)