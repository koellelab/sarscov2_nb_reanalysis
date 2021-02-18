function void = mainPlotNbresults(void)

clear all; close all; clc;

load('Popa_data_reanalysis_1percent_NbAnalysis')
data_onePercent = data; data = [];

load('Popa_data_reanalysis_3percent_NbAnalysis')
data_threePercent = data; data = [];

load('Popa_data_reanalysis_6percent_NbAnalysis')
data_sixPercent = data; data = [];


for i = 1:39
    if data_onePercent(i).n_variants_at_cutoff == 0
        data_onePercent(i).Nb_low = NaN;
        data_onePercent(i).Nb_high = NaN;
    end
    if data_threePercent(i).n_variants_at_cutoff == 0
        data_threePercent(i).Nb_low = NaN;
        data_threePercent(i).Nb_high = NaN;
    end
    if data_sixPercent(i).n_variants_at_cutoff == 0
        data_sixPercent(i).Nb_low = NaN;
        data_sixPercent(i).Nb_high = NaN;
    end
   
   semilogy(i, data_onePercent(i).Nb_MLE, 'ko'); hold on;
   
   semilogy([i i], [data_onePercent(i).Nb_low data_onePercent(i).Nb_high], 'k');
   semilogy(i, data_threePercent(i).Nb_MLE, 'bo'); hold on;
   semilogy([i i], [data_threePercent(i).Nb_low data_threePercent(i).Nb_high], 'b');
   semilogy(i, data_sixPercent(i).Nb_MLE, 'ro'); hold on;
   semilogy([i i], [data_sixPercent(i).Nb_low data_sixPercent(i).Nb_high], 'r');
end
axis([0 40 0.5 11000])

ylabel('Bottleneck size')
xticks([1:39])

%for i = 1:39
%    ts1 = strcat(int2str(data_onePercent(i).donor), '.', int2str(data_onePercent(i).recipient));
%end
%xticklabels({ts(1), ts(2), ts(3), ts(4), ts(5), ts(6), ts(7), ts(8), ts(9), ts(10), ts(11), ts(12), ts(13), ts(14), ts(15), ts(16), ts(17), ts(18), ts(19), ts(20), ts(21), ts(22), ts(23), ts(24), ts(25), ts(26), ts(27), ts(28),ts(29), ts(30), ts(31), ts(32), ts(33), ts(34), ts(35), ts(36), ts(37), ts(38), ts(39)});

xticklabels({'003.143', '003.146', '146.147', '146.172', '146.157', '146.159', '146.150', '146 1056', '159.168', '159.169', '150.162',  '162.161', '1056.166', '1056.176', '1056.171', '1056.1057', '1056.1067', '171.180', '180.197', '1067.183', '1057.175', '1057.177', '1057.1058', '1057.187', '1058.1059', '1058.1064', '1058.1063', '1059.1065', '1059.1062', '1062.218', '1062.219', '1062.217', '217.256', '187.1068', '194.193', '194.195', '194.196', '198.230', '273.271'})
xtickangle(90)
xlabel('Transmission pair')

return;

figure;
for i = 1:39
    semilogy(max(data_onePercent(i).donor_iSNVs), data_onePercent(i).Nb_MLE, 'ko'); hold on;
    semilogy([max(data_onePercent(i).donor_iSNVs) max(data_onePercent(i).donor_iSNVs)], [data_onePercent(i).Nb_low data_onePercent(i).Nb_high], 'k'); 
    semilogy(max(data_threePercent(i).donor_iSNVs), data_threePercent(i).Nb_MLE, 'bo'); hold on;
    semilogy([max(data_threePercent(i).donor_iSNVs) max(data_threePercent(i).donor_iSNVs)], [data_threePercent(i).Nb_low data_threePercent(i).Nb_high], 'b'); 
    semilogy(max(data_sixPercent(i).donor_iSNVs), data_sixPercent(i).Nb_MLE, 'ro'); hold on;
    semilogy([max(data_sixPercent(i).donor_iSNVs) max(data_sixPercent(i).donor_iSNVs)], [data_sixPercent(i).Nb_low data_sixPercent(i).Nb_high], 'r'); 
end
ylabel('Bottleneck size')
xlabel('Max donor iSNV frequency')
axis([0 0.5 0.5 11000])


if 1 % to make spreadsheet:
for i = 1:39
   donor_list(i,1) = data_onePercent(i).donor;
   recipient_list(i,1) = data_onePercent(i).recipient;
   onepercent_MLE(i,1) = data_onePercent(i).Nb_MLE;
   onepercent_lower(i,1) = data_onePercent(i).Nb_low;
   onepercent_upper(i,1) = data_onePercent(i).Nb_high;
   
   threepercent_MLE(i,1) = data_threePercent(i).Nb_MLE;
   threepercent_lower(i,1) = data_threePercent(i).Nb_low;
   threepercent_upper(i,1) = data_threePercent(i).Nb_high;
   
   sixpercent_MLE(i,1) = data_sixPercent(i).Nb_MLE;
   sixpercent_lower(i,1) = data_sixPercent(i).Nb_low;
   sixpercent_upper(i,1) = data_sixPercent(i).Nb_high;
    
end
donor_list
pause
recipient_list
pause
onepercent_MLE
pause
onepercent_lower
pause
onepercent_upper
pause
threepercent_MLE
pause
threepercent_lower
pause
threepercent_upper
pause
sixpercent_MLE
pause
sixpercent_lower
pause
sixpercent_upper
pause
   
end
