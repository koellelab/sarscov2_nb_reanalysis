function void = main_estimateOverallNb(void)

clear all; close all; clc;

infile = 'Popa_data_reanalysis_6percent_NbAnalysis';
outfile = 'Popa_overall_Nb_estimate_6percent';

load(infile);

n_TPs = length(data);

min_Nb = 1; max_Nb = 20;
overall_Nb = min_Nb:max_Nb;
    
for i = 1:n_TPs

    if data(i).n_variants_at_cutoff > 0
        data(i).logL_for_overall = GetLogL_forNb(data(i), var_calling_threshold, overall_Nb);
    else
        data(i).logL_for_overall = NaN*ones(1, length(min_Nb:max_Nb));
    end
    plot(overall_Nb, data(i).logL_for_overall); hold on;
    save(outfile, 'data', 'var_calling_threshold');
end

% zero-truncated Poisson distn
cntr = 1;
%lambda_list = 0.1:0.1:10;
lambda_list = 0.1:0.001:2.5;
for this_lambda = lambda_list
    pmf_poiss = poisspdf(overall_Nb, this_lambda)/(sum(poisspdf(overall_Nb, this_lambda)));
    for i = 1:n_TPs
        logL_allTPS(i) = GetLogLForMeanNb(pmf_poiss, overall_Nb, data(i).logL_for_overall);
        if isnan(logL_allTPS(i))
            logL_allTPS(i) = 0;
        end 
    end
    
    overall_Nb_logL_allTPS(cntr) = sum(logL_allTPS);
    cntr = cntr + 1;
end
mean_Nb = (lambda_list)./(1-exp(-lambda_list));
figure(1); subplot(1,2,1); plot(mean_Nb, overall_Nb_logL_allTPS);
loc_max = find(overall_Nb_logL_allTPS == max(overall_Nb_logL_allTPS));
est_mean_Nb_allTPS = mean_Nb(loc_max)
est_lambda_allTPS = lambda_list(loc_max);

est_pmf_poiss_allTPS = poisspdf(overall_Nb, est_lambda_allTPS)/(sum(poisspdf(overall_Nb, est_lambda_allTPS)));
figure(2); subplot(1,2,1); 
bar(overall_Nb, est_pmf_poiss_allTPS)
xlabel('transmission of N_b virions'); ylabel('probability');
title('all TPs with donor iSNVs > 6%');
save(outfile, 'data', 'var_calling_threshold', 'overall_Nb', 'est_pmf_poiss_allTPS');
overall_Nb'
est_pmf_poiss_allTPS'
pause

% zero-truncated Poisson distn - only low CT
cntr = 1;
%lambda_list = 0.1:0.1:10;
lambda_list = 0.1:0.001:2.5;
for this_lambda = lambda_list
    pmf_poiss = poisspdf(overall_Nb, this_lambda)/(sum(poisspdf(overall_Nb, this_lambda)));
    for i = 1:n_TPs
        logL_lowCT_TPS(i) = GetLogLForMeanNb(pmf_poiss, overall_Nb, data(i).logL_for_overall);
        if isnan(logL_lowCT_TPS(i))
            logL_lowCT_TPS(i) = 0;
        end 
        if ~isnan(logL_lowCT_TPS(i))
            index_donor = find(CT_data.sample_name == data(i).donor);
            this_CT = CT_data.CT_value(index_donor);
            if this_CT > 30
                logL_lowCT_TPS(i) = 0;
            end
        end
    end
    %logL_lowCT_TPS
    overall_Nb_logL_lowCT_TPS(cntr) = sum(logL_lowCT_TPS);
    cntr = cntr + 1;
end

mean_Nb = (lambda_list)./(1-exp(-lambda_list));
figure(1); subplot(1,2,2); plot(mean_Nb, overall_Nb_logL_lowCT_TPS);
loc_max = find(overall_Nb_logL_lowCT_TPS == max(overall_Nb_logL_lowCT_TPS));
est_mean_Nb_lowCT_TPS = mean_Nb(loc_max)
est_lambda_lowCT_TPS = lambda_list(loc_max);

est_pmf_poiss_lowCT_TPS = poisspdf(overall_Nb, est_lambda_lowCT_TPS)/(sum(poisspdf(overall_Nb, est_lambda_lowCT_TPS)));
figure(2); subplot(1,2,2); 
bar(overall_Nb, est_pmf_poiss_lowCT_TPS)
overall_Nb'
est_pmf_poiss_lowCT_TPS'

xlabel('transmission of N_b virions'); ylabel('probability');
title('all TPs with donor iSNVs > 6% and donor CT vals < 30');
save(outfile, 'data', 'var_calling_threshold', 'overall_Nb', 'est_pmf_poiss_allTPS', 'est_pmf_poiss_lowCT_TPS');

