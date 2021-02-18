function void = main_estimateNb_MLE_and_CI(void)

clear all; close all; clc;

load('Popa_data');

var_calling_threshold = 0.06;       % set a variant calling threshold
outfile = 'Popa_data_reanalysis_6percent_NbAnalysis';

% for each TP, find a range of Nb values over which to search for the MLE value of Nb:
for i = 1:n_TPs

    min_Nb = 1; max_Nb = 10000;  % consider bottleneck sizes only between 1 and 10,000

    n_variants = length(find(data(i).donor_iSNVs >= var_calling_threshold));
    data(i).n_variants_at_cutoff = n_variants;
 
    [i n_variants]
    
    if n_variants > 0
        
        while 1
            
            mid_Nb = round((min_Nb + max_Nb)/2);  % find the bottleneck size in the middle of the current min and max values
      
            % get logL at the min, max, and mid Nb values
            min_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, min_Nb);
            max_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, max_Nb);
            mid_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, mid_Nb);
            
            bool_val = 1;
            if isinf(-min_Nb_logL)
                min_Nb = min_Nb + 1;
                bool_val = 0;
            end
            if isinf(-max_Nb_logL)
                max_Nb = round(0.9*max_Nb); %mid_Nb;
                bool_val = 0;
            end
            if bool_val
                break;
            end
        end
        data(i).min_Nb_orig = min_Nb; data(i).min_Nb_logL_orig = min_Nb_logL;
        data(i).max_Nb_orig = max_Nb; data(i).max_Nb_logL_orig = max_Nb_logL;
    end
    save(outfile, 'data', 'var_calling_threshold', 'CT_data');
    [i data(i).min_Nb_orig data(i).max_Nb_orig]
    [data(i).min_Nb_logL_orig data(i).max_Nb_logL_orig]
end

% now, get to Nb MLE for each TP 
for i = 1:n_TPs
    
    min_Nb = data(i).min_Nb_orig;
    max_Nb = data(i).max_Nb_orig;
    
    if data(i).n_variants_at_cutoff > 0
        
        while 1
            mid_Nb = round((min_Nb + max_Nb)/2);
            mid_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, mid_Nb);
            
            midp1_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, mid_Nb+1);
            if mid_Nb > 2
                midm1_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, mid_Nb-1);
            else
                midm1_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, mid_Nb);
            end
            
            if (midp1_Nb_logL >= mid_Nb_logL) % consider only top half
                min_Nb = mid_Nb;
            else % consider only bottom half
                max_Nb = mid_Nb;
            end
            
            if ((max_Nb - min_Nb) == 1)
                min_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, min_Nb);
                max_Nb_logL = GetLogL_forNb(data(i), var_calling_threshold, max_Nb);
                if max_Nb_logL > min_Nb_logL
                    data(i).Nb_MLE = max_Nb;
                    data(i).Nb_MLE_logL = max_Nb_logL;
                else
                    data(i).Nb_MLE = min_Nb;
                    data(i).Nb_MLE_logL = min_Nb_logL;
                end
                break;
            end
            if ((max_Nb - min_Nb) == 0)
                data(i).Nb_MLE = max_Nb;
                data(i).Nb_MLE_logL = max_Nb_logL;
                break;
            end
        end
    else
       data(i).Nb_MLE = NaN;
       data(i).Nb_MLE_logL = NaN; 
    end
    [i data(i).Nb_MLE  data(i).Nb_MLE_logL]
    save(outfile, 'data', 'var_calling_threshold', 'CT_data');
end

% NOW ESTIMATE CONFIDENCE INTERVALS:

for i = 1:n_TPs
    data(i).Nb_low1 = 1;
    data(i).Nb_low2 = data(i).Nb_MLE;
    data(i).Nb_high1 = data(i).Nb_MLE;
    data(i).Nb_high2 = 10000;

    data(i).Nb_low1_logL = GetLogL_forNb(data(i), var_calling_threshold, data(i).Nb_low1);
    
    data(i).Nb_low2_logL = GetLogL_forNb(data(i), var_calling_threshold, data(i).Nb_low2);
    
    data(i).Nb_high1_logL = GetLogL_forNb(data(i), var_calling_threshold, data(i).Nb_high1);
    
    data(i).Nb_high2_logL = GetLogL_forNb(data(i), var_calling_threshold, data(i).Nb_high2);
    
    save(outfile, 'data', 'var_calling_threshold', 'CT_data');
end

display('calculating lower bound');

for i = 1:n_TPs
    
    i
    cutoff_logL = data(i).Nb_MLE_logL - 1.92;
    
    if ~isnan(cutoff_logL)
    
        while 1
            low_mid = round((data(i).Nb_low1 + data(i).Nb_low2)/2);
    
            low_mid_logL = GetLogL_forNb(data(i), var_calling_threshold, low_mid);
        
            if low_mid_logL >= cutoff_logL  % low is still too high
                data(i).Nb_low2 = low_mid;
                data(i).Nb_low2_logL = low_mid_logL;
            else
                data(i).Nb_low1 = low_mid;
                data(i).Nb_low1_logL = low_mid_logL;
            end
            if (data(i).Nb_low2 - data(i).Nb_low1) <= 1
                data(i).Nb_low = data(i).Nb_low2;
                data(i).Nb_low_logL = data(i).Nb_low2_logL;
                break;
            end
        end
    end
    save(outfile, 'data', 'var_calling_threshold', 'CT_data');
end

display('calculating upper bound');

for i = 1:n_TPs
    
    i
    
    cutoff_logL = data(i).Nb_MLE_logL - 1.92;
    
    if ~isnan(cutoff_logL)
            
        while 1
            high_mid = round((data(i).Nb_high1 + data(i).Nb_high2)/2)
    
            high_mid_logL = GetLogL_forNb(data(i), var_calling_threshold, high_mid);
            
            if high_mid_logL >= cutoff_logL  % low is still too high
                data(i).Nb_high1 = high_mid;
                data(i).Nb_high1_logL = high_mid_logL;
            else
                data(i).Nb_high2 = high_mid;
                data(i).Nb_high2_logL = high_mid_logL;
            end
            if (data(i).Nb_high2 - data(i).Nb_high1) <= 1
                data(i).Nb_high = data(i).Nb_high1;
                data(i).Nb_high_logL = data(i).Nb_high1_logL;
                break;
            end
        end
    
    end
    save(outfile, 'data', 'var_calling_threshold', 'CT_data');
end

display('showing results')
for i = 1:n_TPs
    i
    [data(i).Nb_MLE data(i).Nb_low data(i).Nb_high]
    [data(i).Nb_MLE_logL data(i).Nb_low_logL data(i).Nb_high_logL]
end
    
save(outfile, 'data', 'var_calling_threshold', 'CT_data');
    
