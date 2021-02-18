function [logL_vals] = GetLogL_forNb(data, var_calling_threshold, Nb_vals)

locsdonorbelow = find(data.donor_iSNVs < var_calling_threshold);
data.donor_iSNVs(locsdonorbelow) = 0;
locsdonorabove = find(data.donor_iSNVs > (1-var_calling_threshold));
data.donor_iSNVs(locsdonorabove) = 1;

locsrecipientbelow = find(data.recipient_iSNVs < var_calling_threshold);
data.recipient_iSNVs(locsrecipientbelow) = 0;
locsrecipientabove = find(data.recipient_iSNVs > (1-var_calling_threshold));
data.recipient_iSNVs(locsrecipientabove) = 1;

locs = intersect(find(data.donor_iSNVs >= var_calling_threshold), find(data.donor_iSNVs <= (1-var_calling_threshold)));
n_variants = length(locs);
if n_variants == 0
    logL_vals = NaN*ones(size(Nb_vals));
    return;
end

donor_freqs_observed = data.donor_iSNVs(locs);
recipient_freqs_observed = data.recipient_iSNVs(locs);

log_likelihood_matrix = NaN*zeros(n_variants, length(Nb_vals)); % keeps log-likelihood of each var site for each Nb value

for var = 1:n_variants
    if recipient_freqs_observed(var) == 1 
        recipient_freqs_observed(var) = 0; 
        donor_freqs_observed(var) = 1-donor_freqs_observed(var);
    end    
    log_likelihood_matrix(var,:) = GetBetaBinomialApproxLogLikelihoodAtSite(donor_freqs_observed(var), recipient_freqs_observed(var), Nb_vals, var_calling_threshold);
end
logL_vals = sum(log_likelihood_matrix, 1);
