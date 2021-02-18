function logL = GetLogLForMeanNb(pmf_poiss, Nblist, TP_logL_by_Nb)

if isnan(TP_logL_by_Nb)
    logL = NaN;
    return;
end

likelihood_by_Nb = exp(TP_logL_by_Nb);

cntr = 1; overall_likelihood = 0;
for Nb = Nblist
    overall_likelihood = overall_likelihood + pmf_poiss(cntr)*likelihood_by_Nb(cntr);
    cntr = cntr + 1;
end
logL = log(overall_likelihood);

%TP_logL_by_Nb(1)

    