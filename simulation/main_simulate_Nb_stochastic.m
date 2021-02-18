function void = main_simulate_Nb_stochastic(void)

clear all; close all; clc;

Nb = 1000;
xlist = 0.06*rand(1,1000);
cntr = 1;
for x = xlist
    k_variant = binornd(Nb, x);
    freq_variant(cntr,:) = [x betarnd(k_variant, Nb-k_variant)];
    cntr = cntr + 1;
end

subplot(2,1,1); 
plot(freq_variant(:,1), freq_variant(:,2), 'r.'); hold on;
xlabel('iSNV frequency in donor');
ylabel('iSNV frequency in recipient');
axis([0 0.06 0 0.06]);
plot([0 0.06], [0 0.06], 'r--')

subplot(2,1,2); 
xvector = 0:0.0001:0.06;
prob_transmission = 1- (1-xvector).^1000;
plot(xvector, prob_transmission, 'r'); hold on;
xlabel('iSNV frequency in donor');
ylabel('probability of iSNV transmission to recipient');
axis([0 0.06 0 1.05])

freq_variant