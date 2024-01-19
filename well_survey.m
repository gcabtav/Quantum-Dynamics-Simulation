clc;
clear;
% clf;
% Experiment 1:Barrier Survey

tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1; % (boosted Gaussian)
idpar = [0.40, 0.075, 0.0];
vtype = 1; % (rectangular barrier)

x1 = 0.6;
x2 = 0.8;

nx = 2^level + 1;
x = linspace(0.0, 1.0, nx);
vals = 251;
lnvc = linspace(2, 10, vals);

lnFe_array = zeros(1, vals);
for i = 1 : vals
    vc = exp(lnvc(i));
    vpar = [0.6, 0.8, -vc];
    [x, t, psi, psire, psiim, psimod, prob, v] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    prob_mean_j = mean(prob, 1);
    prob_mean_j = prob_mean_j/prob_mean_j(end);
    Fe = (prob_mean_j(411) - prob_mean_j(308))/(x(411)-x(308));
    lnFe_array(1, i) = log(Fe);
end

plot(lnvc, lnFe_array, 'r-o')
xlabel('ln(V0)')
ylabel('ln(F¯e(0.6,0.8))')
title('ln(F¯e(0.6,0.8)) versus ln(V0) for Well Survey')
grid on;
