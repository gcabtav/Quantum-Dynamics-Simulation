clc;
clear;
% Experiment 1:Barrier Survey
tmax = 0.10;
level = 9;
lambda = 0.01;
idtype = 1; % (boosted Gaussian)
idpar = [0.40, 0.075, 20.0];
vtype = 1; % (rectangular barrier)

x1 = 0.8;
x2 = 1.0;

nx = 2^level + 1;
x = linspace(0.0, 1.0, nx);
vals = 251; 
lnvc = linspace(-2, 5, vals);

lnFe_array = zeros(1, vals);
for i = 1 : vals
    vc = exp(lnvc(i));
    vpar = [0.6, 0.8, vc];
    [x, ~, ~, ~, ~, ~, prob, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);
    prob_mean_j = mean(prob, 1);
    prob_mean_j = prob_mean_j / prob_mean_j(end);
    Fe = (prob_mean_j(end) - prob_mean_j(411)) / (x(end) - x(411));
    lnFe_array(1, i) = log(Fe);
end


plot(lnvc, lnFe_array, 'r-o')
xlabel('ln(V0)')
ylabel('ln(F¯e(0.8,1.0))')
title('ln(F¯e(0.8,1.0)) versus ln(V0) for Barrier Survey')
grid on;












