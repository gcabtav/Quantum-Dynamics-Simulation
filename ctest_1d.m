clf;
clf;
clf;
clf;
clear;
clc;
% Convergence Testing 2

% idtype = 1;
% vtype = 0;
% idpar = [0.5 0.075 0.0];
% tmax = 0.01;
% lambda = 0.01;
% lmin = 6;
% lmax = 9;


% Convergence Testing 1

idtype = 0;
vtype = 0;
idpar = 3;
tmax = 0.25;
lambda = 0.1;
lmin = 6;
lmax = 9;

level = lmin; 
[x1, t1, psi_6, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);


level = lmin + 1; 
[~, ~, psi_7, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);


level = lmax - 1; 
[~, ~, psi_8, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);
 
level = lmax; 
[~, ~, psi_9, ~, ~, ~, ~, ~] = sch_1d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);
 

% 2:1 coarsening in space and time dimensions
psi_7 = psi_7(1 : 2 : end, 1 : 2 : end);
psi_8 = psi_8(1 : 4 : end, 1 : 4 : end);
psi_9 = psi_9(1 : 8 : end, 1 : 8 : end);

% d(psi^l)

dpsi_76 = psi_7 - psi_6;
dpsi_76_norm = rms(dpsi_76.');

dpsi_87 = psi_8 - psi_7;
dpsi_87_norm = rms(dpsi_87.');

dpsi_98 = psi_9 - psi_8;
dpsi_98_norm = rms(dpsi_98.');

% re-scaling - to show order accuracy is O(h^2)
dpsi_87_norm_s = 4 * dpsi_87_norm;
dpsi_98_norm_s = 4^2 * dpsi_98_norm;

figure(1);
hold on;
plot(t1, dpsi_76_norm, 'r-o')
plot(t1, dpsi_87_norm, 'b-o')
plot(t1, dpsi_98_norm, 'g-o')
title('||d\psi^l|| vs Time - 1D Case');
ylabel('||d\psi^l||');
xlabel('time')
legend('levels = 6-7', 'levels = 7-8', 'levels = 8-9', location='northwest');
hold off;

figure(2);
hold on;
plot(t1, dpsi_76_norm, 'r-o')
plot(t1, dpsi_87_norm_s, 'b-o')
plot(t1, dpsi_98_norm_s, 'g-o')
title('Scaled ||d\psi^l|| vs Time - 1D Case');
ylabel('Scaled ||d\psi^l||');
xlabel('time');
hold off;



% % E(psi^l)
% 
% % exact solution
% m = idpar;
% psi_exact = @(t,x) sin(m*pi*x) .* exp(-1i * m^2 * pi^2 * t).';
% 
% % Difference
% dpsi_6 = psi_exact(t1, x1) - psi_6;
% dpsi_6_norm = rms(dpsi_6.');
% 
% dpsi_7 = psi_exact(t1, x1) - psi_7;
% dpsi_7_norm = rms(dpsi_7.');
% 
% dpsi_8 = psi_exact(t1, x1) - psi_8;
% dpsi_8_norm = rms(dpsi_8.');
% 
% dpsi_9 = psi_exact(t1, x1) - psi_9;
% dpsi_9_norm = rms(dpsi_9.');
% 
% % re-scaling - to show order accuracy is O(h^2)
% dpsi_7_norm_s = 2^2 * dpsi_7_norm;
% dpsi_8_norm_s = 2^4 * dpsi_8_norm;
% dpsi_9_norm_s = 2^6 * dpsi_9_norm;
% 
% figure(1);
% hold on;
% plot(t1, dpsi_6_norm, 'r-o')
% plot(t1, dpsi_7_norm, 'b-o')
% plot(t1, dpsi_8_norm, 'g-o')
% plot(t1, dpsi_9_norm, 'k-o')
% title('||E(\psi^l)|| vs Time - 1D Case');
% ylabel('||E(\psi^l)||');
% xlabel('time')
% legend('level = 6', 'level = 7', 'level = 8','level = 9', Location='northwest');
% hold off;
% 
% figure(2);
% hold on;
% plot(t1, dpsi_6_norm, 'r-o')
% plot(t1, dpsi_7_norm_s, 'b-o')
% plot(t1, dpsi_8_norm_s, 'g-o')
% plot(t1, dpsi_9_norm_s, 'k-o')
% title('Scaled ||E(\psi^l)|| vs Time - 1D Case');
% ylabel('Scaled  ||E(\psi^l)||');
% xlabel('time')
% hold off;
% 
% 
% 
