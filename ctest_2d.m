clf;
clf;
clf;
clf;
clear;
clc;
% Start the timer
tic;


% Convergence Test 
idtype = 0;
vtype = 0;
idpar = [2, 3];
tmax = 0.05;
lambda = 0.05;
lmin = 6;
lmax = 9;

m1 = idpar(1);
m2 = idpar(2);


level = lmin; 
[x, y, t, psi_6, ~, ~, ~, ~, ~] = sch_2d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);
% disp(1)

level = lmin + 1; 
[x2, y2, ~, psi_7, ~, ~, ~, ~, ~] = sch_2d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);
% disp(2)

level = lmax - 1; 
[x3, y3, ~, psi_8, ~, ~, ~, ~, ~] = sch_2d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);
% disp(3)

level = lmax; 
[x4, y4, ~, psi_9, ~, ~, ~, ~, ~] = sch_2d_cn(tmax, level, lambda, idtype, idpar, vtype, 0);
% disp(4) 

% 2:1 coarsening in space and time dimensions
psi_7 = psi_7(1 : 2 : end, 1 : 2 : end, 1 : 2 : end);
psi_8 = psi_8(1 : 4 : end, 1 : 4 : end, 1 : 4 : end);
psi_9 = psi_9(1 : 8 : end, 1 : 8 : end, 1 : 8 : end);

% ||d(psi^l)||

dpsi_76 = psi_7 - psi_6;
dpsi_76_norm = rms(dpsi_76(:,:).');

dpsi_87 = psi_8 - psi_7;
dpsi_87_norm = rms(dpsi_87(:,:).');

dpsi_8 = psi_9 - psi_8;
dpsi_98_norm = rms(dpsi_8(:,:).');

% re-scaling - to show order accuracy is O(h^2)
dpsi_87_norm_s = 4*dpsi_87_norm;
dpsi_98_norm_s = 4^2*dpsi_98_norm;

figure(1);
hold on;
plot(t, dpsi_76_norm, 'r-o')
plot(t, dpsi_87_norm_s, 'b-o')
plot(t, dpsi_98_norm_s, 'g-o')
title('Scaled ||d\psi^l|| vs Time - 2D Case');
ylabel('Scaled  ||d\psi^l||');
xlabel('time')
hold off;

figure(2);
hold on;
plot(t, dpsi_76_norm, 'r-o')
plot(t, dpsi_87_norm, 'b-o')
plot(t, dpsi_98_norm, 'g-o')
title('||d \psi^l|| vs Time - 2D Case');
ylabel('||d\psi^l||');
xlabel('time')
legend('levels = 6-7', 'levels = 7-8', 'levels = 8-9', location='northwest');
hold off;



% % ||E(psi^l)||
% [Y, T, X] = meshgrid(sin(m2 * pi * y), exp(-1i * (m1^2 + m2^2) * pi^2 * t), sin(m1 * pi * x));
% 
% psi_exact_2d = X.*T.*Y;
% 
% dpsi_6 = psi_exact_2d - psi_6;
% dpsi_6_norm = rms(dpsi_6(:,:).');
% 
% dpsi_7 = psi_exact_2d - psi_7;
% dpsi_7_norm = rms(dpsi_7(:,:).');
% 
% dpsi_8 = psi_exact_2d - psi_8;
% dpsi_8_norm = rms(dpsi_8(:,:).');
% 
% dpsi_9 = psi_exact_2d - psi_9;
% dpsi_9_norm = rms(dpsi_9(:,:).');
% 
% % re-scaling - to show order accuracy is O(h^2)
% dpsi_7_norm_s = 2^2 * dpsi_7_norm;
% dpsi_8_norm_s = 2^4 * dpsi_8_norm;
% dpsi_9_norm_s = 2^6 * dpsi_9_norm;
% 
% figure(1);
% hold on;
% plot(t, dpsi_6_norm, 'r-o')
% plot(t, dpsi_7_norm_s, 'b-o')
% plot(t, dpsi_8_norm_s, 'g-o')
% plot(t, dpsi_9_norm_s, 'k-o')
% title('Scaled ||E(\psi^l)|| vs Time - 2D Case');
% ylabel('Scaled  ||E(\psi^l)||');
% xlabel('time')
% hold off;
% 
% figure(2);
% hold on;
% plot(t, dpsi_6_norm, 'r-o')
% plot(t, dpsi_7_norm, 'b-o')
% plot(t, dpsi_8_norm, 'g-o')
% plot(t, dpsi_9_norm, 'k-o')
% title('||E(\psi^l)|| vs Time - 2D Case');
% ylabel('||E(\psi^l)||');
% xlabel('time')
% legend('level = 6', 'level = 7', 'level = 8','level = 9', Location='northwest');
% hold off;

% Stop the timer and get the elapsed time
elapsed_time = toc;

% Display the elapsed time
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);


