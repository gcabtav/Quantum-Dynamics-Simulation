clf;
clc;
clear;
% Define parameters
level = 9;  % Discretization level
tmax = 0.02 ;  % Maximum integration time
lambda = 0.05;  % dt/dx 
idtype = 1; % Initial data type (for Gaussian initial data)

% EXPERIMENTS

% Scattering off rectangular barrier (quantum tunneling)
% idpar = [0.1, 0.5, 0.08, 0.08, 5, 0]; % Initial data parameters
% vtype = 1;  % Potential type 
% vc = exp(6);  % Potential parameter (raise to exp(10) for almost no tunneling)
% vpar = [0, 1, 0.4, 0.6, vc];

% % Scattering off rectangular well
% idpar = [0.1, 0.5, 0.08, 0.08, 5, 0]; % Initial data parameters
% vtype = 1;  % Potential type 
% vc = exp(8);  % Potential parameter
% vpar = [0, 1, 0.4, 0.6, -vc];

% double slit ( | )
idpar = [0.1, 0.5, 0.08, 0.08, 3, 0];
vtype = 2;
vc = exp(9);  % Potential parameter
slit_width = 0.1; % Adjust as needed
vpar = [0.4 - slit_width / 2, 0.4 + slit_width / 2, 0.6 - slit_width / 2, 0.6 + slit_width / 2, vc];


% Self Exploration (bonus)

% % Circular Well (Deep)
% idpar = [0.5, 0.5, 0.08, 0.08, 0, 0]; % Initial data parameters
% vtype = 4;  % Potential type 
% vc = exp(7);  % Potential parameter
% vpar = [0.3, 0.7, 0.3, 0.7, -vc];

% % Circular Well (Shallow)
% idpar = [0.5, 0.5, 0.08, 0.08, 0, 0]; % Initial data parameters
% vtype = 4;  % Potential type 
% vc = exp(2);  % Potential parameter
% vpar = [0.3, 0.7, 0.3, 0.7, -vc];
% 
% % Gaussian Barrier
% idpar = [0.5, 0.5, 0.08, 0.08, 0, 0]; % Initial data parameters
% vtype = 5;  % Potential type 
% vc = exp(8);  % Potential parameter
% vpar = [0.3, 0.7, 0.3, 0.7, vc];
% 
% % Gaussian Deep Well
% idpar = [0.5, 0.5, 0.08, 0.08, 0, 0]; % Initial data parameters
% vtype = 5;  % Potential type 
% vc = exp(8);  % Potential parameter
% vpar = [0.2, 0.8, 0.2, 0.8, -vc];

% Gaussian Shallow Well
% idpar = [0.5, 0.5, 0.08, 0.08, 0, 0]; % Initial data parameters
% vtype = 5;  % Potential type 
% vc = exp(5);  % Potential parameter
% vpar = [0.2, 0.8, 0.2, 0.8, -vc];

% two double slits  ( |  | )
% idpar = [0.1, 0.5, 0.08, 0.08, 5, 0];
% vtype = 3;
% vc = exp(9);  % Potential parameter
% slit_width = 0.1; 
% vpar = [0.4 - slit_width / 2, 0.4 + slit_width / 2, 0.6 - slit_width / 2, 0.6 + slit_width / 2, vc];






% Solve Schrödinger equation (sch_2d_cn)
[x, y, t, psi, psire, psiim, psimod, prob, v] = sch_2d_cn(tmax, level, lambda, idtype, idpar, vtype, vpar);

% Create contour plots for different time steps and save frames
figure;
for i = 1:length(t)
    contourf(x, y, squeeze((abs(psi(i, :, :))).^2)); 
    title(sprintf('Time: %f', t(i)));
    xlabel('y');
    ylabel('x');
    % Customize plot attributes, add labels, etc. (as per your requirement)
    % Capture frame for the movie
    frame = getframe(gcf);
    % Append each frame to an array (frames)
    frames(i) = frame;
    % pause(0.1)
end

% Create AVI movie
movieFilename = 'scenario8.avi';  % Change the filename accordingly
video = VideoWriter(movieFilename, 'Uncompressed AVI');  % Create video object
video.FrameRate = 10;  % Set frame rate (adjust as needed)
open(video);  % Open video for writing

% Write frames to the video
for i = 1:length(frames)
    writeVideo(video, frames(i));
end

close(video);  % Close video file
