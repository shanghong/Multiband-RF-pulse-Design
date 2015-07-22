% Design a proton water(fat) saturation pulse at 3T
% Compare to a standard maximum-phase SLR pulse

clear; clc; close all;
addpath ./rf_tools          % John Pauly's SLR rf pulse design package

% define parameters
B0 = 3;                    % in T, field strength
dB0ppm = 0.4;              % in ppm, B0 field inhomogeneity, [-0.4 0.4]ppm
cs_wf = 3.5;               % in ppm, chemical shift between water and fat
BW = 0.6;                  % kHz, bandwidth for conventional single-band pulse
n = 200;                   % number of samples of RF pulse/filter
T = 8;                    % in ms, pulse duration
FA = 90;                   % in degree, flip angle
d1 = 0.05;                 % in-slice ripple
d2 = 1e-3;                 % out-of-slice ripple
ptype = 'sat';             % type of pulse, excitation pulse 
ftype = 'ap_minstopripple_cvx';  % type of filter design, arbitrary-phase with minimal stopband ripple
nucleus = 'H-1';           
gamma = 4.2576;            % in kHz/G, gamma/2pi for H-1
flip_zero = 0;             % do not use flip-zero method
downsampling = 2;          % downsampling rate for filter design
Peak = 1e-2;               % control end-spike (Conolly Wing), default 1e-3
dbg = 1;                   % show debug plot if dbg >= 1


% =========================================================================
% no parameters need to be set below
% =========================================================================

% define multiband
f0 = gamma*B0*cs_wf*10;        % in Hz, water/fat frequency shift
mb_cf = [0 f0];        
mb_cf = num2cell(mb_cf*1e-3) ; % in kHz, center frequency for each band
name_cell = {'fat','water'}; % corresponding name
mb_FA = [FA 0];        
mb_ripple = [d1 d2];              
mb_range = gamma*B0*dB0ppm*1e-2*2*ones(1,2); % in kHz, bandwidth for each band


% check sampling rate for implementation
dt = T/n;            % in ms, sampling interval
dt_epic = 4*1e-3;    % in ms, samping interval in epic/GE scanner
if mod(dt/dt_epic, 1) ~= 0  % T cannot be divided by res*dt
    dt = dt_epic * floor(dt/dt_epic);
    T = n * dt;
end
fs = 1/dt;           % in kHz, sampling frequency

% design the pulse
tic;
name1 = 'mb-ap-SLR';
[rf, ~, rf_spec, b_spec] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg);                 
ComTime = toc;

% gold standard of pm(equal-ripple) minimum phase pulse
tic; 
name2 = 'sb-mp-SLR';
rf2 = dzrf(n, T*BW, ptype, 'max', d1, d2);  % design the standard slr pulse 
rf2 = rf2*(FA*pi/180)/sum(rf2);             % scale the pulse
rf2 = rfscaleg(rf2, T, gamma);              % scale the pulse in Gauss
ComTime2 = toc;
dt2 = dt;       % in ms, sampling interval

% compare pulse performance
fprintf('\n');
fprintf('%s vs. %s \n', name1, name2);
fprintf('%s computation time: %0.4f s\n', name1, ComTime)
fprintf('%s computation time: %0.4f s\n', name2, ComTime2)

fprintf('\n');
fprintf('%s pulse duration: %0.3f ms\n', name1, length(rf)*dt);
fprintf('%s pulse duration: %0.3f ms\n', name2, length(rf2)*dt2);

fprintf('\n');
fprintf('%s pulse power: %0.4f G^2*ms\n', name1, sum((abs(rf)).^2)*dt);
fprintf('%s pulse power: %0.4f G^2*ms\n', name2, sum((abs(rf2)).^2)*dt2);
fprintf('power ratio: %1.4f \n', (sum((abs(rf)).^2)*dt) / (sum((abs(rf2)).^2)*dt2) );

fprintf('\n');
fprintf('%s pulse peak amplitude: %0.4f G\n', name1,  max(abs(rf)));
fprintf('%s pulse peak amplitude: %0.4f G\n', name2,  max(abs(rf2)));
fprintf('amplitude ratio: %1.4f \n', (max(abs(rf)) /  max(abs(rf2))) );

% simulate rf pulse profile
rf_cell{1} = rf(:); rf_cell{2} = rf2(:);
rfname_cell{1} = name1; rfname_cell{2} = name2;
sim_rf_spectral(rf_cell, [dt, dt2], rfname_cell, nucleus, ptype, rf_spec.f*(fs/2), rf_spec.a, rf_spec.d, name_cell)

