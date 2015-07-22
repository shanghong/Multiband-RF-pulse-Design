% design a dualband spectral saturation pulse for H-1 MRS at 3T
clear; clc; close all;

% define parameters
n = 260;                     % number of samples
B0 = 127794577/(42.577*1e6); % in T, field strength
T = 26;                      % in ms, duration
d1 = 0.05;                   % in slice ripple
d2 = 0.001;                  % out of slice ripple
ptype = 'sat';               % pulse type 
ftype = 'ap_mintran_cvx';    % FIR filter design type
rfname = 'db-specsat';       % name for saving figures
nucleus = 'H-1';
gamma = 4.2576;              % in kHz/G, gamma/2pi for H-1
flip_zero = 0;               % do not use flip-zero method
downsampling = 1;            % downsampling rate for filter design
Peak = [];                   % control end-spike (Conolly Wing), default 1e-3
dbg = 1;                     % show debug plot if dbg >= 1
min_tran = 0.95;             % in [0 1], to what extent to minimize transition width
shift_f = 1;                 % centralize frequency during RF design
write_pulse = 0;             % write pulse in a file for GE system

% define multiband
mb_cf1 = [1.8 2.5]; % in ppm
mb_cf2 = [3   4.1]; % in ppm
mb_cf3 = [4.8 5.4]; % in ppm, this is water
mb_cfr = mean(mb_cf3);
mb_cf = { (mb_cf1-mb_cfr)*B0*42.577*1e-3   (mb_cf2-mb_cfr)*B0*42.577*1e-3  (mb_cf3-mb_cfr)*B0*42.577*1e-3 };  % in kHz

mb_FA = [120 0 90];             % in degree, flip angle for each band
mb_range = [0.01  0.01  0.01];  % in kHz, frequency range for each band
mb_ripple = [d1 d2 d1];         % ripple of magnetization for each band


% =========================================================================
% no parameters need to be set below
% =========================================================================

% check sampling rate for implementation on GE scanner
dt = T/n;            % in ms, sampling interval
dt_epic = 4*1e-3;    % in ms, samping interval in epic/GE scanner
if mod(dt/dt_epic, 1) ~= 0  % T cannot be divided by res*dt
    dt = dt_epic * floor(dt/dt_epic);
    T = n * dt;
end
fs = 1/dt;           % in kHz, sampling frequency

% pulse design
tic;
[rf, ~, rf_spec, ~] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg, [], min_tran, shift_f);
ComTime = toc;

% display pulse performance
fprintf('\n');
fprintf('computation time: %0.4f s\n',  ComTime);
fprintf('pulse duration:   %0.3f ms\n', length(rf)*dt);
fprintf('total power:      %0.4f G^2*ms\n', sum((abs(rf)).^2)*dt);
fprintf('peak amplitude:   %0.4f G\n',  max(abs(rf)));

% simulation test
sim_rf_spectral(rf(:),dt,rfname,nucleus,ptype, rf_spec.f*(fs/2),rf_spec.a,rf_spec.d)

% write to external waveform
if write_pulse
    fprintf('dual band pulse 1 \n');
    rfwrite(rf, length(rf)*dt*1e-3, mb_FA(3)*pi/180, gamma*1e3, 0);
end





