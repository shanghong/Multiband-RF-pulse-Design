% Design a spectral selective excitation pulse for bSSFP C-13 sequence at 14T.
% changed: 
%   minorder = 43, 
%   d1 = 0.01; 
%   d2 = 0.005
%   FA = 60; 

clear; clc; close all;
addpath ./rf_tools         % John Pauly's SLR rf pulse design package

% define parameters
B0 = 14;                   % in T, field strength
n = 100;                   % number of samples of RF pulse/filter
T = 4;                     % in ms, pulse duration, will be optimized later
FA = 90;                   % in degree, flip angle
d1 = 0.02;                 % in-slice ripple
d2 = 0.02;                 % out-of-slice ripple
ptype = 'ex';              % type of pulse, excitation pulse 
ftype = 'ap_minorder_cvx'; % type of filter design, arbitrary-phase with minimal order/length using cvx
nucleus = 'C-13';          % application on C-13 MRI
gamma = 1.0705;            % in kHz/G, gamma/2pi for C-13
flip_zero = 0;             % do not use flip-zero method
downsampling = 1;          % downsampling rate for filter design, 1 (no downsampling)
Peak = [];                 % control end-spike (Conolly Wing), default 1e-3
dbg = 3;                   % show debug plot if dbg >= 1
min_order = 1;             % specify how to minimize order
                           % 1 find minimal order
                           % n>1, set the order to this fixed value
                           % n = [n1 n2 ...], try different order, compare the performance (peak-amplitude, total-energy)                          
% define multiband
comp_select = 'pyruvate';      % which compound to select
pick_compound = [6 1 3 4 2];  % which compound to consider, with the order of 
                              % 'pyruvate' 'lactate' 'alanine' 'pyruvate hydrate' 'bicarbonate' 'urea'
mb_range = 0.1*ones(1,5);     % in kHz, bandwidth for each band
[mb_cf,name_cell] = spectrum_C13(B0,dbg);  % resonance frequency for C-13 spectrum
mb_cf = mb_cf(pick_compound); % choose those of interest
name_cell = name_cell(pick_compound); % corresponding name

% =========================================================================
% no parameters need to be set below
% =========================================================================

% create multiband spec
switch comp_select
    case 'pyruvate'
        mb_FA = [0 FA 0 0 0];         % in degree, flip angle of each band 
        mb_ripple = [d2 d1 d2 d2 d2]; % ripple of magnetization          
        mb_cf = mb_cf - mb_cf(2);     % shift selected compound on-resonance
        
    case 'lactate'
        mb_FA = [0 0 0 0 FA];         
        mb_ripple = [d2 d2 d2 d2 d1];         
        mb_cf = mb_cf - mb_cf(5);
        
    case 'urea'
        mb_FA = [FA 0 0 0 0];         
        mb_ripple = [d1 d2 d2 d2 d2];              
        mb_cf = mb_cf - mb_cf(1);
        
    case 'alanine'
        mb_FA = [0 0 FA 0 0];        
        mb_ripple = [d2 d2 d1 d2 d2];              
        mb_cf = mb_cf - mb_cf(3);
        
    otherwise
        error(['Option of ',comp_select,' is not available yet']);
end
mb_cf = num2cell(mb_cf*1e-3) ;% in kHz, center frequency for each band

% check sampling rate for implementation
dt = T/n;            % in ms, sampling interval
dt_epic = 4*1e-3;    % in ms, samping interval in epic/GE scanner
if mod(dt/dt_epic, 1) ~= 0  % T cannot be divided by res*dt
    dt = dt_epic * floor(dt/dt_epic);
    T = n * dt;
end
fs = 1/dt;           % in kHz, sampling frequency

% design the optimal pulse
tic;
name1 = 'mb-ap-SLR';
[rf, ~, rf_spec, b_spec] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg, min_order, [], [], name_cell);                 
ComTime = toc;







