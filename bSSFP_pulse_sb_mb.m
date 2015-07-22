% Design a spectral selective excitation pulse for bSSFP C-13 sequence at 14T.
% Compare multiband and single-band profile 

clear; clc; close all;

addpath ./rf_tools          % John Pauly's SLR rf pulse design package

% define parameters
B0 = 14;                   % in T, field strength
n = 100;                   % number of samples of RF pulse/filter
T = 4;                     % in ms, pulse duration, will be optimized later
FA = 60;                   % in degree, flip angle
d1 = 0.01;                 % in-slice ripple
d2 = 0.005;                % out-of-slice ripple
ptype = 'ex';              % type of pulse, excitation pulse 
ftype = 'ap_minorder_cvx'; % type of filter design, arbitrary-phase with minimal order/length using cvx
nucleus = 'C-13';          % application on C-13 MRI
gamma = 1.0705;            % in kHz/G, gamma/2pi for C-13
flip_zero = 0;             % do not use flip-zero method
downsampling = 1;          % downsampling rate for filter design, 1 (no downsampling)
Peak = [];                 % control end-spike (Conolly Wing), default 1e-3
dbg = 2;                   % show debug plot if dbg >= 1
write_pulse = 0;           % write pulse in a file with Varian format

min_order = 43;            % specify how to minimize order
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

% design the pulse
if isscalar(min_order)  
    % optimal pulse design
    tic;
    name1 = 'mb-ap-SLR';
    [rf, ~, rf_spec, b_spec] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg, min_order);                 
    ComTime = toc;
    
    % gold standard of pm(equal-ripple) minimum phase pulse
    tic; 
    BW = 1.3;     % in kHz, passband bandwidth
    T_new = 1.85;     % in ms, pulse duration
    dt2 = T_new/n; % in ms, sampling interval
    name2 = 'sb-mp-SLR';
    rf2 = dzrf(n, T_new*BW, ptype, 'min', d1, d2);  % design the standard slr pulse 
    rf2 = rf2*(FA*pi/180)/sum(rf2);    % scale the pulse
    rf2 = rfscaleg(rf2, T_new, gamma); % scale the pulse in Gauss
    ComTime2 = toc;

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
   
    % write to external waveform with Varian format
    root_fname = ['mbslr_', comp_select,'_',num2str(FA),'deg_',num2str(length(rf)*dt),'ms'];
    if write_pulse       
        rfwrite_varian(rf, length(rf)*dt, [], root_fname);
        rfwrite_varian(rf2, length(rf2)*dt2, [], 'rf_reference');
        
        rf_spec_f = rf_spec.f*(fs/2); % in kHz
        save([root_fname, '_d1_', num2str(d1), '_d2_', num2str(d2), '.mat'],'rf','dt','rf_spec','rf_spec_f','FA','d1','d2');
    end
    
    
% manually choose order considering peak B1 and SAR
else 
   
    peak_b1_array = zeros(size(min_order));
    SAR_array = zeros(size(min_order));
    duration_array = zeros(size(min_order));
    comTime_array = zeros(size(min_order));
    
    for i=1:length(min_order)
        tic;
        [rf, ~, ~, ~] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, 0, min_order(i));                 
        comTime_array(i) = toc;
        
        if ~isempty(rf)
            peak_b1_array(i) = max(abs(rf));
            SAR_array(i) = sum((abs(rf)).^2)*dt;
            duration_array(i) = length(rf)*dt;
        else
            peak_b1_array(i) = 0;
            SAR_array(i) = 0;
            duration_array(i) = min_order(i)*dt;
        end
    end
    
    figure;
    subplot(3,1,1); plot(duration_array, SAR_array, 'r*-');      xlabel('duration (ms)','FontSize',17); ylabel('SAR (G^2*ms)','fontsize',17);         axis tight; 
    subplot(3,1,2); plot(duration_array, peak_b1_array, 'r*-');  xlabel('duration (ms)','FontSize',17); ylabel('Peak B1 (G)','fontsize',17);          axis tight; 
    subplot(3,1,3); plot(duration_array, comTime_array, 'r*-');  xlabel('duration (ms)','FontSize',17); ylabel('computation time (s)','fontsize',17); axis tight; 
      
end






