% for bSSFP RF pulse, compare pulse with different end-spike limit
clear; clc; close all;

addpath ./bloch_simulation

% define parameters
B0 = 14;                   % in T, field strength
n = 100;                   % number of samples of RF pulse/filter
T = 4;                     % in ms, initial pulse duration
FA = 60;                   % in degree, flip angle
d1 = 0.01;                 % in-slice ripple
d2 = 0.005;                % out-of-slice ripple
ptype = 'ex';              % type of pulse, excitation pulse 
ftype = 'ap_minorder_cvx'; % type of filter design, arbitrary-phase with minimal order using cvx
nucleus = 'C-13';          % application on C-13 MRI
gamma = 1.0705;            % in kHz/G, gamma/2pi for C-13
flip_zero = 0;             % do not use flip-zero method
downsampling = 1;          % downsampling rate for filter design, 1 = no downsampling
dbg = 2;                   % show debug plot if dbg >= 1
comp_select = 'lactate';   % which compound to select
min_order = 1;             % find minimal order [default]

% define multiband
pick_compound = [6 1 3 4 2];  % which compound to consider, with the order of 
                              % 'pyruvate' 'lactate' 'alanine' 'pyruvate hydrate' 'bicarbonate' 'urea'
mb_range = 0.1*ones(1,5);     % in kHz, bandwidth for each band
[mb_cf,name_cell] = spectrum_C13(B0,dbg);  % resonance frequency for C-13 spectrum
mb_cf = mb_cf(pick_compound); % choose those of interest
name_cell = name_cell(pick_compound); % corresponding name

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

% design pulse
Peak_array = [1e-2 1e-3 1e-4];
rf_array =  {};
name_array = {};
Num = length(Peak_array);
for i = 1 : Num 
    Peak = Peak_array(i);
    [rf, ~, rf_spec, b_spec] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg, min_order);
    rf_array{i} = rf; 
    name_array{i} = sprintf('Spike limit %0.4f',Peak);
end

% display
% plot rf pulse amplitude in only one plot
disp_color = {'k.-','b.-','r.-','y.-'};
figure; 
subplot(1,2,1); set(gca, 'OuterPosition', [0.02, 0.02,  0.41, 0.96]); 
for i=1:Num
    rf = rf_array{i};
    plot((1:length(rf))*dt, abs(rf), disp_color{i}, 'markersize',13); hold on;
end
hold off; leg = legend(name_array);  set(leg,'FontSize',17, 'position',[0.2 0.77 0.1 0.1]);  
set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0*dt, (length(rf)+1)*dt]); set(gca,'ylim',[0 max(abs(rf_array{1}))*1.01]);  


% simulated profile
disp_color = {'k-','b-','r-','y-'};
f_spec = rf_spec.f*(fs/2);
f_spec = f_spec * 1e3; % in Hz
Wrange = [ f_spec(1)-500, f_spec(end)+500 ];  % in Hz

N_simu = 256*4;
df = transpose(linspace(Wrange(1), Wrange(2), N_simu)); % in Hz
dp = 0; % in cm
mx0 = zeros(length(dp),length(df));
my0 = zeros(length(dp),length(df));
mz0 = ones(length(dp),length(df));

%figure;
subplot(1,2,2); set(gca, 'OuterPosition', [0.47, 0.02,  0.53, 0.96]);
for i=1:Num
    rf = rf_array{i};
    rf = rf(:); 
    G = zeros(length(rf),1);
    [mx,my,mz] = blochC(rf,G,dt*1e-3,1e3,1e3,df,dp,0,mx0,my0,mz0);
    mxy = mx + sqrt(-1)*my;
    
    % plot mxy in log scale  
    semilogy(df,abs(mxy),disp_color{i}, 'linewidth',1); hold on;
end
hold off; leg = legend(name_array); set(leg,'FontSize',17, 'position',[0.81 0.17 0.1 0.1]);  
hold on; plot_spec_log(f_spec, rf_spec.a, rf_spec.d, name_cell); hold off;
xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('|Mxy|','FontSize',18); 

ax = axes('position',[0,0,1,1],'visible','off');  
text(0.06,0.96, '(A)','FontSize',21);
text(0.51,0.96, '(B)','FontSize',21);
    
set(gcf, 'Position', [300, 100, 1100, 400], 'PaperPositionMode', 'auto');
hgexport(gcf, 'bSSFP_pulse_lactate_diff_peak.eps');

