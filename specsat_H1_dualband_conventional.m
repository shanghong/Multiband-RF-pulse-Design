% design a dualband spectral saturation pulse by adding two filters
% followed by inverse SLR transform
% Manually modify ripple to meet design spec

clear; clc; close all;
addpath ./rf_tools

% define parameters
n = 260;        % number of samples
B0 = 3.0015;    % in T, field strength
T = 26;         % in ms, duration
dt = T/n;       % in ms
d1 = 0.0008;     % in slice ripple
d2 = 0.03;      % out of slice ripple
FA1 = 120;      % in degree, FA of NAA
FA2 = 90;       % in degree, FA of water
ptype = 'sat';  % pulse type 
nucleus = 'H-1';
gamma = 4.2576; % in kHz/G, gamma/2pi for H-1
dbg = 1;        % show debug plot if dbg >= 1

% define dualband frequency
f1 = [1.8 2.5]; % in ppm, this is NAA
f2 = [3   4.1]; % in ppm
f3 = [4.8 5.4]; % in ppm, this is water

f1 = f1*B0*42.577*1e-3; % in kHz
f2 = f2*B0*42.577*1e-3; % in kHz
f3 = f3*B0*42.577*1e-3; % in kHz

fr = mean(f3);
f1 = f1 - fr;
f2 = f2 - fr;
f3 = f3 - fr;

BW1 = ( (f1(2)+f2(1))/2 - mean(f1) ) * 2; % in kHz
BW2 = ( mean(f3) - (f2(2)+f3(1))/2 ) * 2; % in kHz

% design two beta polynomial
b1 = dzmp(n, T*BW1, d1, d2);
b1 = sin(0.5*FA1*pi/180)*b1;

b2 = dzmp(n, T*BW2, d1, d2);
b2 = sin(0.5*FA2*pi/180)*b2;

% combine two beta polnomial
shift_f = mean(f1);    % in kHz
t_axis = (0:(n-1))*dt;     % in ms
b1 = b1 .* exp(-sqrt(-1)*2*pi* shift_f * t_axis);
b = b1 + b2;

if dbg >= 1
    Nsamp = 1e4;
    B1 = ifftshift(fft(b1,Nsamp));
    B2 = ifftshift(fft(b2,Nsamp));
    B = ifftshift(fft(b,Nsamp));
    wdisp = linspace(-1,1,Nsamp);
    
    figure; 
    subplot(2,3,1); plot(1:n,real(b1),'b-', 1:n,imag(b1),'r-', 'linewidth',2); leg = legend('real','imag'); set(leg,'FontSize',18); set(gca,'FontSize',18); title('Beta polynomial b1','FontSize',19); axis tight;
    subplot(2,3,2); plot(1:n,real(b2),'b-', 1:n,imag(b2),'r-', 'linewidth',2); leg = legend('real','imag'); set(leg,'FontSize',18); set(gca,'FontSize',18); title('Beta polynomial b2','FontSize',19); axis tight;
    subplot(2,3,3); plot(1:n,real(b),'b-',  1:n,imag(b),'r-', 'linewidth',2);  leg = legend('real','imag'); set(leg,'FontSize',18); set(gca,'FontSize',18); title('b = b1 + b2','FontSize',19);        axis tight;
    
    subplot(2,3,4); plot(wdisp, abs(B1), 'linewidth',1.5); xlabel('Normalized Frequency','FontSize',18); ylabel('|B_N(z)|','FontSize',18); set(gca,'FontSize',18); axis tight;
    subplot(2,3,5); plot(wdisp, abs(B2), 'linewidth',1.5); xlabel('Normalized Frequency','FontSize',18); ylabel('|B_N(z)|','FontSize',18); set(gca,'FontSize',18); axis tight;
    subplot(2,3,6); plot(wdisp, abs(B), 'linewidth',1.5);  xlabel('Normalized Frequency','FontSize',18); ylabel('|B_N(z)|','FontSize',18); set(gca,'FontSize',18); axis tight;
       
    set(gcf, 'Position', [50, 50, 1000, 750], 'PaperPositionMode', 'auto');
end

% design RF pulse
a = b2a(b); 
rf = ab2rf(a,b);
rf = rfscaleg(rf, T, gamma);
rf = rf(:);

% simulation
fs = 1/dt;
rfname = 'sum-two-filter-SLR';
load('dualband_specsat_spec');
sim_rf_spectral(rf, dt, rfname, nucleus, ptype, rf_spec.f*(10/2),rf_spec.a,rf_spec.d);




