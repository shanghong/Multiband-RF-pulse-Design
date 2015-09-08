function sim_rf_scale(rf, dt, rfname, scale, nucleus, ptype, f, a, d)
%
%  sim_rf_scale(rf, dt, rfname, scale, nucleus, ptype, f, a, d)
%
%  Simulate and display the profile of a 1D RF pulse, based on 
%  Brian Hargreaves' bloch simulator, especially when there is an error of 
%  transmit gain calibration (scaled RF pulse)
%  
%  Inputs:
%    rf -- in Gauss, n*1 array, RF pulse
%    dt -- in ms, sampling interval for rf pulse
%    rfname -- name of this RF pulse 
%    scale -- array of scaling factor
%    nucleus -- 'C-13' or 'H-1', different gamma in bloch simulation
%    ptype -- pulse type, options are st, ex, sat, inv
%        below is defined spec for Mxy or Mz profile
%    f -- in kHz, if 1*2m array, frequency band edge
%                 if scalar, single-band pulse passband bandwidth
%    a -- 1*2m array, amplitude at band edges       [optional]
%    d -- 1*m aarray, ripple of each band           [optional]
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Jun 2014

% check input
if nargin == 7
    show_spec = 0;          % do not show multiband spec
    BW = f*1e3;             % in Hz, passband bandwidth 
    Wrange = [-BW*3, BW*3]; % in Hz, frequency range in simulation
    
elseif nargin == 9
    show_spec = 1;                      % show multiband spec
    f = f*1e3;                          % in Hz
    Wrange = [ f(1)-300, f(end)+300 ];  % in Hz, frequency range in simulation
    
else
    error('Number of input should be either 7 or 9');
    
end

if isempty(scale)
    scale = [0.8 0.9 1 1.1 1.2];
end

rf = rf(:);

% choose proper gamma for bloch simulation
switch nucleus
    case 'H-1'
        bloch_fun_handle = @blochH; 
    case 'C-13'
        bloch_fun_handle = @blochC;
    otherwise
        error('No such option for nucleus. Options are H-1 and C-13 \n');
end

% prepare for plot
leg_name = {};
for kk=1:length(scale)
    leg_name{kk} = [num2str( ( scale(kk)-1 )*100 ), '%'];
end

color_name = {'b','k','r','g','m','y','c'};
if length(scale) > length(color_name)
    error('not enough color option');
end

name_save = [rfname,'_scale'];

% prepare for simulation
N_simu = 256*8;
df = transpose(linspace(Wrange(1), Wrange(2), N_simu)); % in Hz
dp = 0; % in cm
mx0 = zeros(length(dp),length(df));
my0 = zeros(length(dp),length(df));
mz0 = ones(length(dp),length(df));

% run simulation
mxy_array = {};
mz_array = {};
for kk=1:length(scale)
    rf1 = rf*scale(kk);
    G1 = zeros(size(rf));
    [mx1,my1,mz1] = bloch_fun_handle(rf1,G1,dt*1e-3,1e3,1e3,df,dp,0,mx0,my0,mz0);
    mxy1 = mx1 + sqrt(-1)*my1;
    mxy_array{kk} = mxy1;
    mz_array{kk} = mz1;
end
    

% plot rf pulse 
figure;
plot((1:length(rf))*dt, real(rf),'r-*',(1:length(rf))*dt, imag(rf),'b-*'); leg = legend('real','imag'); set(leg,'FontSize',17); set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, (length(rf)+1)*dt]); 
set(gcf, 'Position', [50, 100, 700, 400], 'PaperPositionMode', 'auto');
hgexport(gcf, [name_save, '_RF.eps']);

% plot mz
figure;  
subplot(1,2,1); 
for kk=1:length(scale)
    plot(df, mz_array{kk}, color_name{kk}); hold on; 
end
hold off; leg = legend(leg_name); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18);  ylabel('Mz','FontSize',18);  set(gca,'FontSize',18);  axis tight;  
if ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
    hold on; plot_spec(f,a,d); hold off;
end

subplot(1,2,2);
for kk=1:length(scale)
    semilogy(df, 1-mz_array{kk}, color_name{kk}); hold on; 
end
hold off; leg = legend(leg_name); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18); ylabel('1-Mz','FontSize',18);  set(gca,'FontSize',18); axis tight; 
if ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
    hold on; plot_spec_log(f,1-a,d); hold off;
end
set(gcf, 'Position', [150, 100, 1000, 400], 'PaperPositionMode', 'auto');
%ax = axes('position',[0,0,1,1],'visible','off');  
%text(0.07,0.95, '(A)','FontSize',23);
%text(0.51,0.95, '(B)','FontSize',23);
hgexport(gcf, [name_save, '_Mz.eps']);


% plot mxy magnitude
figure;
for kk=1:length(scale)
    plot(df, abs(mxy_array{kk}), color_name{kk}); hold on; 
end
hold off; leg = legend(leg_name);  set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18);  ylabel('|Mxy|','FontSize',18);  set(gca,'FontSize',18); axis tight;  set(gca,'ylim',[0 1+0.04]);
if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
    hold on; plot_spec(f,a,d); hold off;
end
set(gcf, 'Position', [250, 100, 700, 400], 'PaperPositionMode', 'auto');
hgexport(gcf, [name_save, '_Mxy.eps']);


    
% % plot mxy in log scale
% figure;
% for kk=1:length(scale)
%     semilogy(df, abs(mxy_array{kk}), color_name{kk}); hold on; 
% end
% leg = legend(leg_name); set(leg,'FontSize',16);  xlabel('Frequency, Hz','FontSize',17); set(gca,'FontSize',16); axis tight; title('Mxy Magnitude','FontSize',20); 
% if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
%     hold on; plot_spec_log(f,a,d); hold off; 
% end
% set(gcf, 'Position', [650, 500, 600, 400], 'PaperPositionMode', 'auto');
% print(gcf,'-depsc2',[name_save, '_Mxy_log.eps']); 
    
