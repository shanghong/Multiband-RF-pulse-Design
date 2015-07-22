function sim_rf_spectral(rf,dt,rfname,nucleus,ptype, f,a,d,name_cell)
%
%  sim_rf_spectral(rf,dt,rfname,nucleus, f,a,d,name_cell, ptype)
%
%  Simulate and display the profile of a 1D spectral selective RF pulse,
%  based on Brian Hargreaves' bloch.mex
%  If input two rf pulses, simulate both and display them together for
%  comparison
%
%  Inputs:
%    rf -- in Gauss, RF pulse
%          n*1 array if input one RF pulse, 1*2 cell array if input two RF pulses
%          each one could be complex
%    dt -- in ms, sampling interval for rf pulse, scalar if only one pulse, 1*2 array if two
%    rfname -- name of RF pulse, string if only one pulse, 1*2 cell array if two
%    nucleus -- 'C-13' or 'H-1', different gamma in bloch simulation
%    ptype -- pulse type, options are st, ex, sat, inv
%        below is defined spec for Mxy or Mz profile (similar as cfirpm)
%    f -- in kHz, if 1*2m array, frequency band edge
%                 if scalar, single-band pulse passband bandwidth
%    a -- 1*2m array, amplitude at band edges       [optional]
%    d -- 1*m aarray, ripple of each band           [optional]
%    name_cell -- 1*m cell array, name of each band [optional]
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Jun 2014
%
%  Modified by Hong Shang for switch between C13/H1 simulation Jun 2015

% use Brian Hargreaves's block simulation function
addpath /Users/hshang/matlab/rfpulses/mbSLR/bloch_simulation  

if nargin == 6
    show_spec = 0;          % do not show multiband spec
    BW = f*1e3;             % in Hz, passband bandwidth 
    Wrange = [-BW*3, BW*3]; % in Hz, frequency range in simulation
    
elseif nargin >= 8
    show_spec = 1;                      % show multiband spec
    f = f*1e3;                          % in Hz
    Wrange = [ f(1)-500, f(end)+500 ];  % in Hz, frequency range in simulation
    
    if nargin == 8
        name_cell = [];
    end
    
else
    error('Number of input should be either 4 or 6');
    
end

% choose proper gamma for bloch simulation
switch nucleus
    case 'H-1'
        bloch_fun_handle = @blochH; 
    case 'C-13'
        bloch_fun_handle = @blochC;
    otherwise
        error('No such option for nucleus. Options are H-1 and C-13 \n');
end

% define parameters for bloch simulation
N_simu = 256*8;
df = transpose(linspace(Wrange(1), Wrange(2), N_simu)); % in Hz
dp = 0; % in cm
mx0 = zeros(length(dp),length(df));
my0 = zeros(length(dp),length(df));
mz0 = ones(length(dp),length(df));

% =========================================================================
% bloch simulation if only input one RF pulse
% =========================================================================
if ( (~iscell(rf)) && (size(rf,2) == 1) )  
    name_save = [rfname,'-',ptype,'-pulse'];
    G = zeros(size(rf,1),1);
    
    [mx,my,mz] = bloch_fun_handle(rf,G,dt*1e-3,1e3,1e3,df,dp,0,mx0,my0,mz0);
    mxy = mx + sqrt(-1)*my;
    
    % plot rf pulse 
    figure;
    plot((1:length(rf))*dt, real(rf),'r-', (1:length(rf))*dt, imag(rf),'b-','linewidth',2.5); leg = legend('real','imag'); set(leg,'FontSize',18,'position',[0.2 0.8 0.1 0.1]); set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[-2*dt, (length(rf)+3)*dt]);  % title('RF pulse','FontSize',20);
    set(gcf, 'Position', [50, 100, 700, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_RF_real_imag.eps']);
    
    figure;
    subplot(2,1,1); plot((1:length(rf))*dt, abs(rf),'r-*');           leg = legend('magnitude');     set(leg,'FontSize',17); set(gca,'FontSize',17); xlabel('Time (ms)','FontSize',17); ylabel('|RF| (Gauss)','FontSize',17);       title(['RF pulse designed by ', rfname],'FontSize',20); axis tight; set(gca,'xlim',[-2*dt, (length(rf)+3)*dt]); 
    subplot(2,1,2); plot((1:length(rf))*dt, angle(rf)*180/pi,'b-*');  leg = legend('phase');         set(leg,'FontSize',17); set(gca,'FontSize',17); xlabel('Time (ms)','FontSize',17); ylabel('\angle RF (degree)','FontSize',17); title(['RF pulse designed by ', rfname],'FontSize',20); axis tight; set(gca,'xlim',[-2*dt, (length(rf)+3)*dt]); set(gca,'ylim',[-180, 180]);
    set(gcf, 'Position', [50, 100, 700, 600], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_RF_mag_pha.eps']);
 
    % plot mz
    figure; 
    plot(df,mz,'b-','linewidth',1.8); xlabel('Frequency (Hz)','FontSize',18);  axis tight;  set(gca,'FontSize',18); ylabel('Mz','FontSize',18);  set(gca,'ylim',[ min(mz)-0.01,  1+0.01]);
    if ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
        hold on; plot_spec(f,a,d); hold off;
    end
    set(gcf, 'Position', [100, 100, 800, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_Mz.eps']);
    
    % plot mz in log scale
    figure;
    subplot(1,2,1); semilogy(df,abs(mz),'b-');   xlabel('Frequency (Hz)','FontSize',17); set(gca,'FontSize',17); axis tight; ylabel('|Mz|','FontSize',17); 
    subplot(1,2,2); semilogy(df,1-mz,'b-');      xlabel('Frequency (Hz)','FontSize',17); set(gca,'FontSize',17); axis tight; ylabel('1-Mz','FontSize',17); 
    if ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
        subplot(1,2,1); hold on; plot_spec_log(f,abs(a),d); hold off;
        subplot(1,2,2); hold on; plot_spec_log(f,1-a,d); hold off;
    end
    set(gcf, 'Position', [150, 100, 1000, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_Mz_log.eps']);
       
    % plot mz of each band
    if ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
        m_band = length(d);
        figure;
        for i=1:m_band
            idx_i = find( (df > (f(2*i-1)-30)) & (df < (f(2*i)+30)) );
            subplot(1,m_band,i); plot(df(idx_i),mz(idx_i),'b-','linewidth',1.5); xlabel('Frequency (Hz)','FontSize',18);  ylabel('Mz','FontSize',18); set(gca,'FontSize',18); % title(['Mz band ', num2str(i)],'FontSize',19); 
            hold on; plot_spec(f,a,d); hold off; set(gca,'ylim',[ a(2*i-1)-d(i)*2, a(2*i-1)+d(i)*2 ]); set(gca,'xlim',[min(df(idx_i)) max(df(idx_i))]);
        end
        set(gcf, 'Position', [1, 100, min([300*length(a) 1300]), 350], 'PaperPositionMode', 'auto');
        hgexport(gcf, [name_save, '_Mz_each_band.eps']);
    end
    
    % plot mxy magnitude and phase
    figure;
    subplot(1,2,1); plot(df,abs(mxy),'b-');  xlabel('Frequency (Hz)','FontSize',17);  set(gca,'FontSize',17); ylabel('|Mxy|','FontSize',17); axis tight; set(gca,'ylim',[0 max(abs(mxy))+0.02]);
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
        hold on; plot_spec(f,a,d); hold off;
    end
    subplot(1,2,2); plot(df,unwrap(angle(mxy)),'b-'); xlabel('Frequency (Hz)','FontSize',17);  set(gca,'FontSize',17); ylabel('\angle Mxy (radians)','FontSize',17); axis tight; 
    set(gcf, 'Position', [250, 100, 1000, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_Mxy.eps']);
    
    % plot mxy in log scale
    figure;
    subplot(1,2,1); semilogy(df,abs(mxy),'b-');     xlabel('Frequency (Hz)','FontSize',17); set(gca,'FontSize',17); axis tight; ylabel('|Mxy|','FontSize',17); 
    subplot(1,2,2); semilogy(df,1-abs(mxy),'b-');   xlabel('Frequency (Hz)','FontSize',17); set(gca,'FontSize',17); axis tight; ylabel('1-|Mxy|','FontSize',17); 
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
        subplot(1,2,1); hold on; plot_spec_log(f,a,d); hold off;
        subplot(1,2,2); hold on; plot_spec_log(f,1-a,d); hold off;
    end
    set(gcf, 'Position', [300, 100, 1000, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_Mxy_mag_log.eps']);

    % plot mxy of each band
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
        m_band = length(d);
        figure;
        for i=1:m_band
            idx_i = find( (df > (f(2*i-1)-30)) & (df < (f(2*i)+30)) );
            subplot(1,m_band,i); plot(df(idx_i), abs(mxy(idx_i)), 'b-'); xlabel('Frequency (Hz)','FontSize',17);  axis tight;   set(gca,'FontSize',17); title(['|Mxy| band ', num2str(i)],'FontSize',19); 
            hold on; plot_spec(f,a,d); hold off;  set(gca,'xlim',[min(df(idx_i)), max(df(idx_i))]);  set(gca,'ylim',[ a(2*i-1)-d(i)*2, a(2*i-1)+d(i)*2 ]);
        end
        set(gcf, 'Position', [1, 100, min([300*length(a) 1300]), 350], 'PaperPositionMode', 'auto');
        hgexport(gcf, [name_save, '_Mxy_mag_each_band.eps']);
    end
        
      
% =========================================================================
% bloch simulation if input two RF pulses
% =========================================================================
elseif ((iscell(rf)) && (length(rf) == 2))   
    
    rf1 = rf{1};  rf1 = rf1(:); 
    rf2 = rf{2};  rf2 = rf2(:);
    dt1 = dt(1);  dt2 = dt(2);
    name_save = [rfname{1},'_',rfname{2}];
    
    % plot rf pulse 
    range_y_rf = [min([ min(real(rf1)), min(imag(rf1)), min(real(rf2)), min(imag(rf2)) ]) ,   max([ max(real(rf1)), max(imag(rf1)), max(real(rf2)), max(imag(rf2)) ]) ];
    figure;
    subplot(1,2,1); plot((1:length(rf1))*dt1, real(rf1),'r.-',(1:length(rf1))*dt1, imag(rf1),'b.-','markersize',15); leg = legend('real','imag'); set(leg,'FontSize',17,'position',[0.32  0.79 0.1 0.1]); set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, (length(rf1)+1)*dt1]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{1}],'FontSize',19); 
    subplot(1,2,2); plot((1:length(rf2))*dt2, real(rf2),'r.-',(1:length(rf2))*dt2, imag(rf2),'b.-','markersize',13); leg = legend('real','imag'); set(leg,'FontSize',17,'position',[0.76 0.79 0.1 0.1]);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, (length(rf2)+1)*dt2]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{2}],'FontSize',19); 
    
    set(gcf, 'Position', [50, 100, 1000, 400], 'PaperPositionMode', 'auto');
    ax = axes('position',[0,0,1,1],'visible','off');  
    text(0.09,0.95, '(A)','FontSize',23);
    text(0.53,0.95, '(B)','FontSize',23);
    hgexport(gcf, [name_save, '_RF.eps']);
    
    % plot rf pulse amplitude
    figure;
    plot((1:length(rf1))*dt1, abs(rf1),'r.-', (1:length(rf2))*dt2, abs(rf2),'b.-', 'markersize',15); leg = legend(rfname{1},rfname{2});  set(leg,'FontSize',17, 'position',[0.75 0.77 0.1 0.1]);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, max([ (length(rf1)+1)*dt1,  (length(rf2)+1)*dt2 ]) ]);  set(gca,'ylim',[0 max([max(abs(rf1)) max(abs(rf2))])*1.01]);
    set(gcf, 'Position', [100, 100, 700, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_RF_amp.eps']);
    
    % run simulation 
    G1 = zeros(size(rf1,1),1);
    [mx1,my1,mz1] = bloch_fun_handle(rf1,G1,dt1*1e-3,1e3,1e3,df,dp,0,mx0,my0,mz0);
    mxy1 = mx1 + sqrt(-1)*my1;
    
    G2 = zeros(size(rf2,1),1);
    [mx2,my2,mz2] = bloch_fun_handle(rf2,G2,dt2*1e-3,1e3,1e3,df,dp,0,mx0,my0,mz0);
    mxy2 = mx2 + sqrt(-1)*my2;    

    % plot mz
    figure; 
    plot(df,mz1,'b-', df,mz2,'r-'); leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18);  axis tight;  set(gca,'FontSize',18); ylabel('Mz','FontSize',18);  set(gca,'ylim',[ min([min(mz1),min(mz2)]),  1+0.01] );
    if  ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
        hold on; plot_spec(f,a,d); hold off;
    end
    set(gcf, 'Position', [150, 100, 600, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_Mz.eps']);
    
    % plot mz in log scale
    figure;
    subplot(1,2,1); semilogy(df,abs(mz1),'b-', df,abs(mz2),'r-');  leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('|Mz|','FontSize',18); 
    subplot(1,2,2); semilogy(df,1-mz1,'b-', df,1-mz2,'r-');        leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('1-Mz','FontSize',18); 
    if ( (show_spec == 1) && (strcmp(ptype,'sat') || strcmp(ptype,'inv')) )
        subplot(1,2,1); hold on; plot_spec_log(f,a,d); hold off;
        subplot(1,2,2); hold on; plot_spec_log(f,1-a,d); hold off;
    end
    set(gcf, 'Position', [200, 100, 1000, 400], 'PaperPositionMode', 'auto'); 
    hgexport(gcf, [name_save, '_Mz_log.eps']);
    
    % plot mxy
    figure;
    subplot(1,2,1); plot(df,abs(mxy1),'b-', df,abs(mxy2),'r-','linewidth',1.5); leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17, 'position',[0.23 0.8 0.1 0.1]);  xlabel('Frequency (Hz)','FontSize',18);  axis tight;  set(gca,'FontSize',18); ylabel('|Mxy|','FontSize',18); set(gca,'ylim',[0 1]);
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
        hold on; plot_spec(f,a,d); hold off;
    end
    subplot(1,2,2); plot(df,myUnwrap(mxy1),'b-',  df,myUnwrap(mxy2),'r-','linewidth',1.5); leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17, 'position',[0.67 0.8 0.1 0.1]);  xlabel('Frequency (Hz)','FontSize',18);  axis tight;  set(gca,'FontSize',18); ylabel('\angleMxy (radians)','FontSize',18);
    set(gcf, 'Position', [250, 100, 1000, 400], 'PaperPositionMode', 'auto');
    
    ax = axes('position',[0,0,1,1],'visible','off');  
    text(0.07,0.95, '(A)','FontSize',23);
    text(0.51,0.95, '(B)','FontSize',23);
    hgexport(gcf, [name_save, '_Mxy.eps']);
    
%     % plot mxy in log scale
%     figure;
%     subplot(1,2,1); semilogy(df,abs(mxy1),'b-',   df,abs(mxy2),'r-');    leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',17); set(gca,'FontSize',17); axis tight; ylabel('|Mxy|','FontSize',17); 
%     subplot(1,2,2); semilogy(df,1-abs(mxy1),'b-', df,1-abs(mxy2),'r-');  leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',17); set(gca,'FontSize',17); axis tight; ylabel('1-|Mxy|','FontSize',17); 
%     if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
%         subplot(1,2,1); hold on; plot_spec_log(f,a,d); hold off;
%         subplot(1,2,2); hold on; plot_spec_log(f,1-a,d); hold off;
%     end
%     set(gcf, 'Position', [300, 100, 1000, 400], 'PaperPositionMode', 'auto');
%     hgexport(gcf, [name_save, '_Mxy_log.eps']);
    
    % plot mxy in log scale (only show stopband, without 1-|mxy|)
    figure;
    semilogy(df,abs(mxy1),'b-',   df,abs(mxy2),'r-', 'linewidth',1.5);    leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17, 'position',[0.76 0.15 0.1 0.1]);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('|Mxy|','FontSize',18); 
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) && isempty(name_cell) )
        hold on; plot_spec_log(f,a,d); hold off;
    end
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) && (~isempty(name_cell)) )
        hold on; plot_spec_log(f,a,d, name_cell); hold off;
    end
    set(gcf, 'Position', [300, 100, 700, 400], 'PaperPositionMode', 'auto');
    hgexport(gcf, [name_save, '_Mxy_log.eps']);

    
    % plot mxy of each band
    if ( (show_spec == 1) && (strcmp(ptype,'st') || strcmp(ptype,'ex')) )
        m_band = length(d);
        figure;
        for i=1:m_band
            idx_i = find( (df > (f(2*i-1)-30)) & (df < (f(2*i)+30)) );
            subplot(1,m_band,i);  plot(df(idx_i), abs(mxy1(idx_i)), 'b-', df(idx_i), abs(mxy2(idx_i)), 'r-'); % leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',16); 
            xlabel('Frequency (Hz)','FontSize',17);  axis tight;  set(gca,'FontSize',17); %title(['|Mxy| band ', num2str(i)],'FontSize',17); 
            hold on; plot_spec(f,a,d); hold off;  set(gca,'ylim',[ a(2*i-1)-d(i)*2, a(2*i-1)+d(i)*2 ]); set(gca,'xlim',[min(df(idx_i)) max(df(idx_i))]);
        end
        set(gcf, 'Position', [50, 100, min([280*length(a) 1400]), 350], 'PaperPositionMode', 'auto'); 
        hgexport(gcf, [name_save, '_Mxy_each_band.eps']);
    end
    
    
    % plot all in one: RF pulse (real/imag), magnitude profile (log), phase profile
    if (strcmp(ptype,'st') || strcmp(ptype,'ex')) % display Mxy
        figure;
        range_x_rf = max([(length(rf1)+1)*dt1, (length(rf2)+1)*dt2]);
        subplot(2,2,1); set(gca, 'OuterPosition', [0.02, 0.51,  0.41, 0.47]); plot((1:length(rf1))*dt1, real(rf1),'r.-',(1:length(rf1))*dt1, imag(rf1),'b.-','markersize',14); leg = legend('real','imag'); set(leg,'FontSize',17);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, range_x_rf]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{1}],'FontSize',19); 
        subplot(2,2,3); set(gca, 'OuterPosition', [0.02, 0.02,  0.41, 0.47]); plot((1:length(rf2))*dt2, real(rf2),'r.-',(1:length(rf2))*dt2, imag(rf2),'b.-','markersize',14); leg = legend('real','imag'); set(leg,'FontSize',17);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, range_x_rf]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{2}],'FontSize',19); 
        subplot(2,2,2); set(gca, 'OuterPosition', [0.47, 0.51,  0.53, 0.47]); semilogy(df,abs(mxy1),'b-',   df,abs(mxy2),'r-', 'linewidth',1.5);    leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17, 'position',[0.88 0.61 0.01 0.01]);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('|Mxy|','FontSize',18); 
        if ( (show_spec == 1) && isempty(name_cell) ),    hold on; plot_spec_log(f,a,d); hold off;  end;
        if ( (show_spec == 1) && (~isempty(name_cell)) ), hold on; plot_spec_log(f,a,d, name_cell); hold off; end;
        subplot(2,2,4); set(gca, 'OuterPosition', [0.47, 0.02,  0.53, 0.47]); plot(df,myUnwrap(mxy1),'b-',  df,myUnwrap(mxy2),'r-','linewidth',1.5); leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18);  axis tight;  set(gca,'FontSize',18); ylabel('\angleMxy (radians)','FontSize',18);

        set(gcf, 'Position', [50, 100, 1100, 700], 'PaperPositionMode', 'auto');
        ax = axes('position',[0,0,1,1],'visible','off');  
        text(0.04,0.96, '(A)','FontSize',21);
        text(0.49,0.96, '(B)','FontSize',21);
        text(0.04,0.49, '(C)','FontSize',21);
        text(0.49,0.49, '(D)','FontSize',21);
        hgexport(gcf, [name_save, '_all.eps']);
        
       
        % plot FA, phase difference between the two pulses, assuming an excitation pulse
        [FA1, pha1] = myRF2FA(rf1, dt1, df*1e-3, nucleus, 0);
        [FA2, pha2] = myRF2FA(rf2, dt2, df*1e-3, nucleus, 0);
        pha_shift = pha1-pha2; % in degree
        pha_shift = round( mod(pha_shift+179.999,360)-179.999 ) ; % in degree, rewrite phase 

        figure;
        range_x_rf = max([(length(rf1)+1)*dt1, (length(rf2)+1)*dt2]);
        subplot(2,2,1); set(gca, 'OuterPosition', [0.02, 0.51,  0.41, 0.47]); plot((1:length(rf1))*dt1, real(rf1),'r.-',(1:length(rf1))*dt1, imag(rf1),'b.-','markersize',14); leg = legend('real','imag'); set(leg,'FontSize',17);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, range_x_rf]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{1}],'FontSize',19); 
        subplot(2,2,3); set(gca, 'OuterPosition', [0.02, 0.02,  0.41, 0.47]); plot((1:length(rf2))*dt2, real(rf2),'r.-',(1:length(rf2))*dt2, imag(rf2),'b.-','markersize',14); leg = legend('real','imag'); set(leg,'FontSize',17);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, range_x_rf]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{2}],'FontSize',19); 
        subplot(2,2,2); set(gca, 'OuterPosition', [0.47, 0.51,  0.53, 0.47]); semilogy(df,abs(mxy1),'b--',   df,abs(mxy2),'r:', 'linewidth',1.5);    leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17, 'position',[0.88 0.61 0.01 0.01]);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('Flip Angle (degree)','FontSize',18); 
        if ( (show_spec == 1) && isempty(name_cell) ),    hold on; plot_spec_log(f,a,d); hold off;  end;
        if ( (show_spec == 1) && (~isempty(name_cell)) ), hold on; plot_spec_log(f,a,d, name_cell); hold off; end;
        subplot(2,2,4); set(gca, 'OuterPosition', [0.47, 0.02,  0.53, 0.47]); plot(df,pha_shift,'b-', df(1),pha_shift(1),'r.'); leg = legend('total','ROI'); set(leg,'FontSize',17);  xlabel('Frequency (Hz)','FontSize',18);  axis tight;  set(gca,'FontSize',18); ylabel('RF phase shift (degree)','FontSize',18);

        set(gcf, 'Position', [5, 100, 1100, 700], 'PaperPositionMode', 'auto');
        ax = axes('position',[0,0,1,1],'visible','off');  
        text(0.04,0.96, '(A)','FontSize',21);
        text(0.49,0.96, '(B)','FontSize',21);
        text(0.04,0.49, '(C)','FontSize',21);
        text(0.49,0.49, '(D)','FontSize',21);
        hgexport(gcf, [name_save, '_all_phase_shift.eps']);
        
        
    elseif (strcmp(ptype,'sat') || strcmp(ptype,'inv'))  % display Mz
        figure;
        range_x_rf = max([(length(rf1)+1)*dt1, (length(rf2)+1)*dt2]);
        subplot(2,2,1); set(gca, 'OuterPosition', [0.02, 0.51,  0.41, 0.47]); plot((1:length(rf1))*dt1, real(rf1),'r-',(1:length(rf1))*dt1, imag(rf1),'b-','linewidth',2); leg = legend('real','imag'); set(leg,'FontSize',17);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, range_x_rf]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{1}],'FontSize',19); 
        subplot(2,2,3); set(gca, 'OuterPosition', [0.02, 0.02,  0.41, 0.47]); plot((1:length(rf2))*dt2, real(rf2),'r-',(1:length(rf2))*dt2, imag(rf2),'b-','linewidth',2); leg = legend('real','imag'); set(leg,'FontSize',17);  set(gca,'FontSize',18); xlabel('Time (ms)','FontSize',18); ylabel('RF (Gauss)','FontSize',18); axis tight; set(gca,'xlim',[0, range_x_rf]); set(gca,'ylim',range_y_rf); %title(['RF pulse designed by ', rfname{2}],'FontSize',19); 
        subplot(2,2,2); semilogy(df,abs(mz1),'b-', df,abs(mz2),'r-');  leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17,'Position',[0.63 0.62 0.01 0.01]);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('|Mz|','FontSize',18); 
        subplot(2,2,4); semilogy(df,1-mz1,'b-', df,1-mz2,'r-');        leg = legend(rfname{1},rfname{2}); set(leg,'FontSize',17,'Position',[0.63 0.145 0.01 0.01]);  xlabel('Frequency (Hz)','FontSize',18); set(gca,'FontSize',18); axis tight; ylabel('1-Mz','FontSize',18); 
        if ( (show_spec == 1) && isempty(name_cell) )
            subplot(2,2,2); hold on; plot_spec_log(f,a,d); hold off;
            subplot(2,2,4); hold on; plot_spec_log(f,1-a,d); hold off;
        elseif ( (show_spec == 1) && (~isempty(name_cell)) )
            subplot(2,2,2); hold on; plot_spec_log(f,a,d); hold off;
            subplot(2,2,4); hold on; plot_spec_log(f,1-a,d,name_cell); hold off;
        end
        subplot(2,2,2); set(gca,'ylim', [1e-3 1]);
        subplot(2,2,4); set(gca,'ylim', [min([min(1-mz1), min(1-mz2)]) 1]);
        
        set(gcf, 'Position', [50, 100, 1100, 700], 'PaperPositionMode', 'auto');
        ax = axes('position',[0,0,1,1],'visible','off');  
        text(0.04,0.96, '(A)','FontSize',21);
        text(0.49,0.96, '(B)','FontSize',21);
        text(0.04,0.49, '(C)','FontSize',21);
        text(0.49,0.49, '(D)','FontSize',21);
        hgexport(gcf, [name_save, '_all.eps']);
        
    end

end