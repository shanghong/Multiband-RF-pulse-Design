function [rf_pulse,b, rf_spec, b_spec] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg, min_order, min_tran, shift_f, name_cell)
%
%  [rf_pulse,b, rf_spec, b_spec] = dzrf_mb(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, ftype, nucleus, flip_zero, downsampling, Peak, dbg, min_order, min_tran, shift_f, name_cell)
%
%  Design a RF pulse using the SLR algorithm, based on John Pauly's rf_tools 
%  package.  Function conforms with dzrf.m (rf_tools)
%  New feature:
%      (1) multiband profile
%      (2) multiple methods for FIR filter designs with arbitrary phase,
%          quadratic phase. minimize duration, transitino, etc
%      (3) generalized arbitrary flip-angle
%      (4) upsampling after filter design to reduce computation time
%      (5) shift frequency spec to center, design pulse, and shift back, to
%          reduce SLR error for large FA at high frequency
%  cvx is needed for FIR filter design
%
%  Inputs:
%    n -- order of filter, number of hard pulse
%    dt -- in ms, sampling interval of RF pulse/filter
%    mb_cf --    1*m in kHz (cell array), center frequency for each band,
%                    a point (legnth(mb_cf{i}) = 1) or a range (legnth(mb_cf{i}) = 2)
%    mb_range -- 1*m in kHz, width for each band, [-mb_range/2, mb_range/2] 
%    mb_FA --    1*m in degree, flip-angle for each band
%    mb_ripple-- 1*m, absolute ripple of magnetization for each band
%    ptype -- pulse type.  Options are:
%       st  -- small tip angle                    
%       ex  -- excitation pulse
%       se  -- spin-echo pulse
%       sat -- saturation pulse                   [default]
%       inv -- inversion pulse
%    ftype -- FIR filter design method.  Options are:
%       ms  -- Hamming windowed sinc
%       ap_cvx -- arbitrary-phase filter by convex optimization, minimum energy [default]
%       ap_minstopripple_cvx -- ap_cvx and minimum stopband ripple
%       ap_minorder_cvx -- ap_cvx and minimum order/length
%       ap_mintran_cvx -- ap_cvx and minimum transition width
%       ap_mintran_minorder_cvx -- minimize transition and order
%       lp_minorder -- linear-phase filter with minimum order
%       qp_cvx -- quadratic-phase fitler by convex optimization
%    nucleus -- options are 'H-1', 'C-13', Carbon as default            [optional]
%    flip_zero -- if 1, use flip-zero method to minimize peak amplitude [optional]
%    downsampling -- downsampling - filer design - upsampling           [optional]
%    Peak -- constraint on peak amplitude, especially at end            [optional]
%    dbg -- flag to turn on debugging statements/plots                  [optional]
%    min_order -- [0 1] to adjust 'ap_minorder_cvx' filter design       [optional]
%    min_tran -- [0 1] to adjust 'ap_mintran_cvx' filter design         [optional]
%    shift_f -- shift frequency to center - design pulse - shift back   [optional]
%    name_cell -- name of each band                                     [optional]
%
%  Outputs:
%    rf_pulse -- in Gauss, rf pulse 
%    b -- beta polynomial used for rf pulse design
%    rf_spec -- a structure, spec for RF pulse design
%    b_spec -- a structure, spec for beta polynomial design  
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Nov 2014

% check input
if  (nargin < 6),   error('not enough input'); end;
if  (nargin <= 6),  ptype = 'sat';  end;
if  (nargin <= 7),  ftype = 'ap_cvx';  end;
if ((nargin <= 8)  || isempty(nucleus)),  nucleus = 'C-13';  end; 
if ((nargin <= 9)  || isempty(flip_zero)),  flip_zero = 0;  end;
if ((nargin <= 10) || isempty(downsampling)),  downsampling = 1;  end;
if ((nargin <= 11) || isempty(Peak)),  Peak = 1e-3;  end;
if  (nargin <= 12), dbg = 0;  end;
if ((nargin <= 13) || isempty(min_order)), min_order = 0.9;  end
if ((nargin <= 14) || isempty(min_tran)),  min_tran = 0.85;  end
if ((nargin <= 15) || isempty(shift_f)),  shift_f = 0;  end
if ((nargin <= 16) || isempty(name_cell)),  name_cell = [];  end
if  (nargin > 17),  error('too many input'); end;

rf_pulse = [];
b = [];
rf_spec = [];
b_spec = [];

% internal parameters
m_band = length(mb_FA); % number of band
switch nucleus
    case 'H-1'
        gamma = 4.2576; % in kHz/G, gamma/2pi for H-1
    case 'C-13'
        gamma = 1.0705; % in kHz/G, gamma/2pi for C-13
    otherwise
        error('No such option for nucleus. Options are H-1 and C-13 \n');
end

% downsampling
n = n/downsampling;     % new number of samples after downsampling
dt = dt*downsampling;   % in ms, new sampling interval
fs = 1/dt;              % in kHz,new sampling rate
if abs(n-round(n)) > 1e-10  % n is not integer
    error('n/downsampling is not an integer');
end

% spec of frequency range for each band within [-1 1]
f = rf_bandedge(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, dbg-2);
       
% spec of amplitude range for each band
a = zeros(1,2*m_band);
d = zeros(1,m_band);
for i=1:m_band
    [range_B, ~] = rf_ripple_GFA(mb_FA(i), mb_ripple(i), ptype, 0, dbg-2);
    a(2*i-1) = (range_B(2) + range_B(1))/2;
    a(2*i)  =  (range_B(2) + range_B(1))/2;
    d(i)   =   (range_B(2) - range_B(1))/2;
end

% output spec for both b and rf for debug
a_M = zeros(1,2*m_band);
d_M = zeros(1,m_band);
for i=1:m_band
    range_M = rf_Mrange_desired(mb_FA(i), mb_ripple(i), ptype);
    a_M(2*i-1) = (range_M(2) + range_M(1))/2;
    a_M(2*i)  =  (range_M(2) + range_M(1))/2;
    d_M(i)   =   (range_M(2) - range_M(1))/2;
end
b_spec = struct('f',{f/downsampling},'a',{a},'d',{d});
rf_spec = struct('f',{f/downsampling},'a',{a_M},'d',{d_M});

% display spec for debug 
if dbg >= 3
    if isempty(name_cell)
        figure; plot_spec(f*(fs/2), a_M, d_M); xlabel('Frequency (kHz)','FontSize',18);      ylabel('Multiband Spec of |Mxy|','FontSize',18);    set(gca,'FontSize',18);  axis tight;  set(gca,'xlim',[min(f*(fs/2)) - 0.1, max(f*(fs/2)) + 0.1]); set(gca,'ylim',[min(a_M)-max(d_M)-0.02 max(a_M)+max(d_M)+0.02]);  set(gcf, 'Position', [100,100,500,350], 'PaperPositionMode', 'auto'); hgexport(gcf, 'multiband_spec_M.eps');
        figure; plot_spec(f, a, d);            xlabel('Normalized Frequency','FontSize',18); ylabel('Multiband Spec of |B_N(z)|','FontSize',18); set(gca,'FontSize',18);  axis tight;  set(gca,'xlim',[min(f)-0.02, max(f)+0.02]);                 set(gca,'ylim',[min(a)-max(d)-0.02 max(a)+max(d)+0.02]);          set(gcf, 'Position', [200,100,500,350], 'PaperPositionMode', 'auto'); hgexport(gcf, 'multiband_spec_B.eps');
    else
        figure; plot_spec(f*(fs/2), a_M, d_M, name_cell); xlabel('Frequency (kHz)','FontSize',18);      ylabel('Multiband Spec of |Mxy|','FontSize',18);    set(gca,'FontSize',18);  axis tight;  set(gca,'xlim',[min(f*(fs/2)) - 0.1, max(f*(fs/2)) + 0.1]); set(gca,'ylim',[min(a_M)-max(d_M)-0.02 max(a_M)+max(d_M)+0.02]); set(gcf, 'Position', [100,100,500,350], 'PaperPositionMode', 'auto'); hgexport(gcf, 'multiband_spec_M.eps');
        figure; plot_spec(f, a, d, name_cell);            xlabel('Normalized Frequency','FontSize',18); ylabel('Multiband Spec of |B_N(z)|','FontSize',18); set(gca,'FontSize',18);  axis tight;  set(gca,'xlim',[min(f)-0.02, max(f)+0.02]);                 set(gca,'ylim',[min(a)-max(d)-0.02 max(a)+max(d)+0.02]);         set(gcf, 'Position', [200,100,500,350], 'PaperPositionMode', 'auto'); hgexport(gcf, 'multiband_spec_B.eps');
    end
end
    
% shift frequency if necessary
switch shift_f
    case 0  % not shift
        shift_f_back = 0;  % in kHz
        shift_f_forward = 0; 
        
    case 1  % shfit such that center of high FA bands is at f=0
        idx_FA = find( mb_FA>60 );
        idx_FA = [idx_FA(1) idx_FA(end)];
        % shift_f_forward = mean( [f(2*idx_FA-1), f(2*idx_FA)] );
        shift_f_forward = mean( [f(2*idx_FA(1)-1), f(2*idx_FA(2))] );
        f = f - shift_f_forward;
        shift_f_back = shift_f_forward * 0.5*(1/dt);  % in kHz
        
    case 2  % shift such that highest FA band is at f=0
        [~, idx_FA] = max(mb_FA);
        shift_f_forward = mean( [f(2*idx_FA-1), f(2*idx_FA)] );
        f = f - shift_f_forward;
        shift_f_back = shift_f_forward * 0.5*(1/dt);  % in kHz
        
    otherwise
        error(['shift_f = ', num2str(shift_f), ' is not an option. Options are 0,1,2']);
end



status = [];
% different FIR filter design
switch ftype
    case 'ms'          % from John, windowed sinc
        b = msinc(n, TBW/4);

    case 'ap_cvx'
        [b, status] = fir_ap_cvx(n, f, a, d, 1, Peak, dbg);
    
    case 'ap_minstopripple_cvx'
        [b, status] = fir_ap_cvx(n, f, a, d, 1e4, Peak, dbg);
        
    case 'ap_minorder_cvx' 
        if (min_order <= 1)  % min_order in [0 1], which has the original meaning
            [b, status] = fir_ap(n, f, a, d, Peak, min_order, 0, 0, dbg);
        else                 % design with fixed n (=min_order)
            n_new = min_order;
            [b, status] = fir_ap_cvx(n_new, f, a, d, 0.1, Peak);
            fprintf('Finally solved with n = %d \n', n_new);
        end
    
        if dbg >= 1
            fprintf('Minimize order during filter design\n');
            fprintf('Previous duration = %3.3f ms, new duration = %3.3f ms, reduced by %0.4f \n', n*dt, length(b)*dt, (n-length(b))/n );
        end
        
    case 'ap_mintran_cvx'
        [b, status, ~, f_new] = fir_ap(n, f, a, d, Peak, 0, min_tran, 0, dbg);
        f_bw_old = (f(2:2:end) - f(1:2:(end-1))) * (fs*1e3/2);         % in Hz, bandwidth of each band
        f_bw_new = (f_new(2:2:end) - f_new(1:2:(end-1))) * (fs*1e3/2); % in Hz, new bandwidth of each band
        
        if dbg >= 1
            fprintf('\n'); fprintf('Minimize transition width during filter design\n');
            for ii=1:m_band
                fprintf('For Band %d, original BW = %5.2f Hz, new BW = %5.2f Hz, increased by %5.2f Hz \n', ii, f_bw_old(ii), f_bw_new(ii), f_bw_new(ii)-f_bw_old(ii));
            end
        end

        b_spec.f = ( f_new + shift_f_forward ) /downsampling;
        rf_spec.f = ( f_new + shift_f_forward )/downsampling;

    case 'ap_mintran_minorder_cvx'
        [b, status] = fir_ap(n, f, a, d, 0.4, 0.3, 0, dbg);
        
    case 'lp_minorder'
        addpath ./ss
        [b, status] = fir_min_order_linprog(n, f, a, d, [], dbg);
        %[b, status] = fir_min_order(n, f, a, d, [], [], dbg);
        
    case 'qp_cvx'
        k = 120;      % scaling factor for quadratic phase
        obj = 1e6;   % optimization parameter for minimizing peak-amplitude
        [b, status] = fir_qp_cvx(n, f, a, d, k, obj, dbg);
end

if strcmp(status, 'Failed')
    disp('Filter design failed.');
    return;
else
    b = b(end:-1:1);
end

% flip-zero method to reduce peak amplitude with arbitrary phase
if flip_zero
    b = fir_flip_zero(b, dbg);
end       

% upsampling to the original sampling rate
if downsampling >= 2
    b = fir_upsample(b, dt, dt/downsampling, dbg);
end
dt = dt/downsampling;
b = transpose(b(:));

% inverse SLR transform
if strcmp(ptype,'st'),
   rf = b;
else
   a = b2a(b); 
   rf = ab2rf(a,b);
end;

% scaled to Gauss
rf_pulse = rfscaleg(rf,dt*length(rf),gamma);

% display the process of SLR, only for excitation pulse (Mxy)
if ( (dbg >= 3) && strcmp(ptype,'ex') && (downsampling == 1) )
    
    % reverse b and a only for display purpose
    b_reverse = b(end:-1:1);
    a_reverse = b2a(b_reverse);
    
    wdisp = linspace(-1,1,1e4);
    Anz = fftshift(fft(a_reverse, length(wdisp)));
    Bnz = fftshift(fft(b_reverse, length(wdisp)));
    Mxy = 2 * conj(Anz) .* Bnz;
    figure;
    subplot(2,2,1); plot(1:length(b),real(b),'b-', 1:length(b),imag(b),'b--', 1:length(a),real(a)/10,'r-', 1:length(a),imag(a),'r--','linewidth',1.5);     
    leg = legend('\beta Real','\beta Imag', '\alpha Real/10','\alpha Imag'); set(leg,'FontSize',17); 
    set(gca,'FontSize',18); ylabel('\beta & \alpha polynomial','FontSize',18);  axis tight;
    
    subplot(2,2,3); plot(1:length(rf_pulse),real(rf_pulse),'b-', 1:length(rf_pulse),imag(rf_pulse),'b--','linewidth',2); leg = legend('real','imag'); set(leg,'FontSize',17); set(gca,'FontSize',18); ylabel('RF (Gauss)','FontSize',18);  axis tight;

    subplot(2,2,2); plot(wdisp,abs(Bnz),'b-', wdisp,abs(Anz),'r-','linewidth',1.5); leg = legend('|B_N(z)|','|A_N(z)|'); set(leg,'FontSize',17); 
    hold on; plot_spec(b_spec.f, b_spec.a, b_spec.d, name_cell);  hold off;  
    xlabel('Normalized Frequency','FontSize',18);  set(gca,'FontSize',18); axis tight;           
    axis tight;  set(gca,'xlim',[min(b_spec.f)-0.02, max(b_spec.f)+0.02]);  set(gca,'ylim',[min(b_spec.a)-max(b_spec.d)-0.02 1+0.02]);     
    
    subplot(2,2,4); plot(wdisp*(fs/2), abs(Mxy),'b-','linewidth',1.5); 
    hold on; plot_spec(f*(fs/2), a_M, d_M, name_cell); hold off;
    xlabel('Frequency (kHz)','FontSize',18); ylabel('|Mxy|','FontSize',18); set(gca,'FontSize',18); axis tight; set(gca,'xlim',[min(f*(fs/2)) - 0.1, max(f*(fs/2)) + 0.1]); set(gca,'ylim',[min(a_M)-max(d_M)-0.02 max(a_M)+max(d_M)+0.02]); 
    
    set(gcf, 'Position', [150, 150, 750, 600], 'PaperPositionMode', 'auto');
    hgexport(gcf, 'Example_SLR_transform.eps');
    
end

% shift frequency [optional]
t_axis = [1:length(rf)]*dt; % in ms
rf_pulse = rf_pulse .* exp(sqrt(-1) * 2*pi*shift_f_back * t_axis);

