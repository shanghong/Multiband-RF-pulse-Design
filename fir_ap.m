function [h, status, n_op, f_op] = fir_ap(n, f, a, d, Peak, min_order, min_tran, min_peak, dbg) 
%
% Design arbitrary-phase filter that meets multiband frequency respones 
% manitude specification, using fir_ap_cvx.
%
% Additional option:
%     (1) minimize order filter by bisection search
%     (2) minimize transition bandwidth by bisection search
%     (3) minimize both above
%     (4) minimize peak amplitude by changing phase
%
% [h, status] = fir_ap(n, f, a, d, Peak, min_order, min_tran, min_peak, dbg) 
%
% Inputs: --- similar to cfirpm
%   n: initial number of taps
%   f: frequency bands (-1->1)
%   a: amplitude at band edges
%   d: ripple  in bands
%   Peak: constraint on peak amplitude, especially at ends
%   min_order: [0 1], 0 = original order, 1 = minimize order by bisection search, (0 1) somewhere in between
%   min_tran:  [0 1], 0 = original transition band, 1 = minimize transition bandwidth by bisection search, (0 1) somewhere in between
%              the above two are usually not 1 due to RF pulse error when pushing to the limit 
%   min_peak:  1 = minimize peak amplitude by flipping zeros and changing phase
%   dbg: flag to turn on debugging statements/plots
%
% Outputs: 
%   h: filter coefficients
%   status: 'Solved' or 'Failed';
%   n_op: new n
%   f_op: new f
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  June 2014

if (nargin < 4),  error('not enough input'); end;
if (nargin == 4), Peak = 1e-3;    end;
if (nargin <= 5), min_order = 0;  end;
if (nargin <= 6), min_tran = 0;   end;
if (nargin <= 7), min_peak = 0;   end;
if (nargin <= 8), dbg = 0;        end;
if (nargin > 9),  error('too many input'); end;

% internal parameter
df_thre = 0.0005;
lambda = 0.1;  % minimize total energy, instead of stopband ripple
n_op = n;
f_op = f;

% design with original parameters
[h1, status1] = fir_ap_cvx(n, f, a, d, lambda, Peak, dbg);
if strcmp(status1, 'Failed')
    error('original parameters are too tight');
end

if ( (min_tran == 0) && (min_order == 0) )
    h = h1;  
    status = status1;
    return;
end

%  minimize transition bandwidth by bisection search
if (min_tran > 0)  
    
    f_bw = f(2:2:end) - f(1:2:(end-1));  % bandwidth of each band
    if ( norm(diff(f_bw)) > 1e-10 )  % bandwidth is not the same
        disp('Warning: bandwidth is not the same for each band');
    end
    
    df_array = f(3:2:(end-1)) - f(2:2:(end-2)); % transition width, maybe different
    df_min = min(df_array);  % minimum transition width, which is the limit of adding bandwidth
    
    % add bandwidth with a constant value
    f_add_bot = 0;        % filter can be designed
    f_add_top = df_min/2; % filter cannot be designed
    fprintf('filter design Solved with adding BW = %0.4f \n', f_add_bot);
        
    while 1
        f_add_mid = (f_add_bot + f_add_top)/2;
        % expand each band with f_add_mid 
        f_new = f;
        f_new(1:2:(end-1)) = f_new(1:2:(end-1)) - f_add_mid;
        f_new(2:2:end) =  f_new(2:2:end) + f_add_mid;      
        [h0, status0] = fir_ap_cvx(n, f_new, a, d, lambda, Peak);
        
        if strcmp(status0, 'Failed')  % f_add_mid is too large, try something below it
            f_add_top = f_add_mid;
            fprintf('filter design Failed with adding BW = %0.4f \n', f_add_mid);
        else  % f_add_mid works, try something above it
            h = h0;  status = status0;
            f_add_bot = f_add_mid;             
            fprintf('filter design Solved with adding BW = %0.4f \n', f_add_mid);
        end
        
        if dbg >= 2
            figure(101);
            plot_spec(f, a, d,'b-*'); hold on; plot_spec(f_new, a, d,'r-'); hold off;
            xlabel('normalized frequency','FontSize',18); title('Expanded Band Edge','FontSize',20); set(gca,'FontSize',18); axis tight; set(gca,'ylim',[min(a)-0.3 max(a)+0.3]); set(gca,'xlim',[f(1)-0.1, f(end)+0.1]);
            set(gcf, 'Position', [500, 100, 800, 400], 'PaperPositionMode', 'auto'); drawnow;
        end
            
        if (f_add_top - f_add_bot) < df_thre                
            break;
        end
    end 
end


if min_tran  == 0
    fprintf('Finally solved with original frequency spec \n');
elseif ( (min_tran <= 1) && (min_tran > 0) )
    f_add_mid = 0*(1-min_tran) + f_add_bot*min_tran;
    f_new = f;
    f_new(1:2:(end-1)) = f_new(1:2:(end-1)) - f_add_mid;
    f_new(2:2:end) =  f_new(2:2:end) + f_add_mid;
    [h0, status0] = fir_ap_cvx(n, f_new, a, d, lambda, Peak);
    if  strcmp(status0, 'Failed')
        % use the f_add which is already tested
        f_add_mid = f_add_bot;
        f_new = f;
        f_new(1:2:(end-1)) = f_new(1:2:(end-1)) - f_add_mid;
        f_new(2:2:end) =  f_new(2:2:end) + f_add_mid;
    else
        h = h0; status = status0;
    end
    fprintf('Finally solved with adding BW = %0.4f \n', f_add_mid);
    
    if dbg >= 2
        figure(101);
        plot_spec(f, a, d,'b-*'); hold on; plot_spec(f_new, a, d,'r-'); hold off;
        xlabel('normalized frequency','FontSize',18); title('Expanded Band Edge','FontSize',20); set(gca,'FontSize',18); axis tight; set(gca,'ylim',[min(a)-0.3 max(a)+0.3]); set(gca,'xlim',[f(1)-0.1, f(end)+0.1]);
        set(gcf, 'Position', [500, 100, 800, 400], 'PaperPositionMode', 'auto'); drawnow;
    end
    
    f = f_new;
    f_op = f_new;
else
    error('invalid input of min_tran');
end


%  minimize order by bisection search
if (min_order > 0)  
    n_top = n;  % filter can be designed
    n_bot = 2;  % filter cannot be designed
    while 1
        n_mid = ceil((n_top + n_bot)/2);
        [h0, status0] = fir_ap_cvx(n_mid, f, a, d,lambda, Peak, dbg);
        if strcmp(status0, 'Failed')  % n_mid is still too low
            n_bot = n_mid;
            fprintf('filter design Failed with n = %d \n', n_mid);
        else  % n_mid works, try something below it
            h = h0;  status = status0;
            n_top = n_mid;
            fprintf('filter design Solved with n = %d \n', n_mid)
        end
 
        if (n_top - n_bot) == 1                
            break;
        end
    end 
end

if min_order == 1
    fprintf('Finally solved with n = %d \n', n_top);
    n_op = n_top;
elseif min_order  == 0
    fprintf('Finally solved with n = %d \n', n);
elseif ( (min_order < 1) && (min_order > 0) )
    n_new = ceil(n*(1-min_order) + n_top*min_order);
    [h, status] = fir_ap_cvx(n_new, f, a, d,lambda,Peak);
    fprintf('Finally solved with n = %d \n', n_new);
    n_op = n_new;
else
    error('invalid input of min_order');
end

if dbg >= 1
    % display for the min_order mode, compare the original and final filter
    wdisp = linspace(-1, 1, 1e4); 
    H = fftshift(fft(h, length(wdisp))); % final filter
    H1 = fftshift(fft(h1, length(wdisp))); % Initial filter
    
    figure; 
    subplot(1,2,1); plot(1:length(h1),real(h1),'b-',1:length(h1),imag(h1),'b--',   1:length(h),real(h),'r-',1:length(h),imag(h),'r--', 'linewidth',2); 
    leg = legend('Initial Real','Initial Imag', 'Optimal Real','Optimal Imag'); set(leg,'FontSize',17); 
    ylabel('Filter Coefficients','FontSize',18); set(gca,'FontSize',18); axis tight;
    
    subplot(1,2,2); plot(wdisp,abs(H1),'b-', wdisp,abs(H),'r-', 'linewidth',1.5); 
    leg = legend('Initial','Optimal'); set(leg,'FontSize',17);
    hold on; plot_spec(f, a, d); hold off; 
    ylabel('|H(e^j^w)|','FontSize',18); xlabel('Normalized Frequency','FontSize',18); set(gca,'FontSize',18);  set(gca,'xlim',[min(f)-0.02 max(f)+0.02]); set(gca,'ylim',[min(a)-min(d)-0.02 max(a)+max(d)+0.02]);
    set(gcf, 'Position', [100,100,800,350], 'PaperPositionMode', 'auto'); 
    hgexport(gcf, 'FIR_design_minOrder.eps');
    
end
    
% use flip-zero method to reduce peak amplitude
if min_peak
    if dbg >=1    
        fprintf('start to run fir_flip_zero\n');
        tic;
        h = fir_flip_zero(h, dbg);
        T_flip_zero = toc;
        fprintf('computation time %6.4f s\n', T_flip_zero);
    else
        h = fir_flip_zero(h, dbg);
    end
end

                
        
        
   