function [h, status] = fir_qp(n, f, a, d, min_order, min_tran, min_peak, dbg) 
%
% Design quadratic-phase filter using fir_qp_cvx.
%
% Additional option:
%     (1) bisection search for minimizing transition bandwidth
%     (2) bisection search for minimizing order
%
% [h, status] = fir_qp(n, f, a, d, min_order, min_tran, min_peak, dbg) 
%
% Inputs: --- similar to cfirpm
%   n: initial number of taps
%   f: frequency bands (-1->1)
%   a: amplitude at band edges
%   d: ripple  in bands
%   min_order: [0 1], 0 = original order, 1 = minimize order by bisection search, (0 1) somewhere in between
%   min_tran:  [0 1], 0 = original transition band, 1 = minimize transition bandwidth by bisection search, (0 1) somewhere in between
%              the above two are usually not 1 due to RF pulse error when pushing to the limit 
%   min_peak:  1 = minimize peak amplitude by flipping zeros and changing phase
%   dbg: flag to turn on debugging statements/plots
%
% Outputs: 
%   h: filter coefficients
%   status: 'Solved' or 'Failed'; 
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  June 2014

df_thre = 0.001;
lambda = 1e5;

if (nargin < 4),  error('not enough input'); end;

if (nargin == 4), min_order = 0;  end;

if (nargin <= 5), min_tran = 0;   end;

if (nargin <= 6), min_peak = 0;   end;

if (nargin <= 7), dbg = 0;        end;
    
if (nargin > 8), error('too many input'); end;


% design with original parameters
[h1, status1] = fir_ap_cvx(n, f, a, d,lambda);
if strcmp(status1, 'Failed')
    error('original parameters are too tight');
end

if (min_tran == 0) & (min_order == 0)
    h = h1;  status = status1;
end

% determine transition bandwidth
if (min_tran > 0)  %  minimize transition bandwidth by bisection search
    f_tran_center = (f(3)+f(2))/2;
    df_top = (f(3)-f(2))/2;  % filter can be designed
    df_bot = 0;              % filter cannot be designed   
    fprintf('filter design Solved with df = %6.3f \n', df_top);
        
    while 1
        df_mid = (df_top + df_bot)/2;
        fp = f_tran_center - df_mid;
        fs = f_tran_center + df_mid;  
        f_new = [-fp, fp, fs, f(4)];
        [h0, status0] = fir_ap_cvx(n, f_new, a, d,lambda);
        if strcmp(status0, 'Failed')  % n_mid is still too low
            df_bot = df_mid;
            fprintf('filter design Failed with df = %6.3f \n', df_mid);
        else  % n_mid works, try something below it
            h = h0;  status = status0;
            df_top = df_mid;             
            fprintf('filter design Solved with df = %6.3f \n', df_mid);
        end
            
        if (df_top - df_bot) < df_thre                
            break;
        end
    end 
end


if min_tran  == 0
    fprintf('Finally solved with df = %6.3f \n', (f(3)-f(2))/2);
elseif (min_tran <= 1) & (min_tran > 0)
    df_new = ((f(3)-f(2))/2)*(1-min_tran) + df_top*min_tran;
    fp = f_tran_center - df_new;
    fs = f_tran_center + df_new;  
    f_new = [-fp, fp, fs, f(4)];
    [h, status] = fir_ap_cvx(n, f_new, a, d,lambda);
    fprintf('Finally solved with df = %6.3f \n', df_new);
    f = f_new;
else
    error('invalid input of min_tran');
end


fprintf('\n');
% determine order
if (min_order > 0)  %  minimize order by bisection search
    n_top = n;  % filter can be designed
    n_bot = 2;  % filter cannot be designed
    while 1
        n_mid = ceil((n_top + n_bot)/2);
        [h0, status0] = fir_ap_cvx(n_mid, f, a, d,lambda);
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
elseif min_order  == 0
    fprintf('Finally solved with n = %d \n', n);
elseif (min_order < 1) & (min_order > 0)
    n_new = ceil(n*(1-min_order) + n_top*min_order);
    fprintf('Finally solved with n = %d \n', n_new);
    [h, status] = fir_ap_cvx(n_new, f, a, d,lambda);
else
    error('invalid input of min_order');
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

                
        
        
        
            