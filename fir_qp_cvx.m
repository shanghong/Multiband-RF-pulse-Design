function [h, status] = fir_qp_cvx(n, f, a, d, k, obj, dbg)
% fir_qp_cvx - FIR filter design with quadratic phase using cvx
%   
% Design n-tap quadratic-phase filter that meets multiband frequency
% respones manitude specification. 
% Additional function: 
%   (1) multi-objective function of total energy and peak amplitude
%
% [h, status] = fir_qp_cvx(n, f, a, d, k, obj, dbg)
%
% Inputs: --- similar to cfirpm
%   n: number of taps returned
%   f: frequency bands (-1->1)
%   a: amplitude at band edges
%   d: ripple  in bands
%   k: specify the amount of quadratic phase
%   obj: parameter to determine trade-off between two objective functions
%   dbg: flag to turn on debugging statements/plots
%
% Outputs: 
%   h: filter coefficients
%   status: 'Solved' or 'Failed'; 
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  June 2014

if nargin < 4,  error('not enough input');  end;
if nargin == 4,  k = 100;  end;
if nargin <= 5,  obj = 0; end;
if nargin <= 6,  dbg = 0; end;

% Create optimization arrays
f = f * pi;				% Scale to +/- pi
oversamp = 10;
m = n * oversamp;
w = linspace(-pi,pi,m); 
w = sort([w f]);  % Add explicit samples to w at the edge of each specified band

% Find indices to passbands/stopbands, and fill in upper/lower bounds
idx_band = []; 
U_band = [];  L_band = [];  
M_band = [];  D_band = [];
nband = length(f)/2;
for band = 1:nband, 
	idx = find( (w >= f(band*2-1)) & (w <= f(band*2)) );
	% Get amplitude from linear interpolation on band
	idx_band = [idx_band idx];
	if (f(band*2-1) == f(band*2))
	    amp = a(band*2-1);
	else
	    amp = a(band*2-1) + (a(band*2)-a(band*2-1)) * ((w(idx) - f(band*2-1))/(f(band*2)-f(band*2-1)));
	end;
	U_band = [U_band (amp + d(band))];
	L_band = [L_band (amp - d(band))];
    M_band = [M_band  amp];
    D_band = [D_band  d(band)*ones(size(amp))];
end; 

% Get transition indices
idx_tmp = ones(1,length(w));
idx_tmp(idx_band) = 0;
idx_tran = find(idx_tmp == 1);


% Add transition band limits to be between the max specification on each 
% band and min of (0,min(L_band))
if ~isempty(idx_tran)
	U_amp_tran = max(U_band);
	U_tran = U_amp_tran*ones(1,length(idx_tran));
	L_amp_tran = min(0, min(L_band));
	L_tran = L_amp_tran*ones(1,length(idx_tran));
else
	U_tran = [];
	L_tran = [];
end;
    
% Update w, idx_band
wband = w(idx_band);
wtran = w(idx_tran);

if dbg >= 3, 
	figure;  hold on;
    plot(w(idx_band),U_band,'b.');
	plot(w(idx_band),L_band,'b.');
    plot(w(idx_band),M_band,'k.');
    plot(w(idx_tran),U_tran,'r.'); 
    plot(w(idx_tran),L_tran,'r.');
	hold off; title('upper and lower bound on magnitude');
end;

% matrix A used for approximating desired frequency response in minmax sense
% A is R^(2*2N) at each frequency
% use a cell array to store matrix A at each frequency
% only frequency in band is considered for calculation
A_array = {};
for i=1:length(wband)
    wi = wband(i);
    Ai = [ [cos(wi*[0:(n-1)]), sin(wi*[0:(n-1)])];  [-sin(wi*[0:(n-1)]), cos(wi*[0:(n-1)])] ];
    A_array{i} = Ai;
end

% limit transition band also
A_array_tran = {};
for i=1:length(wtran)
    wi = wtran(i);
    Ai = [ [cos(wi*[0:(n-1)]), sin(wi*[0:(n-1)])];  [-sin(wi*[0:(n-1)]), cos(wi*[0:(n-1)])] ];
    A_array_tran{i} = Ai;
end


% desired frequency response
Hd_array = {};
Hd_array_complex = [];
j = sqrt(-1);
for i=1:length(wband)
    wi = wband(i);
    Hdi = M_band(i) * exp( j*(k*wi^2 - wi*(n-1)/2) );
    Hd_array{i} = [real(Hdi); imag(Hdi)];
    Hd_array_complex = [Hd_array_complex, Hdi];
end


if dbg >= 3
    figure;  
    subplot(2,2,1); plot(wband, abs(Hd_array_complex),'b.');           xlabel('w'); title('desired manitude');
    subplot(2,2,3); plot(wband, unwrap(angle(Hd_array_complex)),'b.'); xlabel('w'); title('desired phase');
    subplot(2,2,2); plot(wband, real(Hd_array_complex),'b-');          xlabel('w'); title('desired real part');
    subplot(2,2,4); plot(wband, imag(Hd_array_complex),'b-');          xlabel('w'); title('desired imag part');
end

% Matrix F to pick out real and imag part of h(i) from x
F_array = {};
for i=1:n
    Fi = zeros(2,2*n);
    Fi(1,i) = 1;
    Fi(2,n+i) = 1;
    F_array{i} = Fi;
end


% Call minimization routine
if length(obj) == 1  % direct constraint on frequency response 
    
    cvx_begin %quiet
        variables x(2*n)  E_total   Peak
        minimize E_total + obj*Peak
        subject to
            % costraint on pass/stop band
            for i=1:length(wband)
                norm(A_array{i}*x - Hd_array{i},2) <= D_band(i);
            end
            
            % also limit transition band magnitude
            for i=1:length(wtran)
                norm(A_array_tran{i}*x,2) <= 1 + max(d)*5;
            end
            
            % constraint on peak amplitude
            for i=1:n
                norm(F_array{i}*x,2) <= Peak;
            end
            
            % constraint on total energy
            norm(x,2) <= E_total;
    cvx_end
    
elseif length(obj) == 2  % minimax on frequency response
    
    cvx_begin %quiet
        variables x(2*n)  delta  E_total   Peak
        minimize delta + obj(1)*E_total + obj(2)*Peak
        subject to
            % costraint on pass/stop band
            for i=1:length(wband)
                norm(A_array{i}*x - Hd_array{i},2) <= D_band(i)*delta;
            end
            
            % also limit transition band magnitude
            for i=1:length(wtran)
                norm(A_array_tran{i}*x,2) <= 1 + 0.1;
            end
            
            % constraint on peak amplitude
            for i=1:n
                norm(F_array{i}*x,2) <= Peak;
            end
            
            % constraint on total energy
            norm(x,2) <= E_total;
    cvx_end
    
    
else
    error('invalid input of obj');
end
    
    

if ( isequal(cvx_status, 'Solved') | isequal(cvx_status,'Inaccurate/Solved') )
    status = 'Solved';
else
    status = 'Failed'; 
    h = [];
    return;
end;


h = x(1:n) + sqrt(-1)*x((n+1):end);
% check filter coefficients
if dbg >= 2
    figure;
    subplot(1,2,1); plot(1:length(h), abs(h),'r-*');           legend('magnitude'); title('filter coefficients');
    subplot(1,2,2); plot(1:length(h), real(h),'r-*',1:length(h), imag(h),'b-*'); legend('real','imag'); title('filter coefficients');
end

% check zero-pole
if dbg >= 3
    Z1 = roots(h);
    P1 = zeros(length(Z1),1);
    figure;   zplane(Z1, P1); title('zero-pole of filter');
end

% check frequency response
if dbg >= 2	
    figure; 
    wdisp = linspace(-pi,pi,1e3); H = fftshift(fft(h, length(wdisp)));
    H_phase = unwrap(angle(H));
    H_phase = H_phase + transpose(wdisp*(n-1)/2);  % compensate the linear phase due to shifted h(n)
    H_phase = H_phase - H_phase(ceil(end/2));      % define zero phase at w = 0
    subplot(2,1,1); plot(wdisp, abs(H));  hold on; plot_spec(f, a, d); hold off; title('frequency response magnitude'); axis tight; 
    subplot(2,1,2); plot(wdisp, H_phase );  title('frequency response phase'); axis tight;
    
    figure; 
    idx = find( (wdisp<(f(5)+0.1)) & (wdisp>(f(2)-0.1)) );
    subplot(2,1,1); plot( wdisp(idx), abs(H(idx)) );  title('frequency response magnitude'); axis tight;  hold on;  plot_spec(f, a, d); hold off; 
    subplot(2,1,2); plot( wdisp(idx), H_phase(idx), 'r-', wdisp(idx), k*(wdisp(idx).^2) ,'b-'); legend('actual','desired'); title('frequency response phase');  axis tight;  hold on;  plot_spec(f, zeros(size(a)), zeros(size(d))); hold off; 
    
end


return;

