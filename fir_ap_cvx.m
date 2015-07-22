function [h, status] = fir_ap_cvx(n, f, a, d, obj, Peak, dbg)
% fir_ap_cvx - FIR filter design using cvx
%
% function [h, status] = fir_ap_cvx(n, f, a, d, obj, Peak, dbg)
%   
% Design n-tap arbitrary-phase filter that meets multiband frequency
% respones manitude specification. 
% Additional function: 
%   (1) constraint of spikes(peaks) at two ends
%   (2) multiobjective function of stopband attenuation and total energy
%
% Based on "FIR Filter Design via Spectral Factorization and Convex 
% Optimization" by S.-P. Wu, S. Boyd, and L. Vandenberghe
%
% Inputs: --- similar to cfirpm
%   n: number of taps returned
%   f: frequency bands (-1->1)
%   a: amplitude at band edges
%   d: ripple  in bands
%   obj: parameter to determine trade-off between two objective functions
%   Peak: constraint on peak amplitude, especially at ends
%   dbg: flag to turn on debugging statements/plots
%
% Outputs: 
%   h: filter coefficients
%   status: 'Solved' or 'Failed'; 
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  June 2014

if nargin < 4,   error('not enough input');  end;
if nargin <= 4,  obj = 0;     end;
if nargin <= 5,  Peak = 1e-3; end;
if nargin <= 6,  dbg = 0;     end;


% internal parameters for constraint
% positive constraint on spectrum, to make sure S>=0 at all frequency
epsilon = 1e-10;


% Create optimization arrays
f = f * pi;				% Scale to +/- pi
oversamp = 15;
m = 2 * n * oversamp;
w = linspace(-pi,pi,m); 
w = sort([w f]);  % Add explicit samples to w at the edge of each specified band

% Find indices to passbands/stopbands, and fill in upper/lower bounds
idx_band = []; U_band = []; L_band = [];
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
idx_band = [1:length(wband)];                                                                 
wtran = w(idx_tran);                                                                          
idx_tran = [1:length(wtran)] + length(wband);                                                 
w = [wband(:).' wtran(:).'];                                                                  
m = size(w,2);   

if dbg >= 3, 
	figure; plot(wband,U_band,'b.'); hold on;
	plot(wband,L_band,'b.'); plot(wtran,U_tran,'r.'); plot(wtran,L_tran,'r.');
	hold off; title('upper and lower bound on magnitude');
end;

% create matrix A used to compute the power spectrum
A = [ones(m,1), 2*cos(kron(w',[1:n-1])), 2*sin(kron(w',[1:n-1]))];

% Build matrix for upper bound constraints
A_U = [A(idx_band,:); A(idx_tran,:)];
U_b = [U_band  U_tran];
U_b = U_b.^2;;

% Build matrices for lower bound constraints
A_L = [A(idx_band, :); A(idx_tran,:)];
L_b = [L_band  L_tran];
idx_lb = find( L_b < 0);
L_b(idx_lb) = 0;
L_b = L_b.^2;

% positive constraint
idx_lb = find(L_b < epsilon^2);
L_b(idx_lb) = epsilon^2;

% Combine matrices
A_b = [A_U; -A_L];
b = [U_b -L_b];

% Set linear objective function, add up spectrum energy in transition band 
% and/or stopband, to minimize total energy and stopband ripple respectively.
fmin_tran = sum(A(idx_tran,:), 1);
idx_stop = find( sqrt(U_b) < (min(sqrt(U_b))+1e-2) );

if dbg >= 3
    hold on;
    plot(wtran,0.8*ones(size(wtran)),'g.');  
	plot(w(idx_stop),0.2*ones(size(w(idx_stop))),'k.');
    legend('transition band for objective function 1','stopband for objective function 2' ); hold off;
end
  
% matrix to pick out r(n-1) for spike constraint
Famp_array = {};
Fi = zeros(1,2*n-1); Fi(1) = 1; Famp_array{1} = Fi;
for i=2:n
    Fi = zeros(2,2*n-1);
    Fi(1,i) = 1;
    Fi(2,n+i-1) = 1;
    Famp_array{i} = Fi;
end

% Call minimization routine
if obj >= 0  % multi-objective 
    
    % minimize energy at transition bands
%     cvx_begin quiet
%         variables  x(2*n-1)  ripple_stop    % Peak is actually Peak^2
%         minimize fmin_tran*x  + obj*ripple_stop  %+ obj*Peak   
%         subject to
%             A_b * x <= transpose(b)
%             A_U(idx_stop,:) * x <= ripple_stop
%             for i=1:n    
%                 norm(Famp_array{i} * x,2) <= (n-i+1)*Peak
%             end
%     cvx_end
    
    % minimize total energy
    cvx_begin quiet
        variables  x(2*n-1)  ripple_stop   
        minimize x(1)  + obj*ripple_stop 
        subject to
            A_b * x <= transpose(b)
            A_U(idx_stop,:) * x <= ripple_stop
            for i=1:n    
                norm(Famp_array{i} * x,2) <= (n-i+1)*Peak
            end
    cvx_end
    
else
    error('invalid input of obj');
end
    
    
if ( isequal(cvx_status, 'Solved') | isequal(cvx_status, 'Inaccurate/Solved') )
    status = 'Solved';
else
    status = 'Failed'; 
    h = [];
    return;
end;

% reshape variable x to autocorrelation r 
r = [x(1);  x(2:n) + sqrt(-1) * x((n+1):(2*n-1))];
r = [conj(r(end:-1:2)); r];
if dbg >= 3
    figure;
    plot(-(n-1):(n-1), real(r),'r-',-(n-1):(n-1), imag(r),'b-'); legend('real','imag'); title('Autocorrelation sequence');
end

% check calculation of |r(n)|
if dbg >= 3
    for i=1:n
        r_test(i) = norm(Famp_array{i} * x,2);
    end
    figure;
    plot(-(n-1):(n-1), abs(r),'r-*',0:(n-1), r_test,'b-s', 0:(n-1), [n:-1:1]*Peak,'k-.'); legend('|r(n)|','F*x','Constraint'); title('Autocorrelation magnitude');
end

% mp filter by spectral factorization
h = fmp2(r);
if dbg >= 2
    figure;
    plot(1:length(h), real(h),'r-*',1:length(h), imag(h),'b-*'); legend('real','imag'); title('filter coefficients');
end

% check zero-pole
if dbg >= 3
    Z1 = roots(r);
    P1 = zeros(length(Z1),1);
    Z2 = roots(h);
    P2 = zeros(length(Z2),1);
    figure;  
    subplot(1,2,1); zplane(Z1, P1); title('zero-pole of autocorrelation');
    subplot(1,2,2); zplane(Z2, P2); title('zero-pole of filter');
end

% check magnitude constraint
if dbg >= 2
	S = A * x;  % spectrum, |H(e^(jw))|^2  
    figure;
    [wsort, sidx] = sort(w);  plot(w(sidx), U_b(sidx), 'r-',w(sidx), L_b(sidx), 'r-');  hold on;
    plot(w(sidx), S(sidx)); title('lower/upper bound with Spectrum');  hold off; set(gca,'xlim',[-pi pi]);
    
	figure; plot_spec(f, a, d); hold on;
	[wsort, sidx] = sort(w); plot(w(sidx), sqrt(S(sidx)));  title('square root of Spectrum');  hold off; set(gca,'xlim',[-pi pi]);
     
    figure; 
    wdisp = linspace(-pi,pi,1e3); H = fftshift(fft(h, length(wdisp)));
    subplot(2,1,1); plot_spec(f, a, d); hold on; plot(wdisp, abs(H));  title('frequency response magnitude');  hold off; set(gca,'xlim',[-pi pi]);
    subplot(2,1,2); plot(wdisp, unwrap(angle(H)));  title('frequency response phase'); set(gca,'xlim',[-pi pi]);
end


return;
end





%  fft wrt the center of the array, instead of the first sample
%  written by John Pauly, 1992
%  (c) Board of Trustees, Leland Stanford Junior University
function y=fftc(x)
y = fftshift(fft(fftshift(x)));
end


%  Generate a minimum phase filter by Spectral factorization
%  written by John Pauly, 1992
%  (c) Board of Trustees, Leland Stanford Junior University
%  modified by Hong Shang, withnot adding delta2 as in equal-ripple pm
%  filter design
function hmp = fmp2(h)
h = transpose(h(:));

l = length(h);
if rem(l,2) == 0,
   disp('filter length must be odd');
   return;
end;
lp = 8*exp(ceil(log(l)/log(2))*log(2));
hp = [zeros(1,ceil((lp-l)/2)) h zeros(1,floor((lp-l)/2))];
hpf = fftc(hp);
% hpfs = hpf-min(real(hpf))*1.000001;  
% add ripple of stopband for equal-ripple pm filter design, but this is not
% necessary in my case when the spectrum is already positive, which makes 
% the spectral factorization result less accurate

hpfs = hpf;
hpfmp = mag2mp(sqrt(abs(hpfs)));
hpmp = ifft(fftshift(conj(hpfmp)));
hmp = hpmp(1:(l+1)/2);
end


% mag2mp - take the magnitude of the fft of a signal, and return the
%   fft of the analytic signal.
% as = mag2mp(ms)
%   as - fft of analytic signal
%   ms - magnitude of analytic signal fft.
%  Written by John Pauly, Oct 1989
% (c) Board of Trustees, Leland Stanford Junior University
function [a] = mag2mp(x)
n = length(x);
xl = log(x);                        % log of mag spectrum
xlf = fft(xl);                      % 
xlfp(1) = xlf(1);                   % keep DC the same
xlfp(2:(n/2)) = 2*xlf(2:(n/2));     % double positive freqs
xlfp((n/2+1)) = xlf((n/2+1));       % keep half Nyquist the same,too
xlfp((n/2+2):n) = 0*xlf((n/2+2):n); % zero neg freqs
xlaf = ifft(xlfp);                  % 
a = exp(xlaf);                      % complex exponentiation
end

