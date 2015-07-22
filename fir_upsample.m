function h_new = fir_upsample(h, dt1, dt2, dbg)
%
%  h_new = fir_upsample(h, dt1, dt2, dbg)
%
%  upsampling filter 
%
%  Inputs:
%    h: filter coefficient
%    dt1: in ms, previous sampling interval
%    dt2: in ms, new sampling interval
%    dbg: show debug plot
%
%  Outputs:
%    h_new: filter after upsampling
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  June 2014

n = round(dt1/dt2);

% ======================== method 1 ========================
% use matlab function 'resample'
% easy-to-use, which does not require you to supply a filter or compensate 
% for the signal delay introduced by filtering. The filter is designed 
% using FIRLS. Large deviations from zero at the end points of the sequence
% can cause inaccuracies.

h_new = resample(h,n,1,floor(length(h)/2));
h_new = h_new / n;

% ======================== method 2 ========================
%truncate frequency response directly
% h = transpose(h(:));
% h_us = upsample(h,n);
% N = 1e4;
% N_zeros = round((N-length(h_us))/2);
% h_zp = [zeros(1, N_zeros), h_us, zeros(1, N_zeros)];
% w = linspace(-pi,pi,length(h_zp));
% H_hr = fftc(h_zp);
% idx = find( (w>=(-pi/n)) & (w<=(pi/n)) );
% H_trunc = zeros(size(H_hr));
% H_trunc(idx) = H_hr(idx);
% h_new = ifftc(H_trunc);
% h_new(1:N_zeros) = [];
% h_new((end-N_zeros+1):end) = [];

%h_new = h_new .* transpose(hann(length(h_new)));

% display for debug
if dbg >= 2
    figure;
    subplot(1,2,1); plot((1:length(h))*dt1,real(h),'b-s', (1:length(h_new))*dt2, real(h_new),'r-*'); leg = legend('original','upsampled'); set(leg,'FontSize',17); xlabel('time, ms','FontSize',18); title('filter coefficients Real part','FontSize',19);      axis tight; set(gca,'FontSize',17); 
    subplot(1,2,2); plot((1:length(h))*dt1,imag(h),'b-s', (1:length(h_new))*dt2, imag(h_new),'r-*'); leg = legend('original','upsampled'); set(leg,'FontSize',17); xlabel('time, ms','FontSize',18); title('filter coefficients Imaginary part','FontSize',19); axis tight; set(gca,'FontSize',17);    
    set(gcf, 'Position', [50, 50, 1100, 400], 'PaperPositionMode', 'auto');
    
    w = linspace(-pi,pi,1e3);
    H = fftshift(fft(h, length(w)));
    H_new = fftshift(fft(h_new, length(w)));
    
    figure;
    subplot(1,2,1); plot(w, abs(H),'b-', w, abs(H_new),'r-'); leg = legend('original','upsampled'); set(leg,'FontSize',17); xlabel('w','FontSize',18);  title('H(w) magnitude','FontSize',19); axis tight;  set(gca,'FontSize',17); 
    subplot(1,2,2); plot(w, unwrap(angle(H)),'b-', w, unwrap(angle(H_new)),'r-'); leg = legend('original','upsampled'); set(leg,'FontSize',17);  xlabel('w','FontSize',18);  title('H(w) phase','FontSize',19); axis tight; set(gca,'FontSize',17); 
    set(gcf, 'Position', [50, 50, 1100, 400], 'PaperPositionMode', 'auto');
    
end
    
end
  
function F=fftc(f)
    F = ifftshift(fft(fftshift(f)));
end

function f=ifftc(F)
    f = fftshift(ifft(ifftshift(F)));
end
