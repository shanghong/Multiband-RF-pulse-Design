function h_new = fir_flip_zero(h, dbg)
%
%  h_new = fir_flip_zero(h, dbg)
%         
%  Flip zeros of a filter to get all the possible filters with same 
%  frequency respones magnitude but different phase, among which choose 
%  the one with min peak value
%    
%  Inputs: 
%    h: original filter
%    dbg: flag to turn on debugging statements/plots
%    
%  Outputs:
%    h: new filter 
%
%   modified to implement random sampling in case that order of filter is too high
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  July 2013


% zeros and poles
N = length(h);         % order of filter
Z = roots(h);
P = zeros(N-1,1);

idx_pb = find( abs(Z)>(1+1e-2) | abs(Z)<(1-1e-2) );  % select the passband zero
N_z = length(idx_pb);            % number of passband zero
Z_pb = Z(idx_pb);  
P_pb = P(idx_pb);
Z_flip = flip_zero(Z); % flip all the zeros
Z_pb_flip = flip_zero(Z_pb); % just flip passband zeros


if dbg >= 2    
    figure; 
    subplot(2,2,1); zplane(Z, P); title('diagram of zero-pole');
    subplot(2,2,2); zplane(Z_pb, P_pb); title('diagram of only passband zero-pole'); 
    subplot(2,2,3); zplane(Z_flip, P);    title('flipped all the zeros');
    subplot(2,2,4); zplane(Z_pb_flip, P_pb); title('flipped passband zero');
end


if N_z <= 19
    mask = combination_2power(N_z); % each column corresponds to one unique combination
                                    % cover all the possible combinations
                                    % 1 flip, 0 not flip                                
    % limit number of different combinations to be tried. If N_z is too 
    % large, random sampling to reduce computation and get suboptimal result 
    if N_z <= 12
        Num = 2^N_z;   % try all the cases
    else
        Num = 2^12;    % try only part of the cases
        idx_rand = randperm(2^N_z);
        idx_rand = idx_rand(1:Num);
        idx_rand = sort(idx_rand);
        mask = mask(:,idx_rand);
    end
else  % N_z is too large that full mask cannot be calculated, use monte-carlo method
    Num = 2^12;   
    mask = combination_MC(N_z,Num); % each column may not be unique
end

h_array = zeros(N,Num);       % array to store reproduced filter
Power_Total = zeros(1,Num);   % array to store filter sum-of-square
peak = zeros(1,Num);          % array to store filter peak amplitude 

for i=1:Num;
     
    Z_each = Z;
    Z_each(idx_pb) = Z_pb.*(1-mask(:,i)) + Z_pb_flip.*mask(:,i);  % set passband zero
    h_each = poly(Z_each);  % non linear phase filter coefficient
    h_each = h_each*sum(h)/sum(h_each); % scale   
    h_array(:,i) = transpose(h_each); 
    
    Power_Total(i) = sum((abs(h_each)).^2);
    peak(i) = max(abs(h_each));
    
    if dbg >= 3
        h_zp = [transpose(h_each(:)),zeros(1,length(h_each)*127)]; % zero padding for spectrum interpolation
        H = fftshift(fft(h_zp));          % FFT to get frequency response, which is periodic with 2*pi, or BW
        w = linspace(-0.5,0.5,length(H)); % relative frequency, f/fs 

        figure(300);
        subplot(1,3,1); zplane(Z_each, P); title('zero pole'); set(gca,'xlim',[-3,3]); set(gca,'ylim',[-3,3]);
        subplot(1,3,2); plot(1:N,real(h_each),'b*-',1:N,imag(h_each),'r*-'); legend('real','imag'); xlabel('sample'); title('reproduced filter coefficient'); axis tight; 
        subplot(1,3,3); plot(w,abs(H),'b-'); title('magnitude frequency response'); axis tight; 
        drawnow; %pause(0.5);
    end
end;

% optimize peak amplitude
[~,min_idx] = min(peak); % min of max amplitude, as the optimal case
Z_min = Z;
Z_min(idx_pb) = Z_pb.*(1-mask(:,min_idx)) + Z_pb_flip.*mask(:,min_idx);  % set passband zero
h_new = h_array(:,min_idx);

if dbg >= 1
    fprintf('reduce peak amplitude from %6.4f to %6.4f by %6.4f\n', max(abs(h)), max(abs(h_new)), (max(abs(h))-max(abs(h_new)))/ max(abs(h)) );
end

if dbg >= 2
    figure;
    subplot(2,1,1); plot(Power_Total); xlabel('type'); title('total power');   axis tight;
    subplot(2,1,2); plot(peak);        xlabel('type'); title('peak amplitude');axis tight
    
    figure;
    subplot(2,2,1); zplane(Z, P); title('initial zero-pole'); %set(gca,'xlim',[-3,3]); set(gca,'ylim',[-3,3]);
    subplot(2,2,2); plot(1:N,real(h),'b*-',1:N,imag(h),'r*-'); legend('real','imag'); xlabel('sample'); title('initial filter coefficient'); axis tight; 
    subplot(2,2,3); zplane(Z_min, P); title('flipped zero-pole'); %set(gca,'xlim',[-3,3]); set(gca,'ylim',[-3,3]);
    subplot(2,2,4); plot(1:N,real(h_new),'b*-',1:N,imag(h_new),'r*-'); legend('real','imag'); xlabel('sample'); title('optimized filter coefficient'); axis tight; 
end

end

function z_output = flip_zero(z_input)
% z_output = flip_zero(z_input)
% flip an array of zeros
i = sqrt(-1);
z_output = (1./abs(z_input)) .* exp(i.*angle(z_input));
end
         
function A = combination_2power(n)
	% function A = combination_2power(n) 
    % to search for all the 2^n combination with backtracking algorithm
    %
    % input:
    %   n - 2^n combination
    %
    % output:
    %   A - n*2^n matrix, each column corresponding to a type of combination
    %
    % Written by Hong Shang, May 2031, in John Pauly's RF class 
    
    if n == 1
        A = [1 0];
    else
        A_lower = combination_2power(n-1);
        A = [ [ones(1,size(A_lower,2)); A_lower], [zeros(1,size(A_lower,2)); A_lower]];
    end
end


  
function A = combination_MC(n,m)
	% function A = combination_MC(n) 
    %
    % To generate samples from all 2^n combination with monte-carlo random
    % sampling. Reduce computation especially when n is large. Each sample 
    % may not be unique, but when n is large, this probability is very low.
    %
    % input:
    %   n - 2^n combination totally
    %   m - choose m from them
    % output:
    %   A - n*m matrix, each column corresponding to a type of combination,
    %   which may not be unique
    %
    % Written by Hong Shang, June 14 
    
    A = zeros(n,m);
    for i=1:m
        a = rand(n,1);
        a = round(a);
        A(:,i) = a;
    end
end