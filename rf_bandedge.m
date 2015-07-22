function fn = rf_bandedge(n, dt, mb_cf, mb_range, mb_FA, mb_ripple, ptype, dbg)
%
%  For multiband SLR RF pulse design, figure out the band edge if not 
%  specified, with the following steps
%     (1)Ideal band edge is the middle between two center_frequency.
%     (2)Insert transition band, which is determined by analytical equation
%        in pm filter design, to make it comparable to conventional methods
%  This function is used in dzrf_mb.m, with similar input
%
%  Inputs:
%    n -- order of filter, number of hard pulse
%    dt -- in ms, sampling interval of RF pulse/filter
%    mb_cf --    1*m in kHz (cell array), center frequency for each band,
%                    a point (legnth(mb_cf{i}) = 1) or a range (legnth(mb_cf{i}) = 2) 
%    mb_range -- 1*m in kHz, width for each band, [-mb_range/2, mb_range/2]
%                if empty, figure out automatically
%    mb_FA --    1*m in degree, flip-angle for each band
%    mb_ripple-- 1*m, absolute ripple of magnetization for each band
%    ptype -- pulse type.  Options are: st, ex, se, sat, inv                    
%    dbg -- debug level
%
%  Outputs:
%    fn -- spec of frequency range for each band, [-1 1]
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Sep 2014


T = n*dt;               % in ms, pulse duration
m_band = length(mb_FA); % number of band
fs = 1/dt;              % in kHz, sampling rate
fn = [];

if isempty(mb_range)  
% ================================================================
% band edge is NOT specified
    % figure out the transition bandwidth df for pm filter design
    % for each type of pulse, assume the standard FA exist, like FA = pi/2,
    % which is used in dzrf.m with analytical equation. But such standard
    % FA may not actually exist. Use the closest one approximately.
    switch ptype
    case {'st','ex'}
        % pick out in-slice approximately, FA = pi/2
        [~, i] = min(abs(mb_FA - 90));
        d1 = mb_ripple(i);

        % pick out out-of-slice approximately, FA = 0
        [~, i] = min(abs(mb_FA));
        d2 = mb_ripple(i);
        
        delta1 = sqrt(d1/2);
        delta2 = d2/sqrt(2);
               
    case 'inv'
        % pick out in-slice approximately, FA = pi
        [~, i] = min(abs(mb_FA - 180));
        d1 = mb_ripple(i);

        % pick out out-of-slice approximately, FA = 0
        [~, i] = min(abs(mb_FA));
        d2 = mb_ripple(i);
        
        delta1 = d1/8;
        delta2 = sqrt(d2/2);
  
    case 'sat'
        % pick out in-slice approximately, FA = pi/2
        [~, i] = min(abs(mb_FA - 90));
        d1 = mb_ripple(i);

        % pick out out-of-slice approximately, FA = 0
        [~, i] = min(abs(mb_FA));
        d2 = mb_ripple(i);
        
        delta1 = d1/2;
        delta2 = sqrt(d2);
   
    case 'se'
        % pick out in-slice approximately, FA = pi
        [~, i] = min(abs(mb_FA - 180));
        d1 = mb_ripple(i);

        % pick out out-of-slice approximately, FA = 0
        [~, i] = min(abs(mb_FA));
        d2 = mb_ripple(i);
        
        delta1 = d1/4;
        delta2 = sqrt(d2);
        
    otherwise
        disp(['Unrecognized Pulse Type -- ',ptype]);
        disp('Recognized types are st, ex, se, inv, and sat');
        return;
    end

    df = dinf(delta1, delta2)/T; % in kHz, transition width

    % initial band center
    f = zeros(1,2*m_band); % in kHz, frequency band edge, each band monotonically increasing 
    for i=1:m_band
        if length(mb_cf{i}) == 1 % define the band around one point
            f(2*i-1) = mb_cf{i};
            f(2*i) = mb_cf{i};
        else                     % define the band around one range
            mb_cf_range = mb_cf{i};
            f(2*i-1) = mb_cf_range(1);
            f(2*i) = mb_cf_range(2);
        end
    end

    % expand to ideal band edge
    f1 = f; % in kHz
    for i=1:(m_band-1)
        f1(2*i)  =  (f(2*i) + f(2*i+1))/2;
        f1(2*i+1) = (f(2*i) + f(2*i+1))/2;
    end
    f1(1) = f(1) - (f1(2)-f(2));
    f1(end) = f(end) + (f(end-1) - f1(end-1));

    % insert transition band
    f2 = f1; % in kHz
    for i=1:m_band
        f2(2*i-1) = f1(2*i-1) + df/2;
        f2(2*i) = f1(2*i) - df/2;   
    end
    
    f = f2;% in kHz

else  
% ================================================================
% band edge is specified
    f = zeros(1,2*m_band); % in kHz, frequency band edge, each band monotonically increasing 
    for i=1:m_band
        if length(mb_cf{i}) == 1 % define the band around one point
            f(2*i-1) = mb_cf{i} - mb_range(i)/2;
            f(2*i) = mb_cf{i} + mb_range(i)/2;
        else                     % define the band around one range
            mb_cf_range = mb_cf{i};
            f(2*i-1) = mb_cf_range(1) - mb_range(i)/2;
            f(2*i) = mb_cf_range(2) + mb_range(i)/2;
        end
    end 
end


% check overlap
[f_sort, i_f_sort] = sort(f,'ascend');
if norm(i_f_sort - [1:2*m_band]) > 1e-10  % 
    disp('Incompatible spec of frequency range.');
    disp('f is not monotonically increasing.');
    error('Try reducing mb_range (width for each band) or reorder bands. ');
end

% check sampling rate
if ( (f(1) < (-fs/2)) || (f(end) > (fs/2)) )
    error('the sampling rate is not enough, increase n');
end

% normalized frequency in [-1 1]
fn = f/(fs/2); 

% debug display
if dbg >= 2
    if isempty(mb_range) 
        figure;
        subplot(3,1,1); plot_spec(f1, ones(1,2*m_band), linspace(0.1,0.2,m_band));  xlabel('frequency, kHz','FontSize',18);       title('Ideal Band Edge','FontSize',20);  set(gca,'FontSize',18); axis tight; set(gca,'ylim',[1-0.3 1+0.3]); set(gca,'xlim',[min(f1)-0.1, max(f1)+0.1]);
        subplot(3,1,2); plot_spec(f2, ones(1,2*m_band), linspace(0.1,0.2,m_band));  xlabel('frequency, kHz','FontSize',18);       title('Actual Band Edge','FontSize',20); set(gca,'FontSize',18); axis tight; set(gca,'ylim',[1-0.3 1+0.3]); set(gca,'xlim',[min(f2)-0.1, max(f2)+0.1]);
        subplot(3,1,3); plot_spec(fn, ones(1,2*m_band), linspace(0.1,0.2,m_band));  xlabel('normalized frequency','FontSize',18); title('Actual Band Edge','FontSize',20); set(gca,'FontSize',18); axis tight; set(gca,'ylim',[1-0.3 1+0.3]); set(gca,'xlim',[-1, 1]);
        set(gcf, 'Position', [500, 50, 800, 800], 'PaperPositionMode', 'auto');

    else
        figure;
        subplot(2,1,1); plot_spec(f, ones(1,2*m_band), linspace(0.1,0.2,m_band));   xlabel('frequency, kHz','FontSize',18); title('Specified Band Edge','FontSize',20);       set(gca,'FontSize',18); axis tight; set(gca,'ylim',[1-0.3 1+0.3]); set(gca,'xlim',[min(f)-0.1, max(f)+0.1]);
        subplot(2,1,2); plot_spec(fn, ones(1,2*m_band), linspace(0.1,0.2,m_band));  xlabel('normalized frequency','FontSize',18); title('Specified Band Edge','FontSize',20); set(gca,'FontSize',18); axis tight; set(gca,'ylim',[1-0.3 1+0.3]); set(gca,'xlim',[-1, 1]);
        set(gcf, 'Position', [500, 50, 800, 600], 'PaperPositionMode', 'auto');
    
    end
end

