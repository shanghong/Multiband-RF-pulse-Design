function [range_B, ripple_B] = rf_ripple_GFA(FA, ripple_M, ptype, appro, dbg)
%
%  For SLR algorithm with generalized flip-angle, calculate spec of ripple 
%  for Beta polynomial given flip-angle, ripple of magnetization, pulse type
%
%  Inputs:
%    FA -- in degree, flip angle
%    ripple_M -- absolute ripple of magnetization, instead of relative 
%                ripple, or ripple for flip-angle   
%    ptype -- pulse type, options are:
%       st  -- small tip angle                    
%       ex  -- excitation pulse
%       sat -- saturation pulse 
%       inv -- inversion pulse
%       se  -- spin-echo pulse
%    appro -- 1 use quadratic equation approximation, 0 directly solve by asin
%    dbg -- debug level, 1 run self-test of the range of Magnetization profile
%
%  Outputs:
%    range_B -- range of amplitude of Beta 
%    ripple_B -- ripple of amplitude of Beta 
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Sep 2014

if nargin <= 3
    appro = 0;
end

if nargin <= 4
    dbg = 0;
end

if appro
    [range_B, ripple_B] = rf_ripple_quad(FA, ripple_M, ptype);
else
    [range_B, ripple_B] = rf_ripple_asin(FA, ripple_M, ptype);
end
    
% test the result
if dbg >= 1
    N0 = 1e3;
    switch ptype
        case 'st'
           
        case 'ex'
            B_sample = linspace(range_B(1), range_B(2), N0);
            Mxy_measure = 2 * sqrt(1-B_sample.^2) .* B_sample;
            M_range_measure = [min(Mxy_measure), max(Mxy_measure)];
            M_range_desire = [ sin(FA*pi/180)-ripple_M, min([sin(FA*pi/180)+ripple_M, 1])];
            fprintf(' desired Mxy [%0.4f, %0.4f] = [%0.4f-%0.4f, %0.4f+%0.4f]\n', M_range_desire(1), M_range_desire(2), sin(FA*pi/180), ripple_M, sin(FA*pi/180), ripple_M);
            fprintf('measured Mxy [%0.4f, %0.4f] \n',M_range_measure(1), M_range_measure(2));
            fprintf('       error [%0.4f, %0.4f] \n',M_range_measure(1)-M_range_desire(1), M_range_measure(2)-M_range_desire(2));

        case {'sat','inv'}
            B_sample = linspace(range_B(1), range_B(2), N0);
            Mz_measure = 1 - 2 * B_sample.^2 ;
            M_range_measure = [min(Mz_measure), max(Mz_measure)];
            M_range_desire = [ max([cos(FA*pi/180)-ripple_M, -1]), min([cos(FA*pi/180)+ripple_M, 1])];
            fprintf(' desired Mz [%0.4f, %0.4f] = [%0.4f-%0.4f, %0.4f+%0.4f]\n', M_range_desire(1), M_range_desire(2), cos(FA*pi/180), ripple_M, cos(FA*pi/180), ripple_M);
            fprintf('measured Mz [%0.4f, %0.4f] \n',M_range_measure(1), M_range_measure(2));
            fprintf('      error [%0.4f, %0.4f] \n',M_range_measure(1)-M_range_desire(1), M_range_measure(2)-M_range_desire(2));
            
        case 'se'
            B_sample = linspace(range_B(1), range_B(2), N0);
            Mxy_measure = B_sample.^2;
            M_range_measure = [min(Mxy_measure), max(Mxy_measure)];
            M_range_desire = [ max([sin((FA/2)*pi/180)^2-ripple_M, 0]), min([sin((FA/2)*pi/180)^2+ripple_M, 1])];
            fprintf(' desired Mxy [%0.4f, %0.4f] \n', M_range_desire(1), M_range_desire(2));
            fprintf('measured Mxy [%0.4f, %0.4f] \n',M_range_measure(1), M_range_measure(2));
            fprintf('       error [%0.4f, %0.4f] \n',M_range_measure(1)-M_range_desire(1), M_range_measure(2)-M_range_desire(2));

        otherwise
            disp(['Unrecognized Pulse Type -- ',ptype]);
            disp('Recognized types are st, ex, se, inv, and sat');
            return;
    end
end
    

end


function [range_B, ripple_B] = rf_ripple_quad(FA, ripple_M, ptype)
%  
% step 1: figure out range of flip-angle by solving quadratic equation
% step 2: calcualte range of Bn(z)
%
range_B = [];
ripple_B = [];
   
% check whether FA in the range
if ( (FA < 0) || (FA > 180) )
   disp(['Flip angle of ', num2str(FA), ' degree is out of range']);
   disp('Flip angle should be in the range of [0 180] degree');
   return;
end
FA = FA*pi/180;  % in radians


% figure out ripple of flip-angle
switch ptype
    case 'st'
        mid_M = sin(FA);
        max_M = min(1,mid_M+ripple_M);
        min_M = mid_M-ripple_M;        
        max_B = max_M;
        min_B = min_M;
        mid_B = mid_M;
        range_B = [min_B, max_B];
        ripple_B = abs( range_B - mid_B );
        
    case 'ex'
        if ((sin(FA)+ ripple_M) >= 1) 
            C_r = [-0.5*sin(FA), cos(FA),  ripple_M];
            C_l = [-0.5*sin(FA), -cos(FA), ripple_M];
        elseif (FA <= (pi/2))
            C_r = [-0.5*sin(FA), cos(FA), -ripple_M];
            C_l = [-0.5*sin(FA), -cos(FA), ripple_M];
        else
            C_r = [-0.5*sin(FA), cos(FA), ripple_M];
            C_l = [-0.5*sin(FA), -cos(FA), -ripple_M];
        end       
        
        % solve the quadratic equation will give two roots, only the positive, small-tip-angle root is the right solution
        dfa_r = roots(C_r);
        dfa_l = roots(C_l);
        dfa_r = dfa_r( find( (dfa_r > 0) ) );
        dfa_l = dfa_l( find( (dfa_l > 0) ) );
        dfa_r = min(dfa_r);
        dfa_l = min(dfa_l);
        
        % calculate ripple of |Bn(z)|
        [range_B, ripple_B] = rf_ripple_FA2Beta(FA + dfa_r, FA-dfa_l, FA);

    case {'sat','inv'}
        if ( (cos(FA)+ripple_M) > 1 )
            C_r = [-0.5*cos(FA), sin(FA),  ripple_M];
            C_l = [-0.5*cos(FA), -sin(FA), ripple_M];
        elseif ( (cos(FA)-ripple_M) < -1 )
            C_r = [-0.5*cos(FA), sin(FA),  -ripple_M];
            C_l = [-0.5*cos(FA), -sin(FA), -ripple_M];
        else
            C_r = [-0.5*cos(FA), sin(FA),  -ripple_M];
            C_l = [-0.5*cos(FA), -sin(FA),  ripple_M];
        end
        
        % solve the quadratic equation will give two roots, only the positive, small-tip-angle root is the right solution
        dfa_r = roots(C_r);
        dfa_l = roots(C_l);
        dfa_r = dfa_r( find( (dfa_r > 0) ) );
        dfa_l = dfa_l( find( (dfa_l > 0) ) );
        dfa_r = min(dfa_r);
        dfa_l = min(dfa_l);
        
        % calculate ripple of |Bn(z)|
        [range_B, ripple_B] = rf_ripple_FA2Beta(FA + dfa_r, FA-dfa_l, FA);

    case 'se'
        mid_M = (sin(FA/2))^2;
        max_M = max(0,min(1,mid_M+ripple_M));
        min_M = max(0,min(1,mid_M-ripple_M));        
        max_B = sqrt(max_M);
        min_B = sqrt(min_M);
        mid_B = sin(FA/2);
        range_B = [min_B, max_B];
        ripple_B = abs( range_B - mid_B );
        
    otherwise
        disp(['Unrecognized Pulse Type -- ',ptype]);
        disp('Recognized types are st, ex, se, inv, and sat');
        return;
end

  
end


function [range_B, ripple_B] = rf_ripple_asin(FA, ripple_M, ptype)
%  
% step 1: figure out range of flip-angle by asin 
%         (this should be more accurate than quadratic equation approximation)
% step 2: calcualte range of Bn(z)
% 

range_B = [];
ripple_B = [];
   
% check whether FA in the range
if ( (FA < 0) || (FA > 180) )
   disp(['Flip angle of ', num2str(FA), ' degree is out of range']);
   disp('Flip angle should be in the range of [0 180] degree');
   return;
end
FA = FA*pi/180;  % in radians


% figure out ripple of flip-angle
switch ptype
    case 'st'
        max_M = min(1,mid_M+ripple_M);
        min_M = mid_M-ripple_M;        
        max_B = max_M;
        min_B = min_M;
        mid_B = mid_M;
        range_B = [min_B, max_B];
        ripple_B = abs( range_B - mid_B );
        
    case 'ex'
        if ((sin(FA)+ ripple_M) >= 1) 
            rfa = asin(sin(FA) - ripple_M);
            rfa_l = rfa;
            rfa_r = pi - rfa;
        elseif (FA <= (pi/2))
            rfa_l = asin(sin(FA) - ripple_M);
            rfa_r = asin(sin(FA) + ripple_M);
        else
            rfa_r = pi - asin(sin(FA) - ripple_M);
            rfa_l = pi - asin(sin(FA) + ripple_M);
        end  
        
        % calculate ripple of |Bn(z)|
        [range_B, ripple_B] = rf_ripple_FA2Beta(rfa_r, rfa_l, FA);
        
    case {'sat','inv'}
        if ( (cos(FA)+ripple_M) > 1 )
            rfa = acos(cos(FA) - ripple_M);
            rfa_l = -rfa;
            rfa_r = rfa;
        elseif ( (cos(FA)-ripple_M) < -1 )
            rfa = acos(cos(FA) + ripple_M);
            rfa_l = rfa;
            rfa_r = 2*pi-rfa;
        else
            rfa_r = acos(cos(FA) - ripple_M);
            rfa_l = acos(cos(FA) + ripple_M);
        end
        
        % calculate ripple of |Bn(z)|
        [range_B, ripple_B] = rf_ripple_FA2Beta(rfa_r, rfa_l, FA);
        
    case 'se'
        mid_M = (sin(FA/2))^2;
        max_M = max(0,min(1,mid_M+ripple_M));
        min_M = max(0,min(1,mid_M-ripple_M));        
        max_B = sqrt(max_M);
        min_B = sqrt(min_M);
        mid_B = sin(FA/2);
        range_B = [min_B, max_B];
        ripple_B = abs( range_B - mid_B );
    
    otherwise
        disp(['Unrecognized Pulse Type -- ',ptype]);
        disp('Recognized types are st, ex, se, inv, and sat');
        return;
end


end


function [range_B, ripple_B] = rf_ripple_FA2Beta(rfa_r, rfa_l, FA)
%
% calculate the range of Bn(z) polynomial given the range of Flip-angle.
% |Bn(z)| = sin(theta/2) is an increasing function for theta in [0 pi),
% except around theta = pi. 
%
% Inputs:
%    rfa_r -- in radians, right bound of flip-angle
%    rfa_l -- in radians, left bound of flip-angle
%    FA -- in radians, desired flip-angle
%
% Outputs:
%    range_B -- range of magnitude of Beta profile
%    ripple_B -- ripple of magnitude of Beta profile
%
%  written by Hong Shang, Sep 2014
%  JGGB UCSF/UC Berkeley

if (rfa_r) > pi
    min_B = min([ sin(rfa_l/2), sin(rfa_r/2) ]);
    max_B = 1;
    range_B = [min_B, max_B];
    ripple_B = 1 - min_B;

else
    min_B = sin(rfa_l/2);
    max_B = sin(rfa_r/2);
    mid_B = sin(FA/2);
    range_B = [min_B, max_B];
    ripple_B = abs( range_B - mid_B );
    
end


end

