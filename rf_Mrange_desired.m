function range_M = rf_Mrange_desired(FA, ripple_M, ptype)
%
%  For SLR algorithm with generalized flip-angle, calculate desired range 
%  of amplitude of magnetization
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
%
%  Outputs:
%    range_M -- range of amplitude of magnetization
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Sep 2014

range_M = [];
% check whether FA in the range
if ( (FA < 0) || (FA > 180) )
   disp(['Flip angle of ', num2str(FA), ' degree is out of range']);
   disp('Flip angle should be in the range of [0 180] degree');
   return;
end
FA = FA*pi/180;  % in radians

switch ptype
    case {'st','ex'} 
        range_M = [ sin(FA)-ripple_M, min([sin(FA)+ripple_M, 1]) ];
          
    case {'sat','inv'}
        range_M = [ max([cos(FA)-ripple_M, -1]), min([cos(FA)+ripple_M, 1]) ];
             
    case 'se'
        range_M = [ max([(sin(FA/2))^2-ripple_M, 0]), min([(sin(FA/2))^2+ripple_M, 1]) ];
           
    otherwise
        disp(['Unrecognized Pulse Type -- ',ptype]);
        disp('Recognized types are st, ex, se, inv, and sat');
        return;
end

