function y = myUnwrap(x)
%
%  y = myUnwrap(x)
%
%  Unwrap phase angle. Different from matlab function unwrap(angle()), 
%  here considering zero-crossing point to get rid of pi jump
%  Not optimized in terms of speed and threshold value
%
%  Inputs:
%    x -- N*1 vector, complex value 
%
%  Outputs:
%    y -- N*1 vector, phase in radians
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Jun 2015

thre_pha = pi*0.9;
thre_mag = 1e-3;

% basic phase unwrap
x = x(:);
pha = unwrap(angle(x));

% compensate phase jump of pi where amplitude is close to zero
for i=2:length(pha)
    if ( ((pha(i) - pha(i-1)) > thre_pha) && (abs(x(i)) < thre_mag) )
        pha(i:end) = pha(i:end) - pi;
    elseif ( ((pha(i) - pha(i-1)) < (-thre_pha)) && (abs(x(i)) < thre_mag) )
        pha(i:end) = pha(i:end) + pi;
    end
end

y = pha;
