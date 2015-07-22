function rfs = rfscaleg(rf,t,gamma)

%  rfs = rfscaleg(rf,t,gamma)
%
%    rf  -- rf waveform, scaled so sum(rf) = flip angle
%    t   -- duration of RF pulse in ms
%    rfs -- rf waveform scaled to Gauss
%
%  written by John Pauly, 1992
%  (c) Board of Trustees, Leland Stanford Junior University


gamma = 2*pi*gamma; % kHz/G
dt = t/length(rf);
rfs = rf/(gamma*dt);

