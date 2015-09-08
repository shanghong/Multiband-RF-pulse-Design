function plot_spec(f, a, d, name_array)
% PLOT_SPEC - Utility to plot frequency specifications
%   
%  plot_spec(f, a, d, type)
%
%  f - frequency band edges
%  a - band amplitudes
%  d - band ripple specs
%  type - line/plotting type (see S in 'help plot')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectral-Spatial RF Pulse Design for MRI and MRSI MATLAB Package
%
% Authors: Adam B. Kerr and Peder E. Z. Larson
%
% (c)2007-2011 Board of Trustees, Leland Stanford Junior University and
%	The Regents of the University of California. 
% All Rights Reserved.
%
% Please see the Copyright_Information and README files included with this
% package.  All works derived from this package must be properly cited.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% modified by Hong Shang, July 2015, to display spec as well as name


type = 'k-';
type_lw = 1.7;
nband = length(f)/2;

if nargin == 3
    for band = 1:nband,
        idx = [band*2-1:band*2];
        hp = plot(f(idx), a(idx)+d(band)*ones(1,2), sprintf('%s',type));
        set(hp,'linewidth',type_lw);
        hold on;
        hp = plot(f(idx), a(idx)-d(band)*ones(1,2), sprintf('%s',type));
        set(hp,'linewidth',type_lw);
    end;
    
elseif nargin == 4
    text_pos = 1;
    for band = 1:nband,
        idx = [band*2-1:band*2];
        hp = plot(f(idx), a(idx)+d(band)*ones(1,2), sprintf('%s',type));
        set(hp,'linewidth',type_lw);
        hold on;
        hp = plot(f(idx), a(idx)-d(band)*ones(1,2), sprintf('%s',type));
        set(hp,'linewidth',type_lw);

        text_pos = min([text_pos mean(a(idx))+d(band)]);
    end

    text_pos = text_pos + 0.1;
    for band = 1:nband
        idx = [band*2-1:band*2];
        text(mean(f(idx)), text_pos, name_array{band}, 'FontSize',18, 'rotation',90);
    end

end

