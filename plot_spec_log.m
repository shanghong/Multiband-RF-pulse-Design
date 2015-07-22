function plot_spec_log(f, a, d, name_array)
%   
%  plot_spec_log(f, a, d, name_array)
%  
%  Utility to plot frequency specifications with logarithmic scale
%
%  f - frequency band edges
%  a - band amplitudes
%  d - band ripple specs
%  name_array - name of each band [optional]
%
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  June 2014
%
% modified by Hong Shang, June 2015, to display spec as well as name

type = 'k-';
type_lw = 2;
nband = length(f)/2;

if (nargin == 3)
    
    for band = 1:nband,
        idx = [band*2-1:band*2];
        hp = semilogy(f(idx), ((a(idx)+d(band)*ones(1,2))>0) .* (a(idx)+d(band)*ones(1,2)), sprintf('%s',type));
        set(hp,'linewidth',type_lw);
        hold on;
        hp = semilogy(f(idx), ((a(idx)-d(band)*ones(1,2))>0) .* (a(idx)-d(band)*ones(1,2)), sprintf('%s',type));
        set(hp,'linewidth',type_lw);
    end;
    
elseif (nargin == 4)
    
    text_pos = 1;
    for band = 1:nband,
        idx = [band*2-1:band*2];
        hp = semilogy(f(idx), ((a(idx)+d(band)*ones(1,2))>0) .* (a(idx)+d(band)*ones(1,2)), sprintf('%s',type));
        set(hp,'linewidth',type_lw);
        hold on;
        hp = semilogy(f(idx), ((a(idx)-d(band)*ones(1,2))>0) .* (a(idx)-d(band)*ones(1,2)), sprintf('%s',type));
        set(hp,'linewidth',type_lw);

        text_pos = min([text_pos mean(a(idx))+d(band)]);
    end

    text_pos = text_pos*1.3;
    for band = 1:nband
        idx = [band*2-1:band*2];
        text(mean(f(idx)), text_pos, name_array{band}, 'FontSize',17, 'rotation',90);
    end

end
