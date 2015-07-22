function rfwrite_varian(rf, nompw, nombw, root_fname)
%
% rfwrite_varian(rf, nompw, root_fname)
% 
% Write a RF pulse into a file with the format of varian system
% the parameters in the header are not exactly sure
%
% Inputs:
%    rf -- unitless, can be complex
%    nompw -- in ms, nominal duration
%    nombw -- in kHz, nominal bandwidth
%    root_fname -- name of file
% 
% Output
%   root_fname.RF
%    
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  Oct 2014

if ( (nargin <= 2) || isempty(nombw) )
    nombw = 0.5;  % in kHz, nominal bandwidth, which will not be used actually
end

if nargin <= 3
    root_fname = input('Root file name: ', 's');
    if isempty(root_fname)
        fprintf(1,'Not saving files \n');
        return;
    end;
end

% write dat file    
dat_name = sprintf('%s.RF', root_fname);
fid = fopen(dat_name, 'w');
if fid == -1, 
  fprintf(1, 'Error opening %s \n', dat_name);
  return;
end;


rf_pha = angle(rf)*180/pi;          % in degree
rf_mag = 1024*abs(rf)/max(abs(rf)); % normalized to 1024
rf_int = sum(rf_mag)/(1024*length(rf_mag));

fprintf(fid,'# VERSION   100\n');           % this is not sure
fprintf(fid,'# TYPE    selective\n');
fprintf(fid,'# MODULATION  amplitude\n');
fprintf(fid,'# EXCITEWIDTH   1.8125\n');    % this is not sure
fprintf(fid,'# INVERTWIDTH   0\n');        
fprintf(fid,'# INTEGRAL   %1.5f\n',rf_int);
fprintf(fid,'# T(ms)xBW(KHz)   %1.4fx%2.1f\n', nompw, nombw); 


for i=1:length(rf)
    fprintf(fid,'%3.7f \t %4.7f \t %2.7f\n', rf_pha(i), rf_mag(i), 1);
end

fclose(fid);
    
