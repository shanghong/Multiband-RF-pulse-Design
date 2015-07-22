function rfwrite(rf, nompw, ang, GAMMA, isodelay, g, thk)
% rfwrite(rf, nompw, ang, GAMMA, isodelay, nombw)
%  OR
% rfwrite(rf, nompw, ang, GAMMA, isodelay, g, thk)
% 
% Write RF pulse and rfstat information into prompted filename
% optionally will write phase and gradient
% Files can be read back into matlab by rfread() or signaread()
%   
% rf - in G
% nompw - in sec
% ang - flip angle in radians
% GAMMA (optional) - in Hz/G (leave empty for protons)
% isodelay (optional) - in sec, leave empty for estimate
% nombw - in Hz (can be set to zero if unnecessary)
%   OR
% g - in G/cm
% thk - thickness in cm
%    
% PEZL 9/26/07, 6/16/08
% based on code from Chuck and Adam
%
% modified by Hong Shang July 2014
% reorder parameters to be consistent with our system


GAMMA_H1 = 4257;
if (nargin < 4) || isempty(GAMMA)
  GAMMA = GAMMA_H1;
end

if nargin > 6
  gwrite = 1;
else
  gwrite = 0;
end

root_fname = input('Root file name: ', 's');
if isempty(root_fname)
  fprintf(1,'Not saving files \n');
  return;
end;

% scale for proton GAMMA for epic
rf = rf * GAMMA / GAMMA_H1;
maxrf = max(abs(rf));
rfn = rf / maxrf;
nrf = length(rf);

% write dat file    
dat_name = sprintf('%s.dat', root_fname);
fid = fopen(dat_name, 'w');
if fid == -1, 
  fprintf(1, 'Error opening %s \n', dat_name);
  return;
end;
    

pon = (rfn >= 0.00001);
temp_pw = 0;
max_pw = 0;
for n=1:nrf
  temp_pw = temp_pw + pon(n);
  if (and(pon(n) == 0, temp_pw ~= 0))
    max_pw = max(max_pw, temp_pw);
    temp_pw = 0;
  end;
end;
max_pw = max_pw / n;
    
dty_cyc = sum(abs(rfn) > 0.2236)/nrf; 
if dty_cyc < max_pw, 
  dty_cyc = max_pw;
end;

   
if nargin < 6
    nombw = 0.0;
elseif gwrite
    maxg = max(abs(g));
    nombw = GAMMA * maxg * thk;
else
    nombw = g;
end



fprintf(fid,'%10d \t\t #extgradfile\n', gwrite);
fprintf(fid,'%10d \t\t #res\n', nrf);
fprintf(fid,'%10d \t\t #pw\n',round(nompw*1e6));
fprintf(fid,'%10.7f \t\t #nom_flip \n',ang*180/pi);

abswidth = sum(abs(rfn))/nrf;
fprintf(fid,'%10.7f \t\t #abswidth \n',abswidth);

effwidth = sum(abs(rfn).^2)/nrf;
fprintf(fid,'%10.7f \t\t #effwidth \n',effwidth);

area = sum(abs(rfn))/nrf;
fprintf(fid,'%10.7f \t\t #area \n',area);

fprintf(fid,'%10.7f \t\t #dtycyc \n',dty_cyc);
fprintf(fid,'%10.7f \t\t #maxpw \n',max_pw);

max_b1 = maxrf;
fprintf(fid,'%10.7f \t\t #max_b1 \n',max_b1);
 
int_b1_sqr = sum(abs(rf).^2 * nompw / nrf * 1e3);
fprintf(fid,'%10.7f \t\t #max_int_b1_sqr \n',int_b1_sqr);
    
rms_b1 = sqrt(sum(abs(rf).^2))/nrf;
fprintf(fid,'%10.7f \t\t #max_rms_b1 \n',rms_b1);

fprintf(fid,'%10.7f \t\t #nom_bw \n',nombw);


if gwrite     
  fprintf(fid,'%10.3f \t\t #a_gzs \n',maxg);
    
  thk_scale = thk * GAMMA / GAMMA_H1 * 10;
  fprintf(fid,'%10.3f \t\t #nom_thk(mm) \n',thk_scale);
end

fclose(fid);
    
% Now write out RF and Gradient 
%
rho_fname = sprintf('%s.rho', root_fname);
if any(imag(rf(:)))
    signa(abs(rfn),rho_fname);

  theta_fname = sprintf('%s.pha', root_fname);
  signa(angle(rfn),theta_fname,1/pi);
else
    signa(rfn,rho_fname);    
end

if gwrite
  g_fname = sprintf('%s.grd', root_fname);
  signa(g,g_fname);
end
