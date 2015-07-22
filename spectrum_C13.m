function [f,name_cell] = spectrum_C13(B0,dbg)
%
%    [f,name_cell] = spectrum_C13(B0,dbg)
%
%    calculate the chemical shift frequency of different 13C labeled substance,
%    including pyruvate, lactate, alanine, pyruvate hydrate, HCO3, urea
%    
%    Inputs:
%    B0: main magnetic filed, in Tesla, 3 or 7 typically
%    dbg: 1 show debug plot, 0 no
%    
%    Outputs:
%    plot of spectrum
%    f, in Hz, resonnace frequency shift
%    name_cell, cell array of name, corresponding to f 
%   
%  (c) 2013-2015 The Regents of the University of California
%  All Rights Reserved.
%  Author: Hong Shang  May 2013
%
%  Modified chemical shift of pyr/lac/urea/ala based on 14T measurements
%  by Hong Shang Feb 2015

if nargin == 1
    dbg = 1;
end;

gamma = 10.705*1e6;  % gamma/2pi for C13, in Hz/T
cs_pyr = 170.60;     % chemical shift for pyruvate, in ppm
cs_lac = 182.98;     % chemical shift for lactate, in ppm
cs_ala = 176.32;     % chemical shift for alanine, in ppm
cs_pyr_H2O = 178.91; % chemical shift for pyruvate hydrate, in ppm
cs_bicarb = 160.9;   % chemical shift for bicarbonate HCO3, in ppm
cs_urea = 163.13;    % chemical shift for urea, in ppm
cs = [cs_pyr,cs_lac,cs_ala,cs_pyr_H2O,cs_bicarb,cs_urea];
name_cell = {'Pyruvate' 'Lactate' 'Alanine' 'Pyruvate-H_2O' 'Bicarbonate' 'Urea'};

f0 = gamma*B0*(1+cs*1e-6); % resonance frequency, in Hz
fr = f0(1);                % reference frequency, in Hz
f = f0-fr;                 % resonance frequency compared to the reference, in Hz

if dbg >= 1 
    figure; stem(f,ones(size(f))); text(f,[1.1 1.1 1.1 0.9 1.1 0.9].*ones(size(f)),name_cell,'FontSize',17,'HorizontalAlignment','center'); 
    xlabel('Relative frequency (Hz)','FontSize',17); title('C-13 NMR spectrum','FontSize',20); set(gca,'ylim',[0, 1.4]); set(gca,'xlim',[min(f)*1.2,max(f)*1.15]);
    leg = legend(['at ',num2str(B0),' T']); set(leg, 'FontSize', 18);  set(gca,'FontSize',17); set(gca,'ytick',[]);
    set(gcf, 'Position', [400, 50, 800, 400], 'PaperPositionMode', 'auto');
    print(gcf,'-depsc2','spectrum_C13.eps');
   
end
