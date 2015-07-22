% function [aca] = b2a(bc);
%
% This function takes a b polynomial, and returns the minimum phase, 
%   minimum power a polynomial
%
% Inputs:
%   bc - beta polynomial coefficients
%
% Outputs:
%   aca - minimum phase alpha polynomial
%

function [aca] = b2a(bc)

n = length(bc);
% calculate minimum phase alpha
bcp = bc;
bl = length(bc);
blp = bl*8;
bcp(blp) = 0;
bf = fft(bcp);
bfmax = max(abs(bf));
if bfmax>=1.0,                % PM can result in abs(beta)>1, not physical
  bf = bf/(1e-8 + bfmax);  %   we scale it so that abs(beta)<1 so that
end;                          %   alpha will be analytic
afa = mag2mp(sqrt(1-bf.*conj(bf)));
aca = fft(afa)/blp;
aca = aca(n:-1:1);


