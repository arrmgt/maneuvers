function [zz,dz]=ZHydro(z0,pmb,TK,varargin)
%ZHYPSO: hypsometric altitude
%function zz=ztrue(z0,pmb,TK,[mr,[grav]]) 
%
%Inputs:
%  z0       Initial height [m]
%  pmb      Static pressure [any units ok]
%  TK       Temperature [K]
%  [mr]     Mixing Ratio [g/g] (optional)
%  [grav]   Gravity [m/s2] (optional)
%
%Outputs:
%
%  ztrue    High from integration of hypsometric equation
%  dz       Differential of height (for debugging)
%

C=phycon(); % get physical constants
g = C.g0;
p = inputParser;
addParameter(p, 'Mixing_Ratio', 0, @(x) isnumeric(x));
addParameter(p, 'Gravity'     , g, @(x) isnumeric(x));
parse(p, varargin{:});
Opts = p.Results;

g = Opts.Gravity;
mr = Opts.Mixing_Ratio;

pPa = pmb*100;
rho = pPa./(C.Rd.*TK);

% Compute density (ideal gas)
rho = pPa ./ (C.Rd .* TK);

% Integrate dz = - dp ./ (rho .* g) with respect to pressure using cumtrapz.
% If pPa is decreasing with index (e.g., top->bottom), integrate accordingly.
% We'll compute z relative to first element pPa(1):
dz_cumulative = cumtrapz(pPa, -1 ./ (rho .* g));   % integral from pPa(1) to pPa(i)
zz = z0 + dz_cumulative;

end
