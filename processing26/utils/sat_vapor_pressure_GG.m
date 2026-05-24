function e=sat_vapor_pressure_GG(Tdp) 
%
% Input: Temperature [K]
% Output: Saturation vapor pressure [millibar]
%

Ts = 373.15;    % [K] standard temperature at steam point on Kelvin scale (Goff, 1965)
Ps = 1013.246;  % standard atmospheric pressure at steam point (hPa)
log10e = -7.90298*(Ts./Tdp - 1) ...
       + 5.02808*log10(Ts./Tdp) ...
       - 1.3816e-7*(10.^(11.344*(1 - Tdp./Ts)) - 1) ...
       + 8.1328e-3*(10.^(-3.49149*(Ts./Tdp - 1)) - 1) ...
       + log10(Ps);

e = 10.^log10e; % hPa
% This is sensitive to numerics
kk = find(~isnan(e) & ~isinf(e));
if ~isempty(kk)
    e = interp1(kk,e(kk),[1:numel(e)]','linear',0);
end

return
