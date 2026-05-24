function [minterp,mdecim]=interp_decim(irate,orate);
% Determine interpolation and decimation factors
%
%  Inputs
%   irate :  Input rate
%   orate :  Output rate

if(irate == orate),
	minterp=1;mdecim=1;
        return
end

% Use matlab "rat" function
%   [N,D] = rat(___) returns two factors
%   N and D, such that N./D approximates X 
r = irate / orate;
[mdecim, minterp] = rat(r,1e-12);

end
