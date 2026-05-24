function [va]=tas2va(tas,alpha,beta)
%TAS2VA:  Transform airspeed to 3 components
%function [va]=tas2va(tas,alpha,beta)
%$Source: /home/cvs/kingair/Sys09/tas2va.m,v $
%Project $Name: trans2am21_qc0 $ ($Revision: 1.1 $)
%$Date: 2010/08/10 09:07:15 $
%

% airspeed in aircraft coordinates [angles in radians]
if(nargin == 3)
        ta=tan(alpha);
        tb=tan(beta);
        tasab=tas./sqrt(1.+ta .^2+tb .^2);
        va=[tasab, tasab .* tb,tasab .* ta];
else
'%function [va]=tas2va(tas,alpha,beta) [alpha and beta in radians]'
	error('tas2tas: 3 arguments needed')
end
