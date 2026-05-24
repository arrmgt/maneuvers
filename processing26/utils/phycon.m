function C=phycon(varargin);
%Physical constants (from Iribarne and Godson (1981) second ed):
% returns
%C.Mv=18.015;
%C.Md=28.964;
%C.Rd=287.05;
%C.Rv=461.51;
%C.Cpd=1005;
%C.Cvd=718;
%C.Cpv=1850;
%C.Cvv=1390;
%C.Cw=4218;
%C.eps=C.Mv/C.Md;
%C.Tzero=273.15;
%C.pstp=1013.25;
%C.Tstp=288.15;
%C.g0=9.80665;

C.Mv=18.015;  % Mol mass of water
C.Md=28.964;  % Mol mass of dry air
C.Mco2=44.01; % Mol mass CO2
C.Rd=287.05;  % gas constant dry air
C.Rv=461.51;  % gas constant water vapor
C.Cpd=1005;   % Spec heat const press for dry air
C.Cvd=718;    % Spec heat const vol for dry air
C.Cpv=1850;   % Spec heat const press for water vapor
C.Cvv=1390;   % Spec heat const vol for water vapor
C.Cw=4218;    % Spec heat for liquid water
C.eps=C.Mv/C.Md;
C.Tzero=273.15;
C.pstp=1013.25;
C.Tstp=288.15;
C.g0=9.80665;
C.Rstar=8.314;% universal gas constant (J/mole)
