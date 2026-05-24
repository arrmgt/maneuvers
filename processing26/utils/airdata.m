function OUT = airdata(Ps_meas, Pt_meas, Tm, recovf, Td, varargin)
% AIRDATA
% Compute moist-air thermodynamic and airspeed quantities from aircraft
% measurements, using MKS units, compressible flow, and moist
% thermodynamics.
% Includes pressure correction, humidity effects, true airspeed, Mach,
% virtual/density altitude
%
% -------------------------------------------------------------------------
% SYNTAX:
%
%   OUT = airdata_moist(Ps_meas, Pt, Tm, r, Td)
%   OUT = airdata_moist(..., 'dPs_corr', dPs_corr)
%
% -------------------------------------------------------------------------
% REQUIRED INPUTS  (arrays or scalars, same size)
%
%   Ps_meas : Measured static pressure from port (hPa) 
%   Pt      : Measured Pitot / total pressure (hPa)
%   Tm      : Measured probe (sensor) temperature (K) 
%   recovf  : Probe recovery factor (dimensionless, required)
%   Td      : Dewpoint temperature (K)
%
% -------------------------------------------------------------------------
% OPTIONAL NAME-VALUE INPUTS
%
%   'dPs_corr' : Static pressure correction (hPa), scalar or array
%                Ps_corr = Ps_meas + dPs_corr is used in all calculations
%                If not entered, inputs assumed to be already corrected.
%   'Z_gps'    : GPS altitude (m) used to initialize Zhydrostatic calc.
%
% -------------------------------------------------------------------------
% OUTPUT:  Structure OUT containing the following fields
%
%  --- PRESSURES ---
%   Ps_meas    : Measured static pressure (hPa)
%   dPs_corr   : Applied pressure correction (hPa)
%   Ps_corr    : Corrected static pressure used in calculations (hPa)
%   Pt         : Total (Pitot) pressure (hPa)
%   q_impact   : Impact pressure = Pt - Ps_corr (hPa)
%   q_dyn      : Dynamic pressure = 0.5 * rho * Vt^2 (hPa)
%
%  --- TEMPERATURES ---
%   Tm         : Measured temperature (K)
%   Ts         : Static (ambient) temperature (K)
%   Tt         : Stagnation (total) temperature (K)
%   Tv         : Virtual temperature (K)
%   Td         : Dewpoint temperature (K)
%
%  --- HUMIDITY / MOISTURE ---
%   e          : Actual vapor pressure (hPa)
%   es         : Saturation vapor pressure at Ts (hPa)
%   w          : Mixing ratio (kg vapor / kg dry air)
%   q          : Specific humidity (kg vapor / kg moist air)
%   RH         : Relative humidity (0–1)
%
%  --- THERMODYNAMIC PROPERTIES (moist air) ---
%   Cp         : Specific heat at constant pressure (J/kg/K)
%   Cv         : Specific heat at constant volume (J/kg/K)
%   Rm         : Gas constant of moist air (J/kg/K)
%   gamma      : Ratio of specific heats, Cp/Cv (–)
%
%  --- DENSITY & SOUND SPEED ---
%   rho        : Air density (kg/m^3)
%
%  --- AIRSPEEDS & MACH ---
%   Mach       : Mach number (–)
%   Vt         : True airspeed (TAS) (m/s)
%   Vc         : Calibrated airspeed (CAS) (m/s)
%   Vi         : Indicated airspeed (IAS ≈ Vc) (m/s)
%
%  --- ALTITUDES (Standard & Derived) ---
%   Zp         : Standard pressure altitude (0–11 km model) [m]
%   Zp_ft      : Pressure altitude [ft]
%   Zp_trop    : Tropopause-aware pressure altitude (to 20 km) [m]
%   Zhydro     : Geometric altitude via hysrostatic equation [m]
%   Zgph       : Geopotential height [m]
%   Zvirt      : Virtual-temperature-based height [m]
%   Zdens      : Density altitude [m]
%   Zdens_ft   : Density altitude [ft]
%
%
% -------------------------------------------------------------------------
% ASSUMPTIONS:
%  • Compressible flow (subsonic, no shock)
%  • Isentropic deceleration at Pitot probe
%  • Moist thermodynamics using WMO/ITS-90 saturation vapor pressure
%  • Static pressure correction deterministic (not part of σ propagation)
%  • All inputs in MKS units
%
% -------------------------------------------------------------------------
% EXAMPLE:
%   OUT = airdata_moist(Ps, Pt, Tm, 0.985, Td, ...
%                       'dPs_corr', 45, 'Z_gps', 1340, ...
%                       'sigmaPs', 15, 'sigmaTd',0.5);
%
%   TAS = OUT.TAS;   Mach = OUT.Mach;  Zp_ft = OUT.Zp_ft;
%
% -------------------------------------------------------------------------
%   Corrected static pressure :
%       Ps_corr = Ps_meas + dPs_corr
%   Optional: name/value input:
%   'dPs_corr' : Static pressure correction (hPa), scalar or array
%
% All humidity variables ar computed from Td 
% Goff-Gratch saturation vapor pressure formulation is used
%
% OUTPUT includes corrected Ps, Ts, Tt, Tv, Mach, TAS, CAS, IAS, Cp, Cv,
% Rm, gamma, density, and humidity variables.
%
% -------------------------------------------------------------------------

% ---------- Parse static pressure error (hPa) ----------
p = inputParser;
addParameter(p, 'dPs_corr', [], @(x)isnumeric(x)); 
addParameter(p, 'Z_gps', 0, @(x)isnumeric(x)); 
parse(p, varargin{:});
Opts = p.Results;

% ---------- Apply static pressure correction ----------
% ---------- And check input data ----------------------
kk0 = 1:numel(Ps_meas);
dPs_corr = Opts.dPs_corr;
if ~isempty(dPs_corr)
    kk = find (~isnan(dPs_corr) & ~isinf(dPs_corr) ...
        & abs(dPs_corr)<15 & abs(gradient(dPs_corr))<10 );
    dPs_corr = interp1(kk,Opts.dPs_corr(kk),kk0',"linear",0);
else
    dPs_corr = zeros(size(Ps_meas));
end

kk = find (~isnan(Ps_meas) & ~isinf(Ps_meas) ...
    & Ps_meas>100 & Ps_meas<1200 & abs(gradient(Ps_meas))<10 );
Ps_meas = interp1(kk,Ps_meas(kk),kk0',"linear",500);

kk = find (~isnan(Pt_meas) & ~isinf(Pt_meas) ...
    & Pt_meas>100 & Pt_meas./Ps_meas>1 & abs(gradient(Pt_meas))<10);
Pt_meas = interp1(kk,Pt_meas(kk),kk0',"linear",550);

kk = find (~isnan(Tm) & ~isinf(Tm) & Tm >200);
Tm = interp1(kk,Tm(kk),kk0',"linear",200);

kk = find (~isnan(Td) & ~isinf(Td) & Td>200 );
Td = interp1(kk,Td(kk),kk0',"linear",200);

% ---------- Input size checks ----------
if ~isequal(size(Ps_meas), size(Pt_meas), size(Tm), size(Td))
    error('Ps, Pt, Tm, and Td must be same size.');
end

% ---------- Constants ----------

C=phycon();
Mv=C.Mv; Md=C.Md;
Rd=C.Rd; Rv=C.Rv;
Cpd=C.Cpd; Cpv=C.Cpv;
Cvd=C.Cvd; Cvv=C.Cvv;
eps=C.eps;
eps = Rd ./ Rv;
gamma = Cpd./Cvd;
% Standard temperature and pressure
P0   = 101325.0; T0 = 288.15; gamma0 = 1.4; R0 = Rd;
P0_hPa = P0/100;
a0   = sqrt(gamma0 * R0 * T0); % speed of sound

% ---------- Run physics ----------

BASE = core_calc(Ps_meas, Pt_meas, Tm, recovf, Td, ...
                 Cpd, Cvd, Rd, Cpv, Cvv, Rv, eps, ...
                 P0_hPa, gamma0, R0, a0, ...
                 "dPs_corr", dPs_corr, "Z_gps", Opts.Z_gps);
OUT = BASE;
OUT.Ps_meas = Ps_meas;   % Keep both
OUT.Ps_corr =  dPs_corr;
OUT.dPs_corr = Opts.dPs_corr;
OUT.DZhydro = OUT.Zhydro - Opts.Z_gps;

end 

% =========================================================================
%   CORE CALCULATION (fully MKS)
% =========================================================================
function S = core_calc(Ps_meas, Pt_meas, Tm, recovf, Td, ...
                       Cpd, Cvd, Rd, Cpv, Cvv, Rv, eps, ...
                       P0_hPa, gamma0, R0, a0, varargin)
%
% Ps_meas and Pt_meas are measured pressures (hPa)
% Corrections applied:
%   Ps = Ps_meas + dPs_corr;
%   Pt = Pt_meas + dPs_corr;
% All internal pressure calculations use Ps and Pt (hPa)

%% --- Name-value pressure corrections ---------------------------------
p = inputParser;
addParameter(p, 'dPs_corr', 0, @(x) isnumeric(x));
addParameter(p, 'Z_gps', 0, @(x) isnumeric(x));
parse(p, varargin{:});

dPs_corr = p.Results.dPs_corr;   % correction for static pressure
Zgps = p.Results.Z_gps; % GPS altitude

% Apply corrections
Ps = Ps_meas - dPs_corr;   % corrected static pressure [hPa]
Pt = Pt_meas; % No static correction needed

%% --- HUMIDITY (all pressures are hPa) ---------------------------------
e  = sat_vapor_pressure_GG(Td);    % hPa
w  = eps .* e ./ (Ps - e);          % mixing ratio [kg/kgDryAir]
q  = w ./ (1 + w);                  % specific humidity [kg/kgAir]
% Recompute full moist thermodynamic properties
    Yd = 1./(1+w);   Yv = w./(1+w);
    Cp = Yd.*Cpd + Yv.*Cpv;
    Cv = Yd.*Cvd + Yv.*Cvv;
    Rm = Yd.*Rd  + Yv.*Rv;
    gamma = Cp ./ Cv;
% Pressure Ratio
    PR = Pt ./ Ps; % Ptot/Ps
    kk = find(PR>1);
    PR = interp1(kk,PR(kk),[1:numel(PR)]','linear',1);
% Temperature calcs
    k = (gamma - 1)./gamma; % k = R/Cp
    M = sqrt((2./(gamma-1)).*(PR.^k - 1));
    % Static air temperature
        %Ts = Tm ./ (1 + recovf*((gamma - 1)/2).*M.^2);
    Ts = Ts_from_Tm(Tm, recovf, Rd, Rv, "gamma",gamma,'mach',M);
    % Stagnation temperature (r=1)
    Tt = Ts .* ( 1 + ((gamma-1)./2) .* M.^2 );
    % Relative humidity
    es = sat_vapor_pressure_GG(Ts);
    RH = e./es;

%% --- Virtual Temperature ---------------------------------------------
Tv = Ts .* (1 + w./eps) ./ (1 + w);

%% --- DENSITY (Pa conversion only here) -------------------------------
rho = (Ps .* 100) ./ (Rm .* Ts);

%% --- AIRSPEEDS --------------------------------------------------------
%   Mach       : Mach number (–)
%   Vt         : True airspeed (TAS) (m/s)
%   Vc         : Calibrated airspeed (CAS) (m/s)
%   Vi         : Indicated airspeed (IAS ≈ Vc) (m/s)
a = sqrt(gamma .* Rm .* Ts);
Vt = M .* a;

x = Pt - Ps;                      % impact pressure [hPa] (corrected)
    x(find(x<=0)) = 10; %sanity check
q_impact = x;
q_dyn = (0.5 .* rho .* Vt.^2)./100;       % dynamic pressure [hPa]

PR0 = 1 + q_impact ./ P0_hPa;
M0  = sqrt((2./(gamma0-1)).*(PR0.^((gamma0-1)/gamma0) - 1));
Vc  = M0 .* a0;
Vi  = Vc;

%% --- POTENTIAL TEMPERATURES -----------------------------------------
Cp_d = 1005;    Lv = 2.5e6;
theta   = Ts .* (P0_hPa./Ps).^(Rd/Cp_d);
theta_d = Td .* (P0_hPa./Ps).^(Rd/Cp_d);
theta_v = theta .* (1 + 0.61*q);
theta_e = theta .* exp((Lv.*q) ./ (Cp_d.*Td));

%% --- LCL Temperature and Pressure (Bolton)   -------------------------
TLCL = 1./(1./(Td-56) + log(Ts./Td)/800) + 56;
PLCL = Ps .* (TLCL ./ Ts).^(Cp_d/Rd);    % hPa

%% --- ALTITUDE CALCULATIONS -------------------------------------------
g  = 9.80665;
T0 = 288.15;
L  = 0.0065;
rho0 = 1.225;

% Standard Atmosphere pressure altitude
Zp = (T0./L) .* (1 - (Ps./P0_hPa).^((Rd*L)/g));

% Standard Atmosphere with tropopause
Zp_trop = zeros(size(Ps));
for i = 1:numel(Ps)
    if Ps(i) >= 226.32
        Zp_trop(i) = (T0/L)*(1 - (Ps(i)/P0_hPa).^((Rd*L)/g));
    else
        Zp_trop(i) = 11000 + (T0 - L*11000)/g * log(226.32/Ps(i));
    end
end

% Hydrostatic altitude (moist) -- integral of hydrostatic equation
    Zhydro = ZHydrostatic(Zgps(1), Ps, Ts, 'Mixing_Ratio', w);
% Virtual altitude
    Zvirt = (Rd .* Tv ./ g) .* log(P0_hPa ./ Ps);
% Geopotential height
    Zgph  = Zhydro;

% Density Altitude
    Zdens = (1 - (rho./rho0).^(1/0.2349)) .* (T0./L);
    Zp_ft = Zp * 3.28084;
    Zdens_ft = Zdens * 3.28084;

%% --- PACK OUTPUT ------------------------------------------------------
S.Ps_meas = Ps_meas;
S.Pt_meas = Pt_meas;
S.dPs_corr = dPs_corr;
S.Ps = Ps;
S.Pt = Pt;
S.PR = PR;
S.PLCL = PLCL;

% Impact and Dynamic pressures
S.q_impact = q_impact;
S.q_dyn = q_dyn;    % Pa

% Air density
S.rho   = rho;

% Airspeeds
S.M  = M;
S.TAS = Vt;  S.Vc = Vc;  S.Vi = Vi;

% Temperatures
S.Ts = Ts;   S.Tt = Tt;   S.Tv = Tv;   S.Td = Td;   S.TLCL = TLCL;

% Humidity
S.w = w;     S.q = q;     S.RH = RH;   S.e = e;   S.es = es;

% Gas properties
S.Cp = Cp;   S.Cv = Cv;   S.Rm = Rm;   S.gamma = gamma;

% Potential temperatures
S.theta   = theta;
S.thetad  = theta_d;
S.theta_v = theta_v;
S.theta_e = theta_e;

% Altitudes
S.Zp       = Zp;
S.Zp_ft    = Zp_ft;
S.Zp_trop  = Zp_trop;
S.Zhydro    = Zhydro;
S.Zvirt    = Zvirt;
S.Zgph     = Zgph;
S.Zdens    = Zdens;
S.Zdens_ft = Zdens_ft;

%% --- UNITS TABLE ------------------------------------------------------
S.units.Ps_meas  = 'hPa';
S.units.Pt_meas  = 'hPa';
S.units.dPs_corr = 'hPa';
S.units.Ps       = 'hPa';
S.units.Pt       = 'hPa';
S.units.q_impact = 'hPa';
S.units.PLCL     = 'hPa';
S.units.q_dyn    = 'hPa';
S.units.rho      = 'kg/m^3';
S.units.PR       = '1';

S.units.Ts       = 'K';
S.units.Tt       = 'K';
S.units.Tv       = 'K';
S.units.Td       = 'K';
S.units.TLCL     = 'K';

S.units.TAS = 'm/s';
S.units.Vc = 'm/s';
S.units.Vi = 'm/s';

S.units.w   = 'kg/kg';
S.units.q   = 'kg/kg';
S.units.e   = 'hPa';
S.units.es  = 'hPa';
S.units.RH  = '1';

S.units.Cp    = 'J/kg/K';
S.units.Cv    = 'J/kg/K';
S.units.Rm    = 'J/kg/K';
S.units.gamma = '1';

S.units.theta   = 'K';
S.units.thetad  = 'K';
S.units.theta_v = 'K';
S.units.theta_e = 'K';

S.units.Zp       = 'm';
S.units.Zp_ft    = 'ft';
S.units.Zp_trop  = 'm';
S.units.Zhydro    = 'm';
S.units.Zvirt    = 'm';
S.units.Zgph     = 'm';
S.units.Zdens    = 'm';
S.units.Zdens_ft = 'ft';

end


% =========================================================================
%   Static temperature from Measured temperature
% =========================================================================
function Ts = Ts_from_Tm(Tm, recovf, Rd, Rv, varargin)
% TS_FROM_TM: Recover static temperature from measured (probe) temperature
%
%   Ts = Ts_from_Tm(Tm, r, Rd, Rv, 'gamma', gamma,'mach',M)
%
% Inputs:
%   Tm   : Measured temperature (K) (usually probe or total temperature sensor)
%   r    : Recovery factor (dimensionless, typically 0.85–1.0)
%   Rd   : Gas constant for dry air   [J/kg/K]  ~ 287.05 (from phycon.m)
%   Rv   : Gas constant for water vapor [J/kg/K] ~ 461.5 (from phycon.m)
%
% Optional Name-Value Input:
%   'gamma' : Ratio of specific heats (Cp/Cv), moist or dry. If not given,
%             we use gamma = 1.4 (dry air value).
%   'mach : Mach number
%
% Output:
%   Ts   : Static air temperature (K)
%
% Notes:
%   This formula is valid for compressible, subsonic flow with
%   isentropic deceleration at the probe head. The gamma used should
%   ideally match the actual moist-air gamma, but if not available,
%   1.4 is acceptable in practice.
%
%   Ts ≈ Tm / [1 + recovf*(gamma - 1)/2*M^2]
%   Here, user must supply M or gamma externally if more precise modeling is needed.

% -- Parse optional gamma input
p = inputParser;
addParameter(p, 'gamma', 1.4, @(x) isnumeric(x));
addParameter(p, 'mach', 0, @(x) isnumeric(x));
parse(p, varargin{:});
gamma = p.Results.gamma;
M = p.Results.mach;

% NOTE:
% This routine assumes that the probe temperature Tm is effectively total
% temperature (Tt), but corrected by r (recovery factor) for probe behavior.
Ts = Tm ./ (1 + recovf*((gamma - 1)/2).*M.^2);

end
