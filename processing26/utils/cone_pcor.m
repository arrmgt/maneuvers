function [pcorc,fcoef,machn,qx,tbx,tax,XXf,betaf,qx0,f0]=cone_pcor(dp1,pb,pa,pr,psm,varargin);
% REVISED for NEW KING AIR
%       [pcorc,fcoef,machn,qx,tbx,tax,XXf,betaf,qx0,f0]=
%               cone_pcor(dp1,pb,pa,pr,psm,varargin);
% Calculate flow angles and static pressure correction 
%      PSM is uncorrected static pressure
%      MR is mixing ratio for humidity correction [g/g] (optional)
%
%   Outputs [pcorc,qx,tbx,tax,f0,XXf,betaf,machn] 
%       pcorc = Static pressure correction 
%       qx = dynamic pressure
%       tbx = tan of sideslip angle beta
%       tba = tan of attack angle alpha
%       f0 = 858 probe sensitivity factor
%       XXF, betaf , machn are for diagnostic checking (see below)
%     
% REVISED FOR NEW KING AIR 20260812
p = inputParser;
addParameter(p, "SOURCE", "SHIP",              @(s)ischar(s)||isstring(s));
addParameter(p, "mr",     zeros(size(dp1)),    @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
SOURCE = p.Results.SOURCE;
mr     = p.Results.mr;

switch SOURCE  % SHIP or BOOM static pressure?
    case "SHIP"
        betaf0 = [ ...
           1.699864444944109; ...
          -0.156929423443038; ...
           0.066325085038090; ...
           0.001254576494439  ...
           ];

    case "BOOM"
        betaf0 = [ ...
           1.699864444944109; ...
          -0.156929423443038; ...
           0.066325085038090; ...
           0.001254576494439  ...
           ];
end
% These are independent of static pressure correction
tax = tanAlpha(pa,pb,pr);
tbx = tanBeta(pb,pr);
abFact = 1 + tax.^2 + tbx.^2;
qx0 = impactPcalc(dp1,pa,pb,pr); %uncorrected
fqx0 = fqCalc(pa,pb,pr);  

dp1_min = 10; %mb
% Sanity check
machn=mach(qx0,psm,mr);

kk = find ( dp1>dp1_min & qx0>20 & qx0<80 & ...
    ((qx0+psm)./psm-1)>0 & psm>200 & psm<1200  & machn<.6);
if ~isempty(kk)
    dp1 = interp1(kk,dp1(kk),[1:numel(dp1)]','linear',0);
    qx0 = interp1(kk,qx0(kk),[1:numel(dp1)]','linear',0);
    psm = interp1(kk,psm(kk),[1:numel(dp1)]','linear',0);
end

onez = ones(size(psm));
% Set default f
f0=1.7.*ones(size(dp1)); % just a guess
%  We need mach number to get f, so we have to iterate
pErr = fqx0./f0 -qx0;  %  Error in q
ptotal = qx0 + psm;
for jj=1:3 % Iterate three times
    machn = mach(qx0 + pErr, psm - pErr, mr);
    XX = [machn machn.^2 pa]; 
    XXf=[onez XX];
    f0=XXf*betaf0; 
    pErr = fqx0./f0 - qx0;
end

% clamp the endpoints
pcorc = pErr; 
qx = qx0 + pErr;
fcoef = f0;
if numel(kk)>10
    k1 = kk(1); k2 = kk(end);
    pErr (1:k1-1) = pErr(k1);   pErr (k2+1:end) = pErr(k2);
    qx   (1:k1-1) = qx(k1);     qx   (k2+1:end) = qx(k2);
    fcoef(1:k1-1) = fcoef(k1);  fcoef(k2+1:end) = fcoef(k2);
    tax  (1:k1-1) = tax(k1);    tax  (k2+1:end) = tax(k2);
    tbx  (1:k1-1) = tbx(k1);    tbx  (k2+1:end) = tbx(k2);
    pcorc(1:k1-1) = pcorc(k1);  pcorc(k2+1:end) = pcorc(k2);
end

return


