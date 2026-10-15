function [pcorc,fcoef,machn,qx,tbx,tax,XXf,betaf,qx0,f0]=cone_pcor(dp1,dpb,dpa,dpr,psm,varargin);
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
addParameter(p, "dpn",    [],                  @(x) isnumeric(x));
addParameter(p, "Rate",   [],                  @(x) isnumeric(x));
p.parse(varargin{:});
SOURCE   = p.Results.SOURCE;
mr       = p.Results.mr;
dpn      = p.Results.dpn;
rate     = p.Results.Rate;

% sanity
if mean(dpa)<0
    dpa=-dpa;
end
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
tax = tanAlpha(dpa,dpb,dpr);
tbx = tanBeta(dpb,dpr);

[ta, tb, K, flags] = solveTaTb(dpa, dpb, dpr, dpn);
S    = 1 + ta.^2 + tb.^2;
fq   = K.*S;

% Model predicts alpha variablility only.  
% Predict differential pressures (+/-45 degrees)
%    along actual flow direction when sideslipping
% Difference with dpa is largest when dpb is large and dpa is small.
theta   = atan(hypot(ta, tb));
x       = fq.*sin(2*theta);                
ype
qx0  = impactPcalc(dp1,dpa,dpb,dpr); %uncorrected
fqx0 = fqCalc(dpa,dpb,dpr);  

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
pcorc = fqx0./f0 -qx0;  %  Error in q
ptotal = qx0 + psm;
for jj=1:3 % Iterate three times
    machn = mach(qx0 + pcorc, psm - pcorc, mr);
    %%%XX = [machn machn.^2 dpa]; 
    XX = [machn machn.^2 dptheta]; 
    XXf=[onez XX];
    f0=XXf*betaf0; 
    q    = fq./f0;
    pcorc = q - fq.*(ta.^2 + tb.^2)./S - dp1;
end

% clamp the endpoints
qx = impactPcalc(dp1+pcorc,dpa,dpb,dpr); %uncorrected
fcoef = f0;
if numel(kk)>10
    k1 = kk(1); k2 = kk(end);
    qx   (1:k1-1) = qx(k1);     qx   (k2+1:end) = qx(k2);
    fcoef(1:k1-1) = fcoef(k1);  fcoef(k2+1:end) = fcoef(k2);
    tax  (1:k1-1) = tax(k1);    tax  (k2+1:end) = tax(k2);
    tbx  (1:k1-1) = tbx(k1);    tbx  (k2+1:end) = tbx(k2);
    pcorc(1:k1-1) = pcorc(k1);  pcorc(k2+1:end) = pcorc(k2);
end

return
end

function [ta, tb, K, flags] = solveTaTb(dpa, dpb, dpr, dpn)
% Branch estimates of K = fq/S
Kb = Kroot(2*dpr + dpb, dpb.^2/4);   % uses dpb, dpr
Ka = Kroot(2*dpn + dpa, dpa.^2/4);   % uses dpa, dpn

badB = (dpb == 0);                   % dropout flags
badA = (dpa == 0);
both = badA & badB;

K = (Ka + Kb)/2;                     % normal case: average the branches
K(badB) = Ka(badB);                  % dpb dead -> trust dpa/dpn branch
K(badA) = Kb(badA);                  % dpa dead -> trust dpb/dpr branch

ta = dpa./(2*K);
tb = dpb./(2*K);
tb(badB) = -1 + sqrt(max(2 - 2*dpr(badB)./K(badB), 0));
ta(badA) = -1 + sqrt(max(2 - 2*dpn(badA)./K(badA), 0));

ta(both) = NaN; tb(both) = NaN; K(both) = NaN;

flags = struct('badA', badA, 'badB', badB, 'both', both, ...
    'relRes', (Ka - Kb)./max(Ka, Kb));
end

function K = Kroot(b, c)
% Positive root of K^2 - b*K - c = 0 (c >= 0), stable for either sign of b
s = sqrt(b.^2 + 4*c);
K = (b + s)/2;
neg = b < 0;
K(neg) = 2*c(neg)./(s(neg) - b(neg));
end
