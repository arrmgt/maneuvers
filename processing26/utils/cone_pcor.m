function [pcorc,fcoef,machn,qx,tbx,tax,XXf,betaf,qx0,f0]=cone_pcor(dp1,pb,pa,pr,psm,varargin);
%     [pcorc0,fcoef,machn0,q0,tbx0,tax0,XXf,betaf,qx0,f0]=
% Calculate flow angles and static pressure correction
%  [pcorc,fcoef,machn,qx,tbx,tax,XXf,betaf,qx0,f
%  [pcorc,qx,tbx,tax,f0,XXf,betaf,machn] = cone_pcor(dp1,pb,pa,pr,psm,varargin);
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


% From Rodi&Leon(2012)
betaf = [ ...
   1.699864444944109; ...
  -0.156929423443038; ...
   0.066325085038090; ...
   0.001254576494439  ...
   ];

if(nargin<6)
    mr = zeros(size(dp1));
else
    mr = varargin{1};
end

% These are independent of static pressure correctitbx = tanBeta(pb,pr);
tax = tanAlpha(pa,pb,pr);
tbx = tanBeta(pb,pr);
qx0 = impactPcalc(dp1,pa,pb,pr); %uncorrected
% fqx us f*q; fqx/f = q; fqx is independent of pcor
fqx = fqCalc(pa,pb,pr);  

dp1_min = 10; %mb
% Sanity check
machn=mach(qx0,psm,mr);

kk = find ( dp1>dp1_min & qx0>20 & qx0<80 & ...
    ((qx0+psm)./psm-1)>0 & psm>200 & psm<1200 );
if ~isempty(kk)
    dp1 = interp1(kk,dp1(kk),[1:numel(dp1)]','linear',0);
    qx0 = interp1(kk,qx0(kk),[1:numel(dp1)]','linear',0);
    psm = interp1(kk,psm(kk),[1:numel(dp1)]','linear',0);
end

onez = ones(size(psm));
% Set default f
f0=1.68.*ones(size(dp1)); % just a guess
%  We need mach number to get f, so we have to iterate
pErr = fqx./f0 -qx0;  %  Error in q
for jj=1:3 % Iterate three times
% We need machn to get pErr, and pErr to get machn
    machn=mach(qx0+pErr,psm-pErr,mr);
    % Rodi & Leon(2012)
    XX=[machn machn.^2 pa]; 
    XXf=[onez XX];
    f0=XXf*betaf; 
    pErr=fqx./f0-qx0;
end

% clamp the endpoints
pcorc = pErr; 
qx = fqx./f0;
fcoef = f0;
k1 = kk(1); k2 = kk(end);
pErr (1:k1-1) = pErr(k1);   pErr (k2+1:end) = pErr(k2);
qx   (1:k1-1) = qx(k1);     qx   (k2+1:end) = qx(k2);
fcoef(1:k1-1) = fcoef(k1);  fcoef(k2+1:end) = fcoef(k2);
tax  (1:k1-1) = tax(k1);    tax  (k2+1:end) = tax(k2);
tbx  (1:k1-1) = tbx(k1);    tbx  (k2+1:end) = tbx(k2);
pcorc(1:k1-1) = pcorc(k1);  pcorc(k2+1:end) = pcorc(k2);

return


