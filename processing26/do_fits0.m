function [Params,fu,fv,fw,mag,dir,resid,jacobian,CI,beta_samp]= ...
    do_fits0(kkks0,kkks1,dp1,dpb,dpa,dpr,q_impact,temp,tas,pmb,mr ...
    ,alpha,beta,roll,pitch,thead,rollr,pitchr,yawr,vew,vns,vz,ARM);

pmb         = pmb;
mr          = zeros(size(pmb));
qcx         = q_impact;
pcorc       = zeros(size(dpb));
recovf      = 0.97;
tzero       = 273.15;
ttotK       = ttotal(temp,recovf,qcx,pmb,mr);

ROL         = roll;
PIT         = pitch;
HED         = thead;
ROLR        = rollr;
PITR        = pitchr;
YAWR        = yawr;
VEW         = vew;
VNS         = vns;
VZ          = vz;

ARM=ARM(ones(length(tas),1),:);

X.att       = [ROL, PIT, HED];
X.attr      = [ROLR, PITR, YAWR];
X.air       = [pmb, ttotK, temp, mr];
X.arm       = ARM;
X.boom      = [dp1, dpb, dpa, dpr, q_impact];
X.flow      = [tas, alpha, beta];
X.arm       = ARM;
X.earthv    = [VEW, VNS, VZ];
X.kkks0     = kkks0;
X.kkks1     = kkks1;

method='lsqnonlin';
switch method
    case 'lsqnonlin'
        options=optimset('lsqnonlin');
    case 'fsolve'
        options=optimset('fsolve');
end
options=optimset(options,'Display','iter');
options=optimset(options,'MaxFunEvals',10000);
options=optimset(options,'TolX',5.e-3,'TolFun',5.e-3);

%%fun=@(x)wcal8(x,X); % pass data X structure to function
fun=@wcal8; % pass data X structure to function
% [x,fu,fv,fw,mag,dir,resid,jacobian,CI,beta_samp]=do_fit1(method,fun,X)
try
   [Params,fu,fv,fw,mag,dir,resid,jacobian,CI,beta_samp]=do_fit1(method,fun,X);
catch ME
    catchME(ME)
end
%