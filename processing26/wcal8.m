function [f,fus,fvs,fws,va,vg,mmm]=wcal8(x,Data);
%
%Project $Name:  $ ($Revision: 1.1.1.1 $)
%
%  Inputs Data  |   Data
%  Inputs x     |   Parameters
%       x(1)    Roll offset angle (radians)
%       x(2)    Pitch offset angle (radians)
%       x(3)    Heading offset angle (radians)
%       x(4)    Attack angle scaling factor
%       x(5)    Attack angle offset
%       x(6)    Sideslip angle scaling factor
%       x(7)    Sideslip angle offset
%       x(8)    Pressure offset (hPa)
%
%   Alpha = x(4)*attack_indicated + x(5)
%   Beta  = x(6)*beta_indicated   + x(7)

try
    att     =   Data.att(:,:);      %[roll,pitch,thead]
    pmb     =   Data.air(:,1);      %[pmb,ttotK,trf,mr];
    ttotK   =   Data.air(:,2);
    trf     =   Data.air(:,3); 
    mr      =   Data.air(:,4);

    attr    =   Data.attr(:,:);     %[rollr,pitchr,yawr];
    tas     =   Data.flow(:,1);     %[tas,alpha,beta];
    alpha   =   Data.flow(:,2);
    beta    =   Data.flow(:,3);
    ARM     =   Data.arm(:,:);
    VEW     =   Data.earthv(:,1);
    VNS     =   Data.earthv(:,2);
    VZ      =   Data.earthv(:,3);
    kkks0   =   Data.kkks0(:);
    kkks1   =   Data.kkks1(:);
    diffp   =   Data.boom(:,1:4);   %[dp1,dpb,dpa,dpr,q_impact];
    q_impact=   Data.boom(:,5);
    tzero   =   273.15;
    
    poffset =   x(8)
    dp1     =   Data.boom(:,1) + poffset;
    dpb     =   Data.boom(:,2);
    dpa     =   Data.boom(:,3);
    dpr     =   Data.boom(:,4);
    ps      =   pmb - poffset;
    q_impact = solve858(dp1,dpa,dpb,'dpr',dpr);
    
    rolloff  =   x(1);
    pitchoff =   x(2);
    headoff  =   x(3);
    ROLL     =   att(:,1) + rolloff;
    PITCH    =   att(:,2) + pitchoff;
    HEAD     =   att(:,3) + headoff;
    att      =   [ROLL,PITCH,HEAD]';
    
    OMEGA   =   attr;
    afactor =   [x(4),x(5)];
    bfactor =   [x(6),x(7)];

    ts          = tstatic(ttotK,0.97,q_impact,ps,mr);
    tas         = tasf(q_impact,ps,ts,mr);
    
    % [vg,va]=get_vg(bet0,alp0,tas,att0,omega,arm,bfactor,afactor,rolloff,pitoff,hedoff);
    [vg,va,mmm] = get_vg(beta,alpha,tas,att,OMEGA,ARM,bfactor,afactor,x(1),x(2),x(3));
    fu      =   vg(:,1) - VEW;
    fv      =   vg(:,2) - VNS;

    % kkks0 are the start indices of each flight segment
    % kkks1 are the end indices
    % Each u and v wind segment is detrended and then concatenated 
    %   resulting in a continous array of detrended blocks
    % Vertical wind is not blocked

    fus = []; fvs = [];
    for zz=1:length(kkks0);
        zzz=kkks0(zz):kkks1(zz);
        fus=[fus;detrend(fu(zzz),'constant')];
        fvs=[fvs;detrend(fv(zzz),'constant')];
    end
    fws = vg(:,3) - VZ;
    f(:,1)=fus;
    f(:,2)=fvs;
    f(:,3)=fws;
catch ME
    catchME(ME)
end
sprintf('std = %5.2f mean = %5.2f',std(fus),mean(fus))

