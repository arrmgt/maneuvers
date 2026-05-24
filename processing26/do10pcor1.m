clear;close all
% Restore factory-default MATLAB path for this session
restoredefaultpath
addpath('c:/users/rodi/Github/maneuvers/processing26/utils')
addpath('c:/users/rodi/Github/kingair-Sys26-work/Sys26')
addpath('c:/users/rodi/Github/kingair-Sys26-work/Sys26/get_vars/utilities')
addpath('c:/users/rodi/Github/kingair-Sys26-work/Sys26/get_vars/vapor_press_functions/')
addpath('c:/users/rodi/Github/kingair-Sys26-work/Sys26/get_vars/mfiles858')
addpath('c:/users/rodi/Github/kingair-Sys26-work/ARRfiles')

close all
zzzz=0;
;
PROJ = 'test26'
FLTs = ["20260408a_arrWork1.c10.nc" "20260408b_arrWork1.c10.nc"]
dataDir = fullfile('P:/MATLAB-DATA2/kingair_data/',PROJ,'work');
file = fullfile(dataDir,FLTs{1});
jrate = get_irate(file,'TASX');
orate = 10;
r = orate/jrate;
% these are for 10 Hz
M.flt(1,1).kk = round([38638*r:50369*r]); %Rodi
M.flt(1,2).kk = round([75808*r:86897*r]); %Rodi
M.flt(1,3).kk = round([56137*r:57087*r]); % betas S/L
M.flt(1,4).kk = round([91029*r:91794*r]); % betas S/L

M.flt(2,1).kk = round([22720*r:34908*r]); %Rodi
M.flt(2,2).kk = round([37567*r:38325*r]); % betas S/L
M.flt(2,3).kk = round([42614*r:43668*r]); % betas S/L
M.flt(2,4).kk = [];

[nn,mm] = size(M.flt);

% The *1 variables are the concatinated series
alpha1      = [];
beta1       = [];
q_impact1   = [];
zhydro1     = [];
zalt1       = [];
Mnum1       = [];
DelP1       = [];
Yhat1       = [];
DP1_1       = [];
DP2_1       = [];
DPA1        = [];
DPB1        = [];
DPR1        = [];
DPN1        = [];
PSA1        = [];
PSB1        = [];
PTB1        = [];

bb = 0;
iiis0=1; % this is to start the accumulation indices of maneuver sections
kkks0=nan(nn*mm,1); % this  will be where the segment indices are saved
kkks1=nan(nn*mm,1);

C=phycon;

PRESSURE = {'ship','boom'};
for pp = 1:numel(PRESSURE);
    PRESSURE{pp}
for jj = 1:numel(FLTs)
    arcfile = fullfile(dataDir,FLTs(jj));
    % For string array
    lastUnd = cellfun(@(s) find(s=='_',1,'last'), cellstr(arcfile));
    rawfile = [ extractBefore(arcfile, lastUnd) +  "_raw.nc"]; 

    ARM = ncreadatt(arcfile,'/','AWinds.MomentArm');
    
    % raw measurements
    irate = get_irate(rawfile,'AALT')
    blurf = get_data(rawfile,'AALT',[],irate,orate); zalt0=blurf(:);
    blurf = get_data(rawfile,'PSA',[],1000,orate); PSA0=blurf(:);
    blurf = get_data(rawfile,'PSB',[],1000,orate); PSB0=blurf(:);
    blurf = get_data(rawfile,'TROSE',[],1000,orate); TROSE=blurf(:) +C.Tzero;
    blurf = get_data(rawfile,'PTB',[],1000,orate); PTB=blurf(:);
    blurf = get_data(rawfile,'DPA',[],1000,orate); DPA=blurf(:);
    blurf = get_data(rawfile,'DPB',[],1000,orate); DPB=blurf(:);
    blurf = get_data(rawfile,'DPR',[],1000,orate); DPR=blurf(:);
    blurf = get_data(rawfile,'DP1',[],1000,orate); DP1=blurf(:);
    blurf = get_data(rawfile,'DP2',[],1000,orate); DP2=blurf(:);
    blurf = get_data(rawfile,'DPN',[],1000,orate); DPN=blurf(:);
    blurf = get_data(rawfile,'AROLL',[],irate,orate); ROLL=blurf(:);
    blurf = get_data(rawfile,'APITCH',[],irate,orate); PITCH=blurf(:);
    blurf = get_data(rawfile,'AHEAD',[],irate,orate,true); HEAD=blurf(:);
    blurf = get_data(rawfile,'AWANDER',[],irate,orate,true); WANDER=blurf(:);
    blurf = get_data(arcfile,'AVthead',[],jrate,orate,true); AVthead=blurf(:);

    % Assume PSB is calibrated
    zz = 750/25*orate:750/25*orate; % 30 seconds data
    PSoffset = mean( PSA0(zz)-PSB0(zz) ); % before takeoff
    PSA = PSA0 - PSoffset;
    PSB = PSB0;

    [MM, MMboom, MMship] = getDerivedVariablesR858( ...
            DPB, DPA, DPR, DPN, DP1, DP2, PSA, PSB);

    mr = zeros(size(DPA));
    [q_ship, ~, ta, tb] = solve858(DP1, DPA, DPB, 'dpr', DPR, 'DPN', DPN);
    [q_boom, ~, ta, tb] = solve858(DP2, DPA, DPB, 'dpr', DPR, 'DPN', DPN);
    alpha = atan(ta);
    beta  = atan(tb);
    Ptot_ship  = q_ship + PSA;
    Ptot_boom  = q_boom + PSB;
    TDPK = -40.*ones(size(PSA)) + C.Tzero;
    ADatShip   = airdata(PSA, Ptot_ship, TROSE, 0.97, TDPK);
    ADatBoom   = airdata(PSB, Ptot_boom, TROSE, 0.97, TDPK);
    Mnum_ship  = ADatShip.M;
    Mnum_boom  = ADatBoom.M;
    TEMPK_ship = ADatShip.Ts;
    TEMPK_boom = ADatBoom.Ts;
    rho_ship   = ADatShip.rho;
    rho_boom   = ADatBoom.rho;

    % correct altitude for boom location
    nnn = numel(blurf);
    ARM = [6.77,0.50,-0.44];
    arm = repmat(ARM,nnn,1);
    attZ = [ROLL,PITCH,AVthead]'.*pi./180;
    Zarm = va2vg(attZ,arm)';
    zalt = zalt0 - Zarm(:,3);
    
    g = C.g0;
    switch PRESSURE{pp}
        case 'ship'
            pZgps  = PSA(1)   + cumtrapz( zalt, -1 .* rho_ship .* g )./100;  
            zhydro = zalt0(1) + cumtrapz( PSA*100, -1 ./ (rho_ship .* g) );
            p = polyfit(pZgps, PSA-pZgps, 1);
            pZerr = polyval(p,pZgps);
            Mnum = Mnum_ship;
            q_impact = q_ship;
            PS = PSA;

        case 'boom'
            pZgps  = PSB(1)   + cumtrapz(zalt0, -1 .* rho_boom .* g)./100;
            zhydro = zalt0(1) + cumtrapz(PSB*100, -1 ./ (rho_boom .* g));
            p = polyfit(pZgps, PSB-pZgps, 1);
            pZerr = polyval(p,pZgps);
            Mnum = Mnum_boom;
            q_impact = q_boom;
            PS =  PSB;
    end
   
    for ii=1:mm;
        clear kk
        ss = sprintf("kk = M.flt(%i,%i).kk;",jj,ii);
        eval(ss);
        if isempty(kk)
            continue
        end

        x1 = - [PS(kk) - pZerr(kk) - pZgps(kk)];
        x2 = Mnum(kk);
        x3 = q_impact(kk);
        x4 = alpha(kk);
        x5 = beta(kk);
        x6 = tan(alpha(kk)).^2 + tan(beta(kk)).^2;
        X = [x1,x2,x3,x4,x5,x6];
        [lm,Yhat,BETA0, R2,T] = do_regress(X);

        Yhat1           = [Yhat1;Yhat];
        alpha1          = [alpha1;alpha(kk)];
        beta1           = [beta1;beta(kk)];
        q_impact1       = [q_impact1;q_impact(kk)];
        zhydro1         = [zhydro1;zhydro(kk)];
        zalt1           = [zalt1;zalt(kk)];
        Mnum1           = [Mnum1;Mnum(kk)];
        DelP1           = [DelP1;x1];
        DP1_1           = [DP1_1;DP1(kk)];
        DP2_1           = [DP2_1;DP2(kk)];
        DPA1            = [DPA1;DPA(kk)];
        DPB1            = [DPB1;DPB(kk)];
        DPR1            = [DPR1;DPR(kk)];
        DPN1            = [DPN1;DPN(kk)];
        PSA1            = [PSA1;PSA(kk)];
        PSB1            = [PSB1;PSB(kk)];
        PTB1            = [PTB1;PTB(kk)];
        
    end; %for ii
end; %for jj
% Now process the concatenated data
alpha    = alpha1;
beta     = beta1;
q_impact = q_impact1;
zhydro   = zhydro1;
zalt     = zalt1;
Mnum     = Mnum1;
DelP     = DelP1;
DP1      = DP1_1;
DP2      = DP2_1;
DPA      = DPA1;
DPB      = DPB1;
DPR      = DPR1;
DPN      = DPN1;
PSA      = PSA1;
PSB      = PSB1;
PTB      = PTB1;

if contains(PRESSURE{pp},'ship')
    [q0, f0, ta0, tb0] = solve858(DP1, DPA, DPB, ...
    'dpr', DPR, 'DPN', DPN);
elseif contains(PRESSURE{pp},'boom')
    [q0, f0, ta0, tb0] = solve858(DP2, DPA, DPB, ...
    'dpr', DPR, 'DPN', DPN);
else
    error('Arrrgggg')
end

x1s = DelP;
x2s = Mnum;
x3s = q0;
x4s = alpha;
x5s = beta;
x6s = ta0.^2 + tb0.^2;
X = [x1s,x2s,x3s,x4s,x5s,x6s];
sprintf('Run %i:  pressure is %s',pp,upper(PRESSURE{pp}))
[lm,yhat,BETA0, R2,Tship] = do_regress(X);
zzzz = zzzz+1;
BETA(:,zzzz) = BETA0;

[q0, f0, ta0, tb0] = solve858(DP1, DPA, DPB, ...
    'dpr', DPR, 'DPN', DPN);
[q, f, ta, tb] = solve858(DP1 + yhat, DPA, DPB, ...
    'dpr', DPR, 'DPN', DPN, 'ps_cor',yhat);


figure(zzzz)
plot(x3s,yhat,'.');
title(sprintf('20260408a-b: %s pcor',PRESSURE{pp}),'fontsize',15)
xlabel('Q\_impact [m/s]','fontsize',15)
ylabel('PCOR [hPa]','fontsize',15)
v=axis;
axis([.9*v(1) v(2), -5, 1])
grid
ss=sprintf("print -djpeg 'figs/%s.jpg'",PRESSURE{pp})
eval(ss)

end; % for pp

return

