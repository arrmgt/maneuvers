clear;close all
% Restore factory-default MATLAB path for this session
restoredefaultpath
addpath('c:/users/rodi/Github/maneuvers/processing26/utils')
addpath('c:/users/rodi/Github/kingair-Sys26/Sys26/get_vars/utilities')
addpath('c:/users/rodi/Github/kingair-Sys26/Sys26/get_vars/mfiles858')
addpath('c:/users/rodi/Github/kingair-Sys26/ARRfiles')

close all
zzzz=0;
yyyy=0;
;
PROJ = 'test26'
FLTs = ["20260408a_arr.c10.nc" "20260408b_arr.c10.nc"]
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

bb = 0;
iiis0=1; % this is to start the accumulation indices of maneuver sections
kkks0=nan(nn*mm,1); % this  will be where the segment indices are saved
kkks1=nan(nn*mm,1);

C=phycon;

BETA_ship = [0 0 0 0 0];
BETA_boom = [0 0 0 0 0];

PRESSURE = {'ship','boom'};

BETAX = zeros(5,numel(PRESSURE));

for pp = 1:numel(PRESSURE);
    PCOR = 0;
    BETAXX = [];
    % The *1 variables are the concatinated series
    alpha1      = [];
    beta1       = [];
    q_impact1   = [];
    zhydro1     = [];
    zalt1       = [];
    Mnum1       = [];
    DelP1       = [];
    Yhat1       = [];

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
        yy = 1000:6000;
        PSoffset = mean(PSB0(yy)-PSA0(yy));
        PSB = PSB0; 
        PSA = PSA0 + PSoffset;
        
        [MM, MMboom, MMship] = getDerivedVariablesR858( ...
                DPB, DPA, DPR, DPN, DP1, DP2, PSA, PSB);
    
        mr = zeros(size(DPA));
    
        ta = MM.ta_beta;
        tb = MM.tb_beta;
        alpha = atan(ta);
        beta = atan(tb);
        
        q_ship =  MMship.q_beta;
        q_boom =  MMboom.q_beta;
        
        [Mnum_ship,TEMPK_ship, TAS_ship]  = calc_mach( ...
            MMship.q_beta, PSA, TROSE, 0.97);
        [Mnum_boom,TEMPK_boom, TAS_boom]  = calc_mach( ...
            MMboom.q_beta, PSB, TROSE, 0.97);
    
        rho_ship = PSA.*100./(C.Rd .* TEMPK_ship);
        rho_boom = PSB.*100./(C.Rd .* TEMPK_boom);
    
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
                for i = 1:3  % iterate to get pcor
                    PS = PSA ;
                    QC = q_ship ;
                    TS = tstatic(TROSE,0.97,QC,PS,mr);
                    Mnum = calc_mach(QC, PS, TROSE, 0.97);
                    rho_ship = PS*100 ./ (C.Rd .* TS);
                    pZgps  = PSA(1)   + cumtrapz( zalt0, -1 .* rho_ship .* g )./100;  
                    zhydro = zalt0(1) + cumtrapz( PS*100, -1 ./ (rho_ship .* g) );
                    p = polyfit(pZgps, PS-pZgps, 1);
                    pZerr = polyval(p,pZgps);
                    x1 = PS -pZerr - pZgps;
                    x1 = PS - Zgps;
                    x2 = Mnum;
                    x3 = QC;
                    x4 = alpha;
                    x5 = beta;
                    X = [x1,x2,x3,x4,x5];
                    variableNames = ["x1","x2","x3","x4","x5"]; 
                    [lm,yhat,BETA0, R2,Tship] = do_regress(X,variableNames);
                    PCOR = calc_pcor(BETA0,X(:,2:end));
                end
    
            case 'boom'
                PCOR = 0;
                for i = 1:3  % iterate to get pcor 
                    PS = PSB ;
                    QC = q_ship ;
                    TS = tstatic(TROSE,0.97,QC,PS,mr);
                    Mnum = calc_mach(QC, PS, TROSE, 0.97);
                    BETA
                    rho_ship = PS*100 ./ (C.Rd .* TS);
                    pZgps  = PSB(1)   + cumtrapz( zalt0, -1 .* rho_ship .* g )./100;  
                    zhydro = zalt0(1) + cumtrapz( PS*100, -1 ./ (rho_ship .* g) );
                    p = polyfit(pZgps, PS-pZgps, 1);
                    pZerr = polyval(p,pZgps);
                    x1 = PS -pZerr - pZgps;
                    x2 = Mnum;
                    x3 = QC;
                    x4 = alpha;
                    x5 = beta;
                    X = [x1,x2,x3,x4,x5];
                    variableNames = ["x1","x2","x3","x4","x5"]; 
                    [lm,yhat,BETA0, R2,Tship] = do_regress(X,variableNames);
                    PCOR = calc_pcor(BETA0,X(:,2:end));
                end
        end
        
        for ii=1:mm;
            clear kk
            ss = sprintf("kk = M.flt(%i,%i).kk;",jj,ii);
            eval(ss);
            if isempty(kk)
                continue
            end
       
            x1 = PS(kk) -pZerr(kk) - pZgps(kk);
            x2 = Mnum(kk);
            x3 = QC(kk);
            x4 = alpha(kk);
            x5 = beta(kk);
            X = [x1,x2,x3,x4,x5];
            varNames = ["x1","x2","x3","x4","x5"];
            [lm,Yhat,BETA, R2,T] = do_regress(X,varNames);

            Yhat1           = [Yhat1;Yhat];
            alpha1          = [alpha1;alpha(kk)];
            beta1           = [beta1;beta(kk)];
            q_impact1       = [q_impact1;QC(kk)];
            zhydro1         = [zhydro1;zhydro(kk)];
            zalt1           = [zalt1;zalt(kk)];
            Mnum1           = [Mnum1;Mnum(kk)];
            DelP1           = [DelP1;x1];
             
        end; %for mm (legs/flight)
    end; %for jj (flights_
    % Now process the concatenated data
    alpha    = alpha1;
    beta     = beta1;
    q_impact = q_impact1;
    zhydro   = zhydro1;
    zalt     = zalt1;
    Mnum     = Mnum1;
    DelP     = DelP1;
    
    x1s = DelP;
    x2s = Mnum;
    x3s = q_impact;
    x4s = alpha;
    x5s = beta;
    X = [x1s,x2s,x3s,x4s,x5s];
    variableNames = ["PCOR","Mach","q","alpha","beta"];
    [lm,yhat,BETAXX, R2,Tship] = do_regress(X,variableNames);
    BETAX(:,pp) = BETAXX;

    figure(pp)
    plot(x3s,yhat,'.');
    title(sprintf('20260408a-b: %s pcor',PRESSURE{pp}),'fontsize',15)
    xlabel('Q\_impact [m/s]','fontsize',15)
    ylabel('PCOR [hPa]','fontsize',15)
    v=axis;
    axis([.9*v(1) v(2), -5, 1])
    grid
    ss=sprintf("print -djpeg 'figs/%s.jpg'",PRESSURE{pp})
    eval(ss)

end; % for pp (pressures)


return

