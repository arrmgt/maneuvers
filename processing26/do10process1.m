clear;close all
restoredefaultpath
addpath('utils');

nnn=0;
;
FLTs = ["20260408a_arrWork.c10.nc" "20260408b_arrWork.c10.nc"];
DATE = extractBefore(FLTs{1},'_');
s1=extractAfter(FLTs{1},'.c');
irate=str2num(extractBefore(s1,'.'));
orate = 10;
r = irate/orate;
% Segment indices are for irate = 10 Hz data. Adjust for other rates
blurf = true; % include straight and level S/L leg?
if blurf
    M.flt(1,1).kk = unique([ceil(38638*r):ceil(50369*r)]); %Rodi
    M.flt(1,2).kk = unique([ceil(75808*r):ceil(86897*r)]); %Rodi
    M.flt(1,3).kk = unique([ceil(56137*r):ceil(57087*r)]); % betas S/L
    M.flt(1,4).kk = unique([ceil(91029*r):ceil(91794*r)]); % betas S/L
    %
    M.flt(2,1).kk = unique([ceil(22720*r):ceil(34908*r)]); %Rodi
    M.flt(2,2).kk = unique([ceil(37567*r):ceil(38325*r)]); % betas S/L
    M.flt(2,3).kk = unique([ceil(42614*r):ceil(43668*r)]); % betas S/L
    M.flt(2,4).kk = []
    % these are for 10 Hz
    M.flt(1,1).kk = unique([round(38638*r):round(50369*r)]); %Rodi
    M.flt(1,2).kk = unique([round(75808*r):round(86897*r)]); %Rodi
    
    M.flt(2,1).kk = unique([round(22720*r):round(34908*r)]); %Rodi
    M.flt(2,2).kk = [];
end
[nn,mm] = size(M.flt);

PROJ = 'test26'
% The *1 variables will have the concatinated series
tas1=[];
alpha1=[];
beta1=[];
temp1=[];
pmb1=[];
pcorc1=[];
dp1_1=[];
dpa1=[];
dpb1=[];
dpr1=[];
dp11=[];
roll1=[];
pitch1=[];
thead1=[];
rollr1=[];
pitchr1=[];
yawr1=[];
vew1=[];
vns1=[];
vz1=[];
norma1=[];
lata1=[];
longa1=[];
ARM1=[];
q_impact1 = [];

bb = 0;
iiis0=1; % this is to start the accumulation indices of maneuver sections
kkks0=nan(nn*mm,1); % this  will be where the segment indices are saved
kkks1=nan(nn*mm,1);

C=phycon();
delete figs/*.jpg
PRES = ["SHIP" "BOOM"];
for pppp=1:numel(PRES)
    PRESSURE = PRES(pppp);

for jj = 1:numel(FLTs)
    dataDir = fullfile('P:/MATLAB-DATA2/kingair_data/',PROJ,'work');
    arcFile = fullfile(dataDir,FLTs(jj));
    lastUnd = cellfun(@(s) find(s=='_',1,'last'), cellstr(arcFile));
    rawFile = [ extractBefore(arcFile, lastUnd) +  "_raw.nc"]; 
    ARM     = ncreadatt(arcFile,'/','AWinds.MomentArm');   
    
    % Get the data
    
    % raw measurements
    jrate = get_irate(rawFile,'AALT')
    blurf = get_data(rawFile,'AALT' ,[],jrate,orate); Zgps=blurf(:);
    blurf = get_data(rawFile,'PSA'  ,[],1000,orate); PSA0=blurf(:);
    blurf = get_data(rawFile,'PSB'  ,[],1000,orate); PSB0=blurf(:);
    blurf = get_data(rawFile,'TROSE',[],1000,orate); TROSE=blurf(:) +C.Tzero;
    blurf = get_data(rawFile,'PTB'  ,[],1000,orate); PTB=blurf(:);
    blurf = get_data(rawFile,'DPA'  ,[],1000,orate); DPA=blurf(:);
    blurf = get_data(rawFile,'DPB'  ,[],1000,orate); DPB=blurf(:);
    blurf = get_data(rawFile,'DPR'  ,[],1000,orate); DPR=blurf(:);
    blurf = get_data(rawFile,'DP1'  ,[],1000,orate); DP1=blurf(:);
    blurf = get_data(rawFile,'DP2'  ,[],1000,orate); DP2=blurf(:);
    blurf = get_data(rawFile,'DPN'  ,[],1000,orate); DPN=blurf(:);
    blurf = get_data(rawFile,'PTB'  ,[],1000,orate); PTB=blurf(:);

    %%%%%%%%%% TEST  Assume DPB is calibrated, adjust PSA
    t1sec = 1000/25;
    zz = [t1sec*orate]:[(t1sec+30)*orate]; % 30 sec at beginning
    PSoffset = mean( PSA0(zz)-PSB0(zz) ); % before takeoff
    PSA = PSA0 - PSoffset;
    PSB = PSB0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Get static pressure corrections
    %%%%%%%%%%    and remove outliers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get fcoef from 2005 trailing cone calibration flight
    [~,ship_fcoef] = cone_pcor(DP1,DPB,DPA,DPR,PSA);
    % Using fcoef, calculate pcor from 858 equations.
    ship_pcor = pcor858(ship_fcoef,DP1,DPA,DPB,DPR,PSA); %pcor858(fcoef,DP1,DPA,DPB,DPR,PSA)
    %  Assume PSB-boom_pcor = PSA-ship_pcor, and estimate boom_pcor.
    boom_pcor       = PSB - PSA + ship_pcor;
    kk=find(DP1>10 & DP1<85 & DP2>10 & DP2<85);
    boom_pcor       = interp1(kk,boom_pcor(kk) ,[1:numel(DPA)]','linear');
    ship_pcor       = interp1(kk,ship_pcor(kk) ,[1:numel(DPA)]','linear');
    ship_fcoef      = interp1(kk,ship_fcoef(kk),[1:numel(DPA)]','linear');
    ship_pcor       = setOutPoints(kk,ship_pcor);
    boom_pcor       = setOutPoints(kk,boom_pcor);
    ship_fcoef      = setOutPoints(kk,ship_fcoef);
    
    %%%%%%%%%% Maneuver results
    PSfactor        = ncreadatt(arcFile,'/','AWinds.PStaticOffset'); 
    PSfactor = 0;
    boom_pcor = boom_pcor + PSfactor;
    ship_pcor = ship_pcor + PSfactor;

    TDPK = -40.*ones(size(DPA)) + C.Tzero;
    
    switch PRESSURE
        case 'SHIP'
            dp1         = DP1 + ship_pcor;
            pmb         = PSA - ship_pcor;
            q_impact    = solve858(dp1, DPA, DPB, 'dpr', DPR);
            pTotal      = q_impact + pmb;
            r           = 0.97; % recovery coefficient for Temp
            ADat        = airdata(pmb, pTotal, TROSE, r, TDPK ...
                            , "Z_gps", Zgps); 
        case 'BOOM'
            dp1         = DP2 + boom_pcor;
            pmb         = PSB - boom_pcor;
            q_impact    = solve858(DP2, DPA, DPB, 'dpr', DPR);
            pTotal      = q_impact + pmb;
            r           = 0.97;
            ADat        = airdata(pmb, pTotal, TROSE, r, TDPK ...
                            , "Z_gps", Zgps); 
    end
    tas     = ADat.TAS;
    temp    = ADat.Ts;
    alpha   = atan(tanAlpha(DPA,DPB,DPR));
    beta    = atan(tanBeta(DPB,DPR));
    dpa     = DPA;
    dpb     = DPB;
    dpr     = DPR;
    ias     = ADat.Vi;
    mr      = zeros(size(dpr));
    pcorc   = zeros(size(dpr));
   
    % mixing ratio wasn't always archived -- but set it to zero for
    % now.
    if(numel(mr)<=1)
        mr=zeros(size(tas));
    end
    
    npitch1=[];npitch2=[];npitch3=[];
    try
        [~,npitch1] = get_data(arcFile,'AVpitch',[],jrate,orate);
    catch
        [~,npitch2] = get_data(arcFile,'avpitch',[],jrate,orate);
    end
    
    % Use the post-processed applanix values if available     
    if isempty(npitch1)
        'av variables used'
        jrate   = get_irate(arcFile,'avpitch')
        pitch   = get_data(arcFile,'avpitch'    ,[],jrate,orate);
        roll    = get_data(arcFile,'avroll'     ,[],jrate,orate);
        thead   = get_data(arcFile,'avthead'    ,[],jrate,orate);
        pitchr  = get_data(arcFile,'avpitchr'   ,[],jrate,orate);
        rollr   = get_data(arcFile,'avrollr'    ,[],jrate,orate);
        yawr    = get_data(arcFile,'avyawr'     ,[],jrate,orate);
        ewvel   = get_data(arcFile,'avewvel'    ,[],jrate,orate);
        nsvel   = get_data(arcFile,'avnsvel'    ,[],jrate,orate);
        zvel    = get_data(arcFile,'avzvel'     ,[],jrate,orate);
        norma   = get_data(arcFile,'avnorma'    ,[],jrate,orate);
        lata    = get_data(arcFile,'avlata'     ,[],jrate,orate);
        longa   = get_data(arcFile,'avlonga'    ,[],jrate,orate);
    else
        'AV variables used'
        jrate   = get_irate(arcFile,'AVpitch')
        pitch   = get_data(arcFile,'AVpitch'    ,[],jrate,orate);
        roll    = get_data(arcFile,'AVroll'     ,[],jrate,orate);
        thead   = get_data(arcFile,'AVthead'    ,[],jrate,orate);
        pitchr  = get_data(arcFile,'AVpitchr'   ,[],jrate,orate);
        rollr   = get_data(arcFile,'AVrollr'    ,[],jrate,orate);
        yawr    = get_data(arcFile,'AVyawr'     ,[],jrate,orate);
        ewvel   = get_data(arcFile,'AVewvel'    ,[],jrate,orate);
        nsvel   = get_data(arcFile,'AVnsvel'    ,[],jrate,orate);
        zvel    = get_data(arcFile,'AVzvel'     ,[],jrate,orate);
        norma   = get_data(arcFile,'AVnorma'    ,[],jrate,orate);
        lata    = get_data(arcFile,'AVlata'     ,[],jrate,orate);
        longa   = get_data(arcFile,'AVlonga'    ,[],jrate,orate);
    end
    pitch   = deg2rad(pitch);
    roll    = deg2rad(roll);
    thead   = deg2rad(thead);

    for ii=1:mm;
        clear kk
        ss = sprintf("kk = M.flt(%i,%i).kk;",jj,ii);
        eval(ss);
        if isempty(kk)
            continue
        end
    
        % save indices for this flight and run the regression
        %  to find the Params individually for the flights.
        KKKS0=1;
        KKKS1=numel(kk);
        [Params,fu,fv,fw,mag,dir,resid,jacobian,CI]= ...
        do_fits0(KKKS0,KKKS1,dp1(kk),dpb(kk),dpa(kk),dpr(kk),q_impact(kk), ...
            temp(kk),tas(kk),pmb(kk),mr(kk),alpha(kk),beta(kk), ...
            roll(kk),pitch(kk),thead(kk),rollr(kk),pitchr(kk),yawr(kk), ...
            ewvel(kk),nsvel(kk),zvel(kk),ARM);
    
        % Save the indices 
        bb = bb  + 1;
        kkks0(bb)=iiis0;
        kkks1(bb)=iiis0+length(kk)-1;
        fname{bb}=FLTs(jj);
        iiis0=iiis0+length(kk);
        
        %  Save things in the structureManeuvers
        nnn=nnn+1; % row number
        structManeuvers(nnn,1).fname=FLTs(jj);
        structManeuvers(nnn,1).pressure=round(mean(pmb(kk)));
        structManeuvers(nnn,1).Params=Params;
        structManeuvers(nnn,1).CIlow=CI(:,1)';
        structManeuvers(nnn,1).CIhi=CI(:,2)';
    
        dp1_1=[dp1_1;dp1(kk)];
        alpha1=[alpha1;alpha(kk)];
        beta1=[beta1;beta(kk)];
        temp1=[temp1;temp(kk)];
        tas1=[tas1;tas(kk)];
        dpa1=[dpa1;dpa(kk)];
        dpb1=[dpb1;dpb(kk)];
        dpr1=[dpr1;dpr(kk)];
        dp11=[dp11;dp1(kk)];
        pmb1=[pmb1;pmb(kk)];
        pcorc1=[pcorc1;pcorc(kk)];
        roll1=[roll1;roll(kk)];
        pitch1=[pitch1;pitch(kk)];
        thead1=[thead1;thead(kk)];
        rollr1=[rollr1;rollr(kk)];
        pitchr1=[pitchr1;pitchr(kk)];
        yawr1=[yawr1;yawr(kk)];
        vew1=[vew1;ewvel(kk)];
        vns1=[vns1;nsvel(kk)];
        vz1=[vz1;zvel(kk)];
        norma1=[norma1;norma(kk)];
        lata1=[lata1;lata(kk)];
        longa1=[longa1;longa(kk)];
        q_impact1=[q_impact1;q_impact(kk)];
    end;    %  end ii (segments)
end;        %  end jj (flights)
x1={structManeuvers.Params};
y1={structManeuvers.CIlow};
z1={structManeuvers.CIhi};
nnn=nnn+2;
structManeuvers(nnn,1).fname='means';
structManeuvers(nnn,1).Params=mean(cell2mat(x1(:)));
structManeuvers(nnn,1).CIlow=mean(cell2mat(y1(:)));
structManeuvers(nnn,1).CIhi=mean(cell2mat(z1(:))); 

nnn=nnn+1;
structManeuvers(nnn,1).fname='stdevs';
structManeuvers(nnn,1).Params=std(cell2mat(x1(:)));
structManeuvers(nnn,1).CIlow=std(cell2mat(y1(:)));
structManeuvers(nnn,1).CIhi=std(cell2mat(z1(:)));

% Now process the concatenated data
kkks0=kkks0(~isnan(kkks0));
kkks1=kkks1(~isnan(kkks1));
dp1=dp1_1;
tas=tas1;
alpha=alpha1;
beta=beta1;
temp=temp1;
pmb=pmb1;
pcorc=pcorc1;
dpa=dpa1;
dpb=dpb1;
dpr=dpr1;
dp1=dp11;
vew=vew1;
vns=vns1;
vz=vz1;
roll=roll1;
pitch=pitch1;
thead=thead1;
rollr=rollr1;
pitchr=pitchr1;
yawr=yawr1;
mr=zeros(size(tas));
q_impact = q_impact1;
pcorc = zeros(size(dpa));

[Params,fu,fv,fw,mag,dir,resid,jacobian,CI,beta_samp] = ...
    do_fits0(kkks0,kkks1,dp1,dpb,dpa,dpr,q_impact,temp,tas,pmb,mr ...
    ,alpha,beta,roll,pitch,thead,rollr,pitchr,yawr,vew,vns,vz,ARM);
f=[fu,fv,fw];

nnn=nnn+2; % row number
structManeuvers(nnn,1).fname='Concat';
structManeuvers(nnn,1).Params=Params;
structManeuvers(nnn,1).CIlow=CI(:,1)';
structManeuvers(nnn,1).CIhi=CI(:,2)';

%  Save the results in a spreadsheet
xlsout=["results_all_" + PRESSURE + ".xlsx"];
delete(xlsout)
clear Tnew
Tnew= [struct2table(structManeuvers)];
writetable(Tnew,xlsout,'Sheet',1,'Range','A1')


% Plot results
nn = 10*pppp;
nn=nn+1;
figure(nn)
h=plot((1:length(fu))/orate,detrend([fw,fu,fv],'constant'))
flight='Concat';
jju=1:length(tas);
set(h(1),'LineWidth',1.5);
set(h(2),'LineWidth',1.5);
set(h(3),'LineWidth',1.5);
ss=sprintf('%s flights concatenated (%s PS)',DATE,PRESSURE)
title(ss)
xlabel('Time [secs]')
ylabel('Wind component [m/s]')
legend('up','east','north','location','northeast')
grid
ss=sprintf('figs/Detrended-comps-%s-PS.jpg',PRESSURE);
saveas(gcf,ss,'jpg')

nn=nn+1;
figure(nn)
h=plot((1:length(jju))/orate,[alpha(jju),beta(jju)].*180./pi)
set(h,'LineWidth',1)
xlabel('Relative time [sec]')
ylabel('Flow angle [deg]')
ss=sprintf('%s flights concatenated (%s PS)',DATE,PRESSURE)
title(ss)
legend('Attack angle','Sideslip angle','location','northwest')
grid
ss=sprintf('figs/Flow angles-%s-PS.jpg',PRESSURE);
saveas(gcf,ss,'jpg')

nn=nn+1;
figure(nn)
edges=linspace(-5,5,100);
centers1=edges2centers(edges);
edges1=edges(ones(3,1),:)';
centers1=centers1(ones(3,1),:)';
[n1,bin]=histcounts(detrend(f(:,1)),edges);
[n2,bin]=histcounts(detrend(f(:,2)),edges);
[n3,bin]=histcounts(detrend(f(:,3)),edges);
n=[n1;n2;n3]';
h=plot(centers1,n)
set(h,'LineWidth',1.5)
ss=sprintf('%s flights concatenated (%s PS)',DATE,PRESSURE)
title(ss)
xlabel('Residual [m/s]')
ylabel('Frequency of occurence');
v=axis;
y=0.9*v(4);
x=v(1)+0.1*abs(v(1));
dy=(v(4)-v(3))/20;
text(x,y,sprintf('std u = %5.3f',std(f(:,1))))
y=y-dy;
text(x,y,sprintf('std v = %5.3f',std(f(:,2))))
y=y-dy;
text(x,y,sprintf('std w = %5.3f',std(f(:,3))))
y=y-dy;
legend('u','v','w')
grid
ss=sprintf('figs/Resids_hist-%s-PS.jpg',PRESSURE);
saveas(gcf,ss,'jpg')


nn=nn+1;
figure(nn)
boxplot([fu,fv,fw],'notch','on','whisker',1,'Symbol','.' ...
    ,'Labels',{'East-component','North-component','Up-component'})
ylabel('Wind residual [m/s]')
ss=sprintf('%s flights concatenated (%s PS)',DATE,PRESSURE)
title(ss)
grid
ss=sprintf('figs/Boxplot-factors-%s-PS.jpg',PRESSURE);
saveas(gcf,ss,'jpg')

cols = {'pit_offset' 'roll_offset' 'alpha_factor1' ...
    'alpha_factor2' 'beta_factor1' 'beta_factor2' 'head_offset' 'PS_fact'};
cols1 = strrep(cols, '_','\_');
ii = contains(cols,'offset')
facts = ii;

sz=size(beta_samp);
bx=beta_samp./mean(beta_samp,1);

nn = pppp*20;
figure(nn)
boxplot(bx,"label",cols,'notch','on');
title(sprintf('R858 factors divided by their mean (%s PS)',PRESSURE))
grid
ss=sprintf('figs/Boxplots_factors_all-%s-PS.jpg',PRESSURE);
saveas(gcf,ss,'jpg')

% Write out results 
fid=fopen(["./maneuvers_" + PRESSURE + ".txt"],'w')
fprintf(fid,                      'AWinds.Maneuver_date = 20260408');
fprintf(fid,            'AWinds.RollOffsetRadians= %g;\n',Params(1));
fprintf(fid,           'AWinds.PitchOffsetRadians= %g;\n',Params(2));
fprintf(fid,            'AWinds.HeadOffsetRadians= %g;\n',Params(3));
fprintf(fid,  'AWinds.AttackFactor= %g, %g ;\n',Params(4),Params(5));
fprintf(fid, 'AWinds.SideslipFactor= %g, %g;\n',Params(6),Params(7));
fprintf(fid,                'AWinds.PStaticOffset= %g;\n',Params(8));


fclose(fid)

xlsout2=["resultsTable_" + PRESSURE + ".xlsx"];

T_final=outputTable(structManeuvers,xlsout2);

end ; %pppp

