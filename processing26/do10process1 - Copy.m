clear;close all
% Restore factory-default MATLAB path for this session
restoredefaultpath
topdir = 'c:/users/rodi/Github/maneuvers/processing26';
addpath(fullfile(topdir,'utils'))
close all
nnn=0;
clear structManeuvers  M
;
FLTs = ["20260408a_arr.c10.nc" "20260408b_arr.c10.nc"]
irate = 10;
orate = 10;
;
M.flt(1,1).kk = [38638:48126];
M.flt(1,2).kk = [75808:84509];
M.flt(2,1).kk = [22720:31776];
M.flt(2,2).kk = [];

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

PRESSURE = 'SHIP';
for jj = 1:numel(FLTs)
    dataDir = fullfile('P:/MATLAB-DATA2/kingair_data/',PROJ,'work');
    file = fullfile(dataDir,FLTs(jj));
    ARM = ncreadatt(file,'/','AWinds.MomentArm');
   
    for ii=1:mm;
        clear kk
        ss = sprintf("kk = M.flt(%i,%i).kk;",jj,ii);
        eval(ss);
        if isempty(kk)
            continue
        end

        % Get the data
        % BOOM or SHIP depends on get_varTAS.m
        tas = get_data(file,'TASX',kk,irate,orate); 
        psm = get_data(file,'PSX',kk,irate,orate);
        dp1 = get_data(file,'DP1X',kk,irate,orate);
        q_impact = get_data(file,'QimpactX',kk,irate,orate);
        alpha = get_data(file,'alpha',kk,irate,orate);
        beta = get_data(file,'beta',kk,irate,orate);
        temp = get_data(file,'TEMPX',kk,irate,orate); % K
        dpa = get_data(file,'dpa',kk,irate,orate);
        dpb = get_data(file,'dpb',kk,irate,orate);
        dpr = get_data(file,'dpr_beta',kk,irate,orate);
        
        ias = get_data(file,'ias',kk,irate,orate);
        mr= zeros(size(dpr));

        pcorc = zeros(size(dpr));
        pmb = psm-pcorc;
        
        % mixing ratio wasn't always archived -- but set it to zero for
        % now.
        if(numel(mr)<=1)
            mr=zeros(size(tas));
        end
        
        npitch1=[];npitch2=[];npitch3=[];
        [~,npitch1] = get_data(file,'AVpitch',kk,irate,orate);
        [~,npitch2] = get_data(file,'avpitch',kk,irate,orate);
        
        % Use the post-processed applanix values if available        
        if(npitch1>1)
            'AV variables used'
            pitch = get_data(file,'AVpitch',kk,irate,orate);
            roll = get_data(file,'AVroll',kk,irate,orate);
            thead = get_data(file,'AVthead',kk,irate,orate,true);
            pitchr = get_data(file,'AVpitchr',kk,irate,orate);
            rollr = get_data(file,'AVrollr',kk,irate,orate);
            yawr = get_data(file,'AVyawr',kk,irate,orate);
            ewvel = get_data(file,'AVewvel',kk,irate,orate);
            nsvel = get_data(file,'AVnsvel',kk,irate,orate);
            zvel = get_data(file,'AVzvel',kk,irate,orate);
            norma = get_data(file,'AVnorma',kk,irate,orate);
            lata = get_data(file,'AVlata',kk,irate,orate);
            longa = get_data(file,'AVlonga',kk,irate,orate);
        elseif(npitch2>1)
            'av variables used'
            pitch = get_data(file,'avpitch',kk,irate,orate);
            roll = get_data(file,'avroll',kk,irate,orate);
            thead = get_data(file,'avthead',kk,irate,orate,true);
            pitchr = get_data(file,'avpitchr',kk,irate,orate);
            rollr = get_data(file,'avrollr',kk,irate,orate);
            yawr = get_data(file,'avyawr',kk,irate,orate);
            ewvel = get_data(file,'avewvel',kk,irate,orate);
            nsvel = get_data(file,'avnsvel',kk,irate,orate);
            zvel = get_data(file,'avzvel',kk,irate,orate);
            norma = get_data(file,'avnorma',kk,irate,orate);
            lata = get_data(file,'avlata',kk,irate,orate);
            longa = get_data(file,'avlonga',kk,irate,orate);
        end
        % save indices for this flight and run the regression
        %  to find the Params individually for the flights.
        KKKS0=1;
        KKKS1=numel(kk);
        [Params,fu,fv,fw,mag,dir,resid,jacobian,CI]= ...
        do_fits0(KKKS0,KKKS1,dp1,dpb,dpa,dpr,q_impact,temp,tas,pmb,pcorc,mr,alpha,beta ...
        ,roll,pitch,thead,rollr,pitchr,yawr,ewvel,nsvel,zvel,ARM);

        % Save the indices 
        bb = bb  + 1;
        kkks0(bb)=iiis0;
        kkks1(bb)=iiis0+length(kk)-1;
        fname{bb}=FLTs(jj);
        iiis0=iiis0+length(kk);
        
        %  Save things in the structureManeuvers
        nnn=nnn+1; % row number
        structManeuvers(nnn,1).fname=FLTs(jj);
        structManeuvers(nnn,1).pressure=round(mean(pmb));
        structManeuvers(nnn,1).Params=Params;
        structManeuvers(nnn,1).CIlow=CI(:,1)';
        structManeuvers(nnn,1).CIhi=CI(:,2)';

        dp1_1=[dp1_1;dp1];
        alpha1=[alpha1;alpha];
        beta1=[beta1;beta];
        temp1=[temp1;temp];
        tas1=[tas1;tas];
        dpa1=[dpa1;dpa];
        dpb1=[dpb1;dpb];
        dpr1=[dpr1;dpr];
        dp11=[dp11;dp1];
        pmb1=[pmb1;pmb];
        pcorc1=[pcorc1;pcorc];
        roll1=[roll1;roll];
        pitch1=[pitch1;pitch];
        thead1=[thead1;thead];
        rollr1=[rollr1;rollr];
        pitchr1=[pitchr1;pitchr];
        yawr1=[yawr1;yawr];
        vew1=[vew1;ewvel];
        vns1=[vns1;nsvel];
        vz1=[vz1;zvel];
        norma1=[norma1;norma];
        lata1=[lata1;lata];
        longa1=[longa1;longa];
        q_impact1=[q_impact1;q_impact];
    end
end
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

% Moment arm from IMU to boom tip
ARM = [6.77,  0.50,   -0.44] ;
[Params,fu,fv,fw,mag,dir,resid,jacobian,CI,beta_samp]= ...
    do_fits0(kkks0,kkks1,dp1,dpb,dpa,dpr,q_impact,temp,tas,pmb,pcorc,mr ...
    ,alpha,beta,roll,pitch,thead,rollr,pitchr,yawr,vew,vns,vz,ARM);
f=[fu,fv,fw];

x1={structManeuvers.Params};
y1={structManeuvers.CIlow};
z1={structManeuvers.CIhi};
x2={structManeuvers.Params};
y2={structManeuvers.CIlow};
z2={structManeuvers.CIhi};

nnn=nnn+2;
structManeuvers(nnn,1).fname='means';
structManeuvers(nnn,1).Params=mean(cell2mat(x1(:)));
structManeuvers(nnn,1).CIlow=mean(cell2mat(y1(:)));
structManeuvers(nnn,1).CIhi=mean(cell2mat(z1(:))); 

nnn=nnn+1;
structManeuvers(nnn,1).fname='stdevs';
structManeuvers(nnn,1).Params=std(cell2mat(x2(:)));
structManeuvers(nnn,1).CIlow=std(cell2mat(y2(:)));
structManeuvers(nnn,1).CIhi=std(cell2mat(z2(:)));

nnn=nnn+2; % row number
structManeuvers(nnn,1).fname='concat';
structManeuvers(nnn,1).Params=Params;
structManeuvers(nnn,1).CIlow=CI(:,1)';
structManeuvers(nnn,1).CIhi=CI(:,2)';

nnn=nnn+1; 
structManeuvers(nnn,1).fname='CI-low';
structManeuvers(nnn,1).Params=CI(:,1)';
nnn=nnn+1; 
structManeuvers(nnn,1).fname='CI-hi';
structManeuvers(nnn,1).Params=CI(:,2)';
nnn=nnn+1; 
structManeuvers(nnn,1).fname='CI-diff';
structManeuvers(nnn,1).Params=CI(:,2)'-CI(:,1)';


%  Save the results in a spreadsheet
xlsout='./results_all.xlsx';
delete(xlsout)
clear Tnew
Tnew= [struct2table(structManeuvers)];
writetable(Tnew,xlsout,'Sheet',1,'Range','A1')


% Plot results
delete figs/*.jpg
figure(1)
h=plot((1:length(fu))/orate,detrend([fw,fu,fv],'constant'))
flight='Concat';
jju=1:length(tas);
set(h(1),'LineWidth',1.5);
set(h(2),'LineWidth',1.5);
set(h(3),'LineWidth',1.5);
ss=sprintf('%i flights concatenated',numel(fname))
title(ss)
xlabel('Time [secs]')
ylabel('Wind component [m/s]')
legend('up','east','north','location','northeast')
grid
print('figs/Detrended-comps.jpg','-djpeg100')


figure(2)
h=plot((1:length(jju))/orate,[alpha(jju),beta(jju)].*180./pi)
set(h,'LineWidth',1)
xlabel('Relative time [sec]')
ylabel('Flow angle [deg]')
ss=sprintf('%i flights concatenated',numel(fname))
title(ss)
legend('Attack angle','Sideslip angle','location','northwest')
grid
print('figs/Flow angles.jpg','-djpeg100');

figure(3)
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
ss=sprintf('%i flights concatenated',numel(fname))
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
print('figs/Resids_hist.jpg','-djpeg100');

figure(4)
boxplot([fu,fv,fw],'notch','on','whisker',1,'Symbol','.' ...
    ,'Labels',{'East-component','North-component','Up-component'})
ylabel('Wind residual [m/s]')
ss=sprintf('%i flights concatenated',numel(fname))
title(ss)
grid
print('figs/Boxplot.jpg','-djpeg100');

cols = {'pit_offset' 'roll_offset' 'alpha_factor1' ...
    'alpha_factor2' 'beta_factor1' 'beta_factor' 'head_offset' 'Q_fact'};
cols1 = strrep(cols, '_','\_');
ii = contains(cols,'offset')
facts = ii;

sz=size(beta_samp);
bx=beta_samp./mean(beta_samp,1);

for i=1:sz(2);
    figure(80+i)
    boxplot(beta_samp(:,i),"label",cols{i},'notch','on');
    ylabel(sprintf('%s',cols1{i}));
    title(sprintf('Maneuver-derived &s factor',cols{i}))
    if(facts(i))
        ylabel(sprintf('%s [degrees]',cols1{i}));
    else
        ylabel(sprintf('%s [multiplier]',cols1{i}))
    end
    ss=sprintf("print('figs/Boxplot-%s.jpg','-djpeg100')",cols{i})
    eval(ss)
end
figure(80+i+1)
boxplot(bx,"label",cols,'notch','on');
title('Boom factors (divided by their mean)')
grid
print('figs/Boxplots_factors_all.jpg','-djpeg100');

% Write out results 
fid=fopen('./maneuvers.txt','w')

fprintf(fid,'                :AWinds.QFactor= %g;\n',Params(8));
fprintf(fid,'                :AWinds.AttackFactor= %g, %g ;\n',Params(3),Params(4));
fprintf(fid,'                :AWinds.SideslipFactor= %g, %g;\n',Params(5),Params(6));
fprintf(fid,'                :AWinds.PitchOffsetRadians= %g;\n',Params(1));
fprintf(fid,'                :AWinds.RollOffsetRadians= %g;\n',Params(2));
fprintf(fid,'                :AWinds.HeadOffsetRadians= %g;\n',Params(7));

fclose(fid)

T_final=outputTable(structManeuvers,'resultsTable.xlsx');