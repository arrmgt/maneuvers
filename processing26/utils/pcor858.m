function pcor = pcor858(fcoef,DP1,DPA,DPB,DPR,PSA)
[~,fcoef] = cone_pcor(DP1,DPB,DPA,DPR,PSA); %cone_pcor(dp1,pb,pa,pr,psm)
fq1 = fqCalc(DPA,DPB,DPR); % fqCalc(dpa,dpb,dpr)
pcor = pcor_beta(DP1,DPA,DPB,DPR,fcoef,fq1); %pcor_beta(dp1a,dpa,dpb,dpr,f,fq)
