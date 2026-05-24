function [ts,Ts]=tstatic(varargin)
%function ts=tstatic    
%TSTATIC: compute static air temp from total temp and pressures (moist version)
%function ts=tstatic(Tm,recovf,Qc,Ps,mr)    
%$Source: /home/cvs/kingair/Sys09/tstatic.m,v $
%Project $Name:  $ ($Revision: 1.3 $)
%$Date: 2012/05/07 15:41:52 $
%c..     
%c.. inputs:  
%c.. Tm (measured temperature, K)
%c.. recovf (recovery factor)    
%c.. Qc     (pitot-static pressure, corrected, mb)     
%c.. Ps     (static pressure, corrected, mb)  
%c.. mr     (mixing ratio, g/g)
%c..     
%c..  output:  tstatic (static air temperature, K)     
%c..    

Tm=varargin{1};
recovf=varargin{2};
Qc=varargin{3};
Ps=varargin{4};
if(length(varargin)>4)
   mr=varargin{5};
   if(isempty(mr)),mr=0;end
   q=mr./(1+mr);
 else
   mr=0;
   q=0;
end

C=phycon();
Mv=C.Mv;
Md=C.Md;
Rd=C.Rd;
Cpd=C.Cpd;
Cvd=C.Cvd;
Cpv=C.Cpv;
eps=C.eps;

R=Rd.*(1+q.*(1/eps-1));
Cp=Cpd.*(1+(Cpv/Cpd-1).*q);
k=R./Cp;
	ts = Tm  ./(1.+recovf.*((1.+Qc./Ps).^k-1.)); 
	ii=find(Qc < 5.);ts(ii)=Tm(ii);   
%NCAR method (Bull 9, Apppendix B)
        M=mach(Qc,Ps);% dry
        Gamma=Cpd/Cvd;
	Ts= Tm./(1+recovf.*M.^2.*(Gamma-1.)/2);
 
