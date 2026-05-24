function [tas,tasn]=tasf(varargin) ;               
%function [tas,tasn]=tasf(Qc,Ps,Ta[,mr]) ;               
%TASFM: Compute true airspeed from press/temps (moist version)
%$Source: /home/cvs/kingair/Sys09/tasf.m,v $
%Project $Name:  $ ($Revision: 1.3 $)
%$Date: 2012/05/07 15:41:52 $
%
%c                    
%c... input: 
%c...   Qc (pitot - static pressure, corrected, mb)          
%c...   Ps (static pressure, corrected, mb)            
%c...   Ta (static air temperature, K)             
%c...   mr (mixing ratio, kg/kg) -- optional for humidity correction
%c...                    
%c... output: tas (true airspeed, m/s) -- optionally humidity corrected              
%c...         tasn(true airspeed, m/s) -- optionally humidity corrected
%c...                                     using NCAR method
%c... 

Qc=varargin{1};
Ps=varargin{2};
Ta=varargin{3};
if(nargin>3)
   mr=varargin{4};
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
Cvd=C.Cvd;
Cpd=C.Cpd;
Cpv=C.Cpv;
eps=C.eps;

R=Rd.*(1+q.*(1/eps-1));
Cp=Cpd.*(1+(Cpv/Cpd-1).*q);
k=R./Cp;
	tas=sqrt(2.*Cp.*Ta.*((1.+Qc./Ps).^k-1.));
	tas(isnan(tas))=0;

%NCAR method (Bull 9, Apppendix B)
	M=mach(Qc,Ps);% dry
	Gamma=Cpd/Cvd;
        S=sqrt(Gamma.*Rd.*Ta);
	Tas=M.*S;
	tasn=Tas.*(1+.000304.*q.*1000);
