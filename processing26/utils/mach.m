function [m,fact]=mach(varargin)
%function m=mach(qc,ps,[[mr],recovf]])
%MACH: Compute mach number from pressures
%$Source: /home/cvs/kingair/Sys09/mach.m,v $
%Project $Name:  $ ($Revision: 1.3 $)
%$Date: 2011/06/27 16:30:01 $
%

qc=varargin{1};
ps=varargin{2};
if(nargin>2)
  mr=varargin{3};
  if(isempty(mr)),mr=0;end;
  q=mr./(1+mr);
 else
  q=0;
end
if(nargin==4)
  recovf=varargin{4};
 else
  recovf=0.6425;
end

C=phycon();
Mv=C.Mv;
Md=C.Md;
Rd=C.Rd;
Cpd=C.Cpd;
Cvd=C.Cvd;
Cpv=C.Cpv;
Cvv=C.Cvv;
eps=C.eps;

R=Rd.*(1+q.*(1/eps-1));
Cp=Cpd.*(1+(Cpv/Cpd-1).*q);
Cv=Cvd.*(1+(Cvv/Cvd-1).*q);
k=R./Cp;
Gamma=Cp./Cv;

fact=(qc./ps+1).^k-1;

m2=2.*Cv./R.*fact;
m2(find(m2<0))=0;
m=sqrt(m2);

return

% Using symbolic toolbox
% syms Ua R Ta gmma r Tr Cp M Tt S Qc Ps k Cv
% Cp=Cv*gmma
%Ta=Tr/(r*((Qc/Ps + 1)^k - 1) + 1)
%Ua2=(2*Cv*Tr*gmma*((Qc/Ps + 1)^k - 1))/(r*((Qc/Ps + 1)^k - 1) + 1)
%S2=(R*Tr*gmma)/(r*((Qc/Ps + 1)^k - 1) + 1)
%M2=(2*Cv*((Qc/Ps + 1)^k - 1))/R
%%
%%       / / Qc     \k     \
%%  2 Cv | | -- + 1 |  - 1 |
%%       \ \ Ps     /      /
%%  ------------------------
%%             R


