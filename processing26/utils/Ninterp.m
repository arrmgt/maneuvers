function [out,nnan]=Ninterp(varargin)
%NINTERP: special function to orchestrate linear interpolation
%$Id: Ninterp.m,v 1.2 2012/02/20 13:53:13 rodi Exp $
%function [out,nnan]=Ninterp(in,ninterp)
%$Source: /home/cvs/kingair/Sys09/Ninterp.m,v $
%Project $Name: trans2am21_qc0 $ ($Revision: 1.2 $)
%$Date: 2012/02/20 13:53:13 $

in=varargin{1};
ninterp=varargin{2};
if(ninterp==1);
    out=in;
    return
end

if(nargin>=3)
  METHOD=varargin{3};
 else
  METHOD='linear';
end

if(nargin==4);
  FillValue=varargin{4};
else
  FillValue=-32767; % default
end

jj=find(isnan(in));
in(jj)=FillValue;% just in case FillValue~=NaN;
kk=find(in~=FillValue); % NOT FillValue
if(~isempty(kk))
    in=in(:);
    nn=length(in);
    iii=1:1/ninterp:nn+1-1/ninterp;
    out=interp1(kk,in(kk),iii,METHOD,FillValue)';
else
    out=interp(in,ninterp);
end

return
