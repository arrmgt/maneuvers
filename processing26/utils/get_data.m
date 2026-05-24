function [y,ny]=get_data(varargin)
% x = get_data(filename,var,kk[,irate,orate])
%   > x = get_data(filename,var,kk,1000,10,false)
%  Inputs
%       filename    file name
%       var         variable name 
%       kk          indices of output []=all
%       irate       variable input rate (optional)
%       orate       variable output rate (optional)
%       flag        true if var is a heading variable
%   If irate and orate are missing, it returns all the data

fname=varargin{1};
vname=varargin{2};
n=length(varargin);
kk=[];
irate=[];
orate=[];
head=false;
ny=0;
y=[];

if n>2
    kk=varargin{3};
else
    kk=[];
end
if n>3
    irate=varargin{4};
    orate=varargin{5};
end
if n>5
    head=true;
end

try
  blurf=ncread(fname,vname);
  blurf=blurf(:);
catch 
  blurf=0;
  sprintf(sprintf('get_data: no data %s',vname))
  return
end

if(length(blurf)>1 & ~isempty(orate)),
    if(~head)
        y=change_rate(blurf,irate,orate);
    else
        s=sin(blurf.*pi/180);
        c=cos(blurf.*pi./180);
        s1=change_rate(s,irate,orate);
        c1=change_rate(c,irate,orate);
        blurf=atan2(s1,c1).*180/pi;
        y=blurf*180./pi;
    end
    if(~isempty(kk))
        y=y(kk);
    end
else
    y=blurf;
end
ny=length(y);

nn=find(~isnan(y));
y = interp1(nn,y(nn),[1:numel(y)]','linear',mean(y(nn)));

end