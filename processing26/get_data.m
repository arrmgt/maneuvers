function [y,ny] = get_data(varargin)
% x = get_data(filename,var,kk[,irate,orate,flag])
%  Inputs
%       filename    file name
%       var         variable name 
%       kk          indices of output []=all
%       irate       variable input rate (optional)
%       orate       variable output rate (optional)
%       flag        heading? (optional)
% Validate the input filename and variable name

fname=varargin{1};
vname=varargin{2};
n=length(varargin);
kk=[];
irate=[];
orate = [];
head =[];
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

if n > 5
    head = varargin{5};
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
    if numel(head)==0 
        y=change_rate(blurf,irate,orate);
    else
        fact = 1;
        if any(y > 2*pi)
            fact = pi/180; % degrees
        end
        s=sin(blurf.*fact);
        c=cos(blurf.*fact);
        s1=change_rate(s,irate,orate);
        c1=change_rate(c,irate,orate);
        blurf=atan2(s1,c1).*180/pi;
        y=blurf/fact;
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


