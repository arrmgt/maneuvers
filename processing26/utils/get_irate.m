function [varargout]=get_irate(ncfile, Var, varargin);
% function [irate,dims,frate]=get_irate(ncfile, Var[,"OutputRate",N]);
%
%  Input:  
%  ncfile           : *_raw.nc full path file name
%  Var              : Variable name
%  "OutputRate"     : name,value pair (optional)
%                   : e.g. ...,"OutputRate", 25,....
%                   : if "OutputRate" absent, output rate = input rate
%
%  Output:
%   irate           : input rate
%   dims            : dimensions
%   frate           : output rate if output rate is given and >= input rate
%   FillValue       : Fill value read from ncfile
%
%  If "OutputRate" not set, then irate = "raw" rate

% If raw input rate < processing output rate
%       then frate = raw input rate.
% Dimensions possible
%   nt = time dimension
%   sps = samples per second
%   ncell = number of cells 
try
    x = ncinfo(ncfile,Var);
    DIMS = cell2mat({x.Dimensions.Length});
    varargout{2} = DIMS;
    if numel(DIMS) == 3 
        ncell = DIMS(1);
        sps = DIMS(2);
        nt = DIMS(3);
    elseif numel(DIMS) == 2
        sps=DIMS(1);
        nt = DIMS(2);
    elseif numel(DIMS) == 1
        nt = DIMS;
        sps = 1;
    end
    irate = sps;
    varargout{1} = irate;    % Has a output rate been specified
    if ~isempty(varargin)
        p = inputParser;
        addParameter(p, 'OutputRate', [], @(x) isnumeric(x));
        parse(p, varargin{:});
        try
            orate = p.Results.OutputRate; % desired orate
        catch
            error('get_irate: unknown input %s',varargin)
        end
        % If orate was entered, change to frate if irate<orate
        if exist('orate','var') && irate <= orate 
            frate=irate;
            varargout{3}=frate;
        else
            varargout{3}=[];
        end
    end

    try
        FillValue = ncreadatt(ncfile,Var,'FillValue');
    catch
        FillValue = [];
    end
    varargout{4} = FillValue;

catch
    error(sprintf('get_irate: %s not available',Var))
    irate =[];
    dims = [];
    frate = [];
end

end