function T_final = outputTable(structManeuvers, outfile)
% do_it Consolidate structManeuvers into XLSX with interleaved CI/Values
%   do_it(structManeuvers, outfile)
%   structManeuvers: N×1 struct with fields fname, pressure, Params, CIlow, CIhi
%   outfile: filename for writetable, e.g., 'table1.xlsx'

if nargin < 2 || isempty(outfile)
    outfile = 'Maneuvers_Params.xlsx';
end

% Parameter names (8)
paramNames = {'Pitch','Roll','Alpha1','Alpha0','Beta1','Beta0','Head','Poff'};
nparam = numel(paramNames);

% Ensure column struct array
S = structManeuvers(:);
N = numel(S);

% Preallocate
fname = cell(N,1);
pressure = cell(N,1);
Params = nan(N,nparam);
CIlow  = nan(N,nparam);
CIhi   = nan(N,nparam);

% Convert each element to 1x8 numeric rows (pad/truncate as needed)
for i = 1:N
    % fname/pressure allow char/string/cell
    fname{i} = safeToChar(S(i).fname);
    pressure{i} = safeToChar(S(i).pressure);

    Params(i,:) = safeRow8(S(i).Params, i, 'Params');
    CIlow(i,:)  = safeRow8(S(i).CIlow,  i, 'CIlow');
    CIhi(i,:)   = safeRow8(S(i).CIhi,   i, 'CIhi');
end

% Build table pieces
CIlowNames = strcat('CIlow', paramNames);
CIhiNames  = strcat('CIhi',  paramNames);

T_base   = table(fname(:), pressure(:), 'VariableNames', {'fname','pressure'});
T_params = array2table(Params, 'VariableNames', paramNames);
T_CIlow  = array2table(CIlow,  'VariableNames', CIlowNames);
T_CIhi   = array2table(CIhi,   'VariableNames', CIhiNames);

% Combine and reorder to interleave per parameter: [CIlow_i, Param_i, CIhi_i]
T_combined = [T_base T_params T_CIlow T_CIhi];
interleaved = cell(1, 3*nparam);
ci = 1;
for k = 1:nparam
    interleaved{ci}   = CIlowNames{k};
    interleaved{ci+1} = paramNames{k};
    interleaved{ci+2} = CIhiNames{k};
    ci = ci + 3;
end
newOrder = [{'fname','pressure'}, interleaved];
T = T_combined(:, newOrder);

% Compute column-wise mean/std for numeric columns (skip fname, pressure)
numericVars = T.Properties.VariableNames(3:end);
dataMat = table2array(T(:, numericVars)); % N x M
colMean = mean(dataMat, 1, 'omitnan');
colStd  = std(dataMat, 0, 1, 'omitnan');

% Create mean and std rows
meanRow = array2table([[{'mean'},{''}], num2cell(colMean)], 'VariableNames', T.Properties.VariableNames);
stdRow  = array2table([[{'std' },{''}], num2cell(colStd) ], 'VariableNames', T.Properties.VariableNames);

% Ensure consistent variable classes across T, meanRow, stdRow (Method C)
tables = {T, meanRow, stdRow};
vars = tables{1}.Properties.VariableNames;

% Decide desired class per variable: prefer numeric if all tables can convert
desiredClass = cell(size(vars));
for vi = 1:numel(vars)
    name = vars{vi};
    canBeNumeric = true;
    for ti = 1:numel(tables)
        Ttmp = tables{ti};
        if ~ismember(name, Ttmp.Properties.VariableNames)
            canBeNumeric = false; break;
        end
        col = Ttmp.(name);
        if isnumeric(col)
            continue;
        elseif iscell(col) && all(cellfun(@(x) isnumeric(x) && isscalar(x), col))
            continue;
        else
            canBeNumeric = false; break;
        end
    end
    if canBeNumeric
        desiredClass{vi} = 'numeric';
    else
        desiredClass{vi} = 'cell';
    end
end

% Convert each table's variables to the desired class
for ti = 1:numel(tables)
    Ttmp = tables{ti};
    for vi = 1:numel(vars)
        name = vars{vi};
        if ~ismember(name, Ttmp.Properties.VariableNames)
            error('Table %d is missing variable %s', ti, name);
        end
        if strcmp(desiredClass{vi}, 'numeric')
            if iscell(Ttmp.(name))
                col = Ttmp.(name);
                if all(cellfun(@(x) isnumeric(x) && isscalar(x), col))
                    Ttmp.(name) = cellfun(@(x) double(x), col);
                else
                    error('Cannot convert variable %s in table %d to numeric.', name, ti);
                end
            end
            % numeric already OK
        else % desired class = cell
            if isnumeric(Ttmp.(name))
                Ttmp.(name) = num2cell(Ttmp.(name));
            end
            % cell already OK
        end
    end
    tables{ti} = Ttmp;
end

% Vertically concatenate and write
T_final = vertcat(tables{:});
writetable(T_final, outfile);

fprintf('Wrote %d maneuvers + 2 summary rows to %s\n', N, outfile);

%% Helper subfunctions

    function s = safeToChar(x)
        % Convert fname/pressure to char (or cellstr) for table cell column
        if ischar(x)
            s = x;
        elseif isstring(x)
            s = char(x);
        elseif iscell(x) && numel(x)==1 && (ischar(x{1}) || isstring(x{1}))
            s = char(x{1});
        else
            % Fallback: convert whatever to char short representation
            try
                s = char(string(x));
            catch
                s = '';
            end
        end
    end

    function v = safeRow8(x, idx, name)
        % Convert various types into 1x8 numeric row. Pad/truncate as needed.
        % Returns 1x8 double row (NaN when padded). Warns on empty or conversion events.
        % idx/name used for warning/error messages.
        % Unwrap 1x1 cell
        if iscell(x) && numel(x) == 1
            x = x{1};
        end
        % Table with one row -> array
        if istable(x) && height(x) == 1
            try, x = table2array(x); catch, end
        end
        % String/char that might encode numbers
        if ischar(x) || (isstring(x) && isscalar(x))
            nums = str2num(char(x)); %#ok<ST2NM>
            if ~isempty(nums)
                x = nums;
            end
        end
        % If still cell of numerics, try to mat
        if iscell(x)
            try, x = cell2mat(x); catch, end
        end

        % If empty produce NaNs
        if isempty(x)
            warning('Element %d: field %s is empty — padding with NaNs.', idx, name);
            v = nan(1, nparam);
            return
        end

        % Flatten to row
        try
            xr = double(x(:).'); % convert to numeric if possible
        catch
            error('Element %d: field %s could not be converted to numeric.', idx, name);
        end

        % Pad or truncate
        if numel(xr) < nparam
            warning('Element %d: field %s length %d < %d — padding with NaNs.', idx, name, numel(xr), nparam);
            xr = [xr, nan(1, nparam - numel(xr))];
        elseif numel(xr) > nparam
            warning('Element %d: field %s length %d > %d — truncating to %d.', idx, name, numel(xr), nparam, nparam);
            xr = xr(1:nparam);
        end
        v = xr;
    end

end
