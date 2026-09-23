function y = changeRate(x, irate, orate)
% CHANGERATGE  Zero-phase resampling: decimation, upsampling, or fractional rate conversion.
%
%   y = changeRate(x, irate, orate)
%
%   x     : input signal or matrix [samples x channels]
%   irate : input sample rate (Hz)
%   orate : output sample rate (Hz)
%
%   orate > irate : upsamples via interp1 (pchip)
%   orate < irate : decimates in stages using Kaiser-windowed filtfilt
%                   Non-integer ratios upsample to LCM first, then decimate.
%   Each decimation stage uses: fc=fout/2, transWidth=fout/2, Ast=40 (~80dB via filtfilt).

%% --- Passthrough ---
if irate == orate
    y = x;
    return
end

% Quick check for NaN
x = fillmissing(x,'nearest');

%% --- Upsampling ---
if orate > irate
    %%%fprintf('changeRate: upsampling %g Hz --> %g Hz via pchip\n', irate, orate);
    y = upSample(x,irate,orate);
    return
end

%% --- Step 1: Upsample if needed to reach integer decimation ratio ---
lcm_rate = lcm(irate, orate);
up       = lcm_rate / irate;    % upsample factor (1 if already integer ratio)
R        = lcm_rate / orate;    % total integer decimation factor

if up > 1
    %%%fprintf('changeRate: upsampling by %d to %g Hz to achieve integer ratio\n', up, lcm_rate);
    x = upSample(x, irate,up*irate);

end

%% --- Step 2: Staged integer decimation ---
if R <= 1
    y = x;
    return
elseif R <= 10
    N = R;
    y = decimate_stage(x, N, orate);
    return
else
    N2 = plan_stages(R, 2);
    if max(N2) > 10
        N = plan_stages(R, 3);
    else
        N = N2;
    end
end

%%fprintf('changeRate: R=%d planned as [%s] (%d stage(s)), %g Hz --> %g Hz\n', ...
%%        R, num2str(N), length(N), lcm_rate, orate);

fout_stages = zeros(1, length(N));
for i = 1:length(N)
    fout_stages(i) = orate * prod(N(i+1:end));
end

y = x;
for i = 1:length(N)
    y = decimate_stage(y, N(i), fout_stages(i));
end

end

% =========================================================================
function stages = plan_stages(r, n)
% Factor r into n balanced stages, sorted ascending.
f      = sort(factor(r), 'descend');
stages = ones(1, n);
for i = 1:length(f)
    [~, idx]    = min(stages);
    stages(idx) = stages(idx) * f(i);
end
stages = sort(stages);

end

function y = upSample(x, irate, orate)
    isRow = isrow(x);
    if isRow, x = x(:); end
    nsamp = size(x, 1);
    t_in  = (1 : nsamp)' / irate;
    t_up  = (1 : ceil(nsamp * orate/irate))' / orate;
    t_up  = t_up(t_up <= t_in(end));   % stay within original time range
    y     = interp1(t_in, x, t_up, 'pchip');
    if isRow, y = y.'; end

end

% =========================================================================
function y = decimate_stage(x, r, fout)
% Single decimation stage: Kaiser lowpass via filtfilt, then downsample.

fs         = fout * r;
fc         = fout / 2;
transWidth = fout / 2;
Ap         = 0.1;
Ast        = 40;

b = kaiser_lp(fs, fc, transWidth, Ap, Ast);

N = length(b) - 1;
if size(x, 1) <= 3*N
    msg = sprintf('decimate1: signal length (%d) must be > 3x filter order (%d) for fout=%g Hz, r=%d.', ...
                  size(x,1), 3*N, fout, r);
    error(msg);
end

isRow = isrow(x);
if isRow, x = x(:); end

filtered = filtfilt(b, 1, x);
nout     = floor(size(x,1) / r);
y        = filtered(1 : r : nout*r, :);

if isRow, y = y.'; end

end


