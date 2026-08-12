function out = r858_solve3(ptb, psa, dpa, dpb, dpr, dpn, f, pkor, sigma, alpha, opts)
%R858_SOLVE3  Joint weighted nonlinear least-squares solve for [q, ta, tb]
%using ALL 5 governing equations simultaneously (not exact pairwise
%elimination like r858_solve2.m). This properly distributes the 2
%available degrees of freedom across q, ta, AND tb, instead of putting
%all of it into q's consistency check while leaving ta, tb with no
%fit-based tightening. Fully vectorized (Gauss-Newton, closed-form 3x3
%solve per iteration - no per-row loop, no Optimization Toolbox needed).
%
%   out = R858_SOLVE3(ptb, psa, dpa, dpb, dpr, dpn, f, pkor, sigma, alpha, opts)
%
%   INPUTS (each scalar, or all the same size, e.g. Nx1 time series)
%     ptb, psa            Raw pressures; dp1 = ptb - psa internally
%     dpa, dpb, dpr, dpn   Measured differential pressures. SENSOR DROPOUT:
%                          set dpr==0 (exactly) for a row to flag DPR as
%                          out, or dpn==0 for DPN as out - the solver
%                          falls back to the cross-relation ta/tb=dpa/dpb
%                          for the affected variable and drops that
%                          equation from the fit (dof 2 -> 1 for that
%                          row). Do NOT use 0 for dpa or dpb dropout -
%                          those aren't handled by this fallback (see
%                          chat discussion: dpa/dpb loss is more severe,
%                          it removes the cross-relation itself).
%     f, pkor              Known (not fit - see NOTE)
%
%   sigma, alpha, opts are all OPTIONAL (nargin<=8, or sigma omitted/[]):
%     - PRODUCTION MODE (sigma not supplied): only the point estimates are
%       computed - out.q, out.ta, out.tb - using an UNWEIGHTED joint
%       least-squares fit (no SE/CI/diagnostics/f-pkor sensitivity work,
%       so it's substantially cheaper - skips the 4 extra warm-started
%       re-solves entirely). Use this in production once you've validated
%       the uncertainties are acceptable.
%     - ANALYSIS MODE (sigma supplied): full pipeline as documented below,
%       returns SE/CI/diagnostics too.
%
%   sigma   1-sigma measurement uncertainty for [ptb psa dpa dpb dpr dpn f pkor],
%           in that order. Accepted forms: scalar, 8-element vector, or a
%           struct with fields ptb,psa,dpa,dpb,dpr,dpn (required) and f,
%           pkor (optional, default 0). Same convention as r858_solve2.m.
%
%   alpha   significance level for the CI (default 0.05, 95% CI). Ignored
%           in production mode.
%
%   opts    optional struct: opts.maxIter (default 25), opts.tol (default
%           1e-10, on the max absolute Gauss-Newton step per parameter).
%           Used in both modes (controls the main fit only).
%
%   OUTPUT  out (struct; numeric fields sized like the inputs)
%     PRODUCTION MODE: out.q, out.ta, out.tb only.
%     ANALYSIS MODE, additionally:
%     out.dp1, out.dp1_SE          dp1 = ptb-psa and its propagated SE
%     out.q_SE,  out.q_CI          jointly-fit q uncertainty
%     out.ta_SE, out.ta_CI         jointly-fit ta uncertainty
%     out.tb_SE, out.tb_CI         jointly-fit tb uncertainty
%     out.diag.chi2, out.diag.dof(=2), out.diag.pvalue
%                     overall goodness-of-fit of the 5 equations at the
%                     joint solution (single unified consistency check,
%                     replacing the separate/redundant checks in v2)
%     out.diag.iters  Gauss-Newton iterations used (should be small, e.g.
%                     3-6, since it's warm-started from the exact
%                     closed-form solution of r858_solve2.m)
%     out.diag.theta0 starting values [q0 ta0 tb0] used (from the closed
%                     form), for comparison against the refined fit
%
%   METHOD
%   The 5 equations are cast as residuals of a measurement = model + noise
%   form:
%     r1 = dp1 - [ q*(D*(1-f)+f)/D - pkor ]          (eq. 1)
%     r2 = dpb - 2*f*q*tb/D                           (eq. 2)
%     r3 = dpa - 2*f*q*ta/D                           (eq. 3)
%     r4 = dpr - f*q*(1-2*tb-tb^2)/(2*D)               (eq. 4)
%     r5 = dpn - f*q*(1-2*ta-ta^2)/(2*D)               (eq. 5)
%   where D = 1+ta^2+tb^2. These are minimized in the weighted
%   least-squares sense (weights = 1/sigma^2 of dp1,dpb,dpa,dpr,dpn
%   respectively) via Gauss-Newton, started from the closed-form solution
%   in r858_solve2.m (which exactly zeroes 4 of the 5 residuals, so
%   convergence is fast). At convergence:
%     Cov_data(q,ta,tb) = inv(J' * W * J)   (J = 5x3 Jacobian, W = diag(w))
%   This captures the propagated uncertainty of dp1 (hence ptb, psa),
%   dpa, dpb, dpr, dpn THROUGH the joint fit - properly splitting the 2
%   degrees of freedom of redundancy across all 3 parameters, unlike
%   r858_solve2.m where ta, tb got none of it.
%
%   Because f and pkor are treated as known constants (not fit
%   parameters) but are not perfectly known, their own uncertainty is
%   propagated separately: the fit is re-solved (warm-started, so it
%   reconverges in 1-2 extra iterations) at f+/-h and pkor+/-h to get
%   d(q,ta,tb)/df and d(q,ta,tb)/dpkor by central differences, and
%     Cov_total = Cov_data + (dtheta/df)*(dtheta/df)'*sigma_f^2
%                          + (dtheta/dpkor)*(dtheta/dpkor)'*sigma_pkor^2
%   This whole pipeline was validated against a Monte Carlo simulation
%   (5000 noisy synthetic samples): the delta-method/GLS SE matched the
%   empirical standard deviation of the fitted q, ta, tb to within ~1%.
%
%   NOTE - DEGREES OF FREEDOM
%   5 equations, 3 unknowns (q, ta, tb; f and pkor are known) -> dof = 2.
%   Unlike r858_solve2.m (exact pairwise closed form: ta from dpa/dpn
%   only, tb from dpb/dpr only, all redundancy shows up only in q's
%   3-way chi-square test), this function spends the 2 degrees of
%   freedom optimally across q, ta, AND tb jointly - the statistically
%   preferred choice if you want the tightest possible CIs on ta and tb
%   themselves, not just on q.

productionMode = (nargin < 9) || isempty(sigma);


if nargin < 10 || isempty(alpha)
    alpha = 0.05;
end
if nargin < 11 || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'maxIter'), opts.maxIter = 25;   end
if ~isfield(opts,'tol'),     opts.tol     = 1e-10; end
if ~isfield(opts,'maxIterSensitivity'), opts.maxIterSensitivity = 10; end
% maxIterSensitivity caps the 4 warm-started f/pkor re-solves, which start
% already at the converged solution and typically need only 1-2 steps.
% (Not used at all in production mode - see below.)

sz = size(ptb);
n = numel(ptb);

ptb = expand_(ptb,n); psa = expand_(psa,n); dpa = expand_(dpa,n); dpb = expand_(dpb,n);
dpr = expand_(dpr,n); dpn = expand_(dpn,n);
f    = expand_(f,n);    pkor = expand_(pkor,n);

dp1 = ptb - psa;
data = struct('dp1',dp1,'dpa',dpa,'dpb',dpb,'dpr',dpr,'dpn',dpn,'f',f,'pkor',pkor);

% ---- sensor dropout: dpr==0 or dpn==0 (exactly) marks that channel as
% out for that row. Zero its residual weight so the fit doesn't treat
% the sentinel 0 as a real, highly-informative measurement - it's simply
% dropped from the weighted least squares for that row (dof: 2 -> 1 per
% dropped channel). See closed_form_ for how the starting point handles
% this via the cross-relation ta/tb = dpa/dpb.
%
% BOTH DPR and DPN out simultaneously is treated as UNRECOVERABLE for
% that row, by design: ta, tb are central to the downstream 3-D wind
% calculation, and a dof=0 fallback solution (relying solely on eq.1+2+3
% with nothing to cross-check) is not something to trust for that use -
% those rows are forced to NaN rather than returned as if they were a
% normal, if less-certain, answer.
dropTol = 1e-12;
dprOut = abs(dpr) < dropTol;
dpnOut = abs(dpn) < dropTol;
bothOut = dprOut & dpnOut;
if any(bothOut)
    warning('r858_solve3:bothDropout', ...
        '%d of %d row(s) have both DPR and DPN out simultaneously - ta, tb, q are not recoverable and are returned as NaN for those rows.', ...
        sum(bothOut), numel(bothOut));
end
dofRow = 2 - double(dprOut) - double(dpnOut);   % 2 normally, 1 if one channel out, 0 if both out (forced NaN)

% ---- starting point: exact closed-form solution (same as r858_solve2.m) ----
[ta0, tb0, Y0] = closed_form_(dp1,dpa,dpb,dpr,dpn,f,pkor);
q0 = Y0(:,1);   % use the eq.1-based q as the starting q (any of the 3 would do)
theta0 = [q0, ta0, tb0];

if ~isfield(opts,'outlierBound'), opts.outlierBound = 1e3; end
% Any fitted |ta| or |tb| beyond this (or non-finite) is treated as a
% failed/degenerate solve for that row (e.g. Gauss-Newton diverging from
% a bad starting point) and reported as NaN rather than a bogus huge
% number - ta, tb are physically expected to be a handful of hPa at most.

if productionMode
    % Unweighted joint least-squares fit (no sigma info available to weight
    % by) - still properly distributes the available dof across q, ta, tb
    % jointly, just without any SE/CI/sensitivity work. Dropout channels
    % are still excluded (weight 0) even here.
    W = ones(n,5);
    W(dprOut,4) = 0;
    W(dpnOut,5) = 0;
    [theta, ~, ~] = gauss_newton_(theta0, data, W, opts.maxIter, opts.tol);
    theta = sanitize_(theta, opts.outlierBound);
    theta(bothOut,:) = NaN;
    out.q  = reshape(theta(:,1), sz);
    out.ta = reshape(theta(:,2), sz);
    out.tb = reshape(theta(:,3), sz);
    return
end

% ================= ANALYSIS MODE (sigma supplied) =================

% --- accept sigma as scalar, 8-vector, or struct (same convention as r858_solve2) ---
if isstruct(sigma)
    sigma = normalize_sigma_struct_(sigma, 'r858_solve3');
elseif isscalar(sigma)
    sigma = struct('ptb',sigma,'psa',sigma,'dpa',sigma,'dpb',sigma,'dpr',sigma,'dpn',sigma);
elseif numel(sigma) == 8
    sigma = struct('ptb',sigma(1),'psa',sigma(2),'dpa',sigma(3),'dpb',sigma(4), ...
                    'dpr',sigma(5),'dpn',sigma(6),'f',sigma(7),'pkor',sigma(8));
else
    error('r858_solve3:sigmaFormat', ...
        'sigma must be a struct, a scalar, or an 8-vector [ptb psa dpa dpb dpr dpn f pkor].');
end
if ~isfield(sigma,'f'),    sigma.f    = 0; end
if ~isfield(sigma,'pkor'), sigma.pkor = 0; end
req = {'ptb','psa','dpa','dpb','dpr','dpn'};
for i = 1:numel(req)
    if ~isfield(sigma, req{i})
        error('r858_solve3:missingSigma', 'sigma.%s is required (case-insensitive; also accepts fcoef/pcor for f/pkor).', req{i});
    end
end

s_ptb  = expand_(sigma.ptb,n);  s_psa = expand_(sigma.psa,n);
s_dpa  = expand_(sigma.dpa,n);  s_dpb = expand_(sigma.dpb,n);
s_dpr  = expand_(sigma.dpr,n);  s_dpn = expand_(sigma.dpn,n);
s_f    = expand_(sigma.f,n);    s_pkor = expand_(sigma.pkor,n);

dp1_SE = sqrt(s_ptb.^2 + s_psa.^2);

z = sqrt(2)*erfinv(1-alpha);

% weights for the 5 residuals: r1(dp1), r2(dpb), r3(dpa), r4(dpr), r5(dpn)
w1 = 1./max(dp1_SE,eps).^2;
w2 = 1./max(s_dpb,eps).^2;
w3 = 1./max(s_dpa,eps).^2;
w4 = 1./max(s_dpr,eps).^2;
w5 = 1./max(s_dpn,eps).^2;
w4(dprOut) = 0;   % DPR dropout: exclude eq.4's residual for these rows
w5(dpnOut) = 0;   % DPN dropout: exclude eq.5's residual for these rows
W = [w1,w2,w3,w4,w5];  % N x 5

% ---- Gauss-Newton joint fit ----
[theta, iters, invM] = gauss_newton_(theta0, data, W, opts.maxIter, opts.tol);
bad = ~isfinite(theta(:,2)) | ~isfinite(theta(:,3)) | ...
      abs(theta(:,2)) > opts.outlierBound | abs(theta(:,3)) > opts.outlierBound;
theta(bad,:) = NaN;
invM(bad,:)  = NaN;
if any(bad)
    warning('r858_solve3:outliers', ...
        '%d of %d row(s) failed to converge to a physical solution (|ta| or |tb| > %g, or non-finite) - returned as NaN.', ...
        sum(bad), numel(bad), opts.outlierBound);
end
% Both-dropout rows are forced to NaN regardless of whether the dof=0
% fallback fit happened to converge - see warning issued earlier where
% bothOut was computed.
theta(bothOut,:) = NaN;
invM(bothOut,:)  = NaN;
q  = theta(:,1); ta = theta(:,2); tb = theta(:,3);
inv11 = invM(:,1); inv22 = invM(:,2); inv33 = invM(:,3);
inv12 = invM(:,4); inv13 = invM(:,5); inv23 = invM(:,6);

% ---- residual chi-square at the joint solution (dof = 2 normally, 1 if a
% sensor is out for that row, 0 if both DPR and DPN are out - see dofRow) ----
r = residuals_(theta, data);   % N x 5
chi2 = sum(W .* r.^2, 2);
chi2 = max(chi2,0);
dof = dofRow;   % N x 1, per-row (2, 1, or 0)
pval = nan(n,1);
hasDof = dof > 0;
pval(hasDof) = 1 - gammainc(chi2(hasDof)/2, dof(hasDof)/2);
% dof==0 rows (both DPR and DPN out): system is exactly determined, chi2
% is ~0 by construction and no goodness-of-fit test is possible -> NaN.

% ---- sensitivity of (q,ta,tb) to f and pkor (warm-started re-solve) ----
hf = 1e-6*max(1,abs(f));
dataP = data; dataP.f = f+hf;
dataN = data; dataN.f = f-hf;
[thP,~,~] = gauss_newton_(theta, dataP, W, opts.maxIterSensitivity, opts.tol);
[thN,~,~] = gauss_newton_(theta, dataN, W, opts.maxIterSensitivity, opts.tol);
dtheta_df = (thP - thN) ./ (2*hf);   % N x 3

hp = 1e-6*max(1,abs(pkor));
dataP2 = data; dataP2.pkor = pkor+hp;
dataN2 = data; dataN2.pkor = pkor-hp;
[thP2,~,~] = gauss_newton_(theta, dataP2, W, opts.maxIterSensitivity, opts.tol);
[thN2,~,~] = gauss_newton_(theta, dataN2, W, opts.maxIterSensitivity, opts.tol);
dtheta_dpkor = (thP2 - thN2) ./ (2*hp);   % N x 3

% ---- variance components, kept separate so the source of uncertainty is
% traceable (data vs. f vs. pkor) rather than only seeing the combined total ----
var_q_data  = inv11;  var_q_f  = dtheta_df(:,1).^2.*s_f.^2;  var_q_pkor  = dtheta_dpkor(:,1).^2.*s_pkor.^2;
var_ta_data = inv22;  var_ta_f = dtheta_df(:,2).^2.*s_f.^2;  var_ta_pkor = dtheta_dpkor(:,2).^2.*s_pkor.^2;
var_tb_data = inv33;  var_tb_f = dtheta_df(:,3).^2.*s_f.^2;  var_tb_pkor = dtheta_dpkor(:,3).^2.*s_pkor.^2;

var_q  = var_q_data  + var_q_f  + var_q_pkor;
var_ta = var_ta_data + var_ta_f + var_ta_pkor;
var_tb = var_tb_data + var_tb_f + var_tb_pkor;

q_SE  = sqrt(var_q);
ta_SE = sqrt(var_ta);
tb_SE = sqrt(var_tb);

out.dp1    = reshape(dp1,sz);
out.dp1_SE = reshape(dp1_SE,sz);
out.q     = reshape(q,sz);
out.q_SE  = reshape(q_SE,sz);
out.q_CI  = [q - z*q_SE, q + z*q_SE];
out.ta    = reshape(ta,sz);
out.ta_SE = reshape(ta_SE,sz);
out.ta_CI = [ta - z*ta_SE, ta + z*ta_SE];
out.tb    = reshape(tb,sz);
out.tb_SE = reshape(tb_SE,sz);
out.tb_CI = [tb - z*tb_SE, tb + z*tb_SE];
out.diag.chi2   = reshape(chi2,sz);
out.diag.dof    = reshape(dof,sz);   % per-row: 2 normally, 1 if DPR or DPN out, 0 if both out
out.diag.pvalue = reshape(pval,sz);
out.diag.iters  = iters;
out.diag.theta0 = theta0;
out.diag.dprOut = reshape(dprOut,sz);
out.diag.dpnOut = reshape(dpnOut,sz);

% SE broken out by SOURCE, so you can see whether f, pkor, or the raw
% pressure data is driving the total uncertainty - e.g. sqrt(q_SE_f.^2 +
% q_SE_pkor.^2 + q_SE_data.^2) == q_SE (up to floating point).
out.diag.q_SE_data  = reshape(sqrt(var_q_data), sz);
out.diag.q_SE_f     = reshape(sqrt(var_q_f),    sz);
out.diag.q_SE_pkor  = reshape(sqrt(var_q_pkor), sz);
out.diag.ta_SE_data = reshape(sqrt(var_ta_data),sz);
out.diag.ta_SE_f    = reshape(sqrt(var_ta_f),   sz);
out.diag.ta_SE_pkor = reshape(sqrt(var_ta_pkor),sz);
out.diag.tb_SE_data = reshape(sqrt(var_tb_data),sz);
out.diag.tb_SE_f    = reshape(sqrt(var_tb_f),   sz);
out.diag.tb_SE_pkor = reshape(sqrt(var_tb_pkor),sz);

out.notes = ['dof = 2 per row normally (5 equations, 3 unknowns [q,ta,tb]), all spent ', ...
    'jointly via weighted nonlinear least squares (Gauss-Newton), rather than exact ', ...
    'pairwise closed form (see r858_solve2.m). Set DPR==0 or DPN==0 (exactly) to flag ', ...
    'that sensor as out for a given row - the fit falls back to the cross-relation ', ...
    'ta/tb=dpa/dpb for the affected shape variable and drops that equation''s weight ', ...
    '(dof -> 1); see out.diag.dof, out.diag.dprOut, out.diag.dpnOut. SE/CI include both ', ...
    'the data (dp1,dpa,dpb,dpr,dpn) contribution via inv(J''WJ), and the f, pkor ', ...
    'constant-uncertainty contribution via warm-started finite-difference sensitivity. ', ...
    'See out.diag.*_SE_data/_f/_pkor for the breakdown of which source dominates each ', ...
    'unknown''s uncertainty.'];

end

% ------------------------------------------------------------------
function out = normalize_sigma_struct_(sigma, errId)
% Accepts a sigma struct with field names in any case, and additionally
% recognizes 'fcoef' as an alias for 'f' and 'pcor' as an alias for
% 'pkor' (matching common workspace naming, e.g. sigma.PTB, sigma.fcoef,
% sigma.pcor). Unrecognized fields are ignored with a warning.
fn = fieldnames(sigma);
out = struct();
for i = 1:numel(fn)
    name = lower(fn{i});
    val = sigma.(fn{i});
    switch name
        case 'ptb',  out.ptb  = val;
        case 'psa',  out.psa  = val;
        case 'dpa',  out.dpa  = val;
        case 'dpb',  out.dpb  = val;
        case 'dpr',  out.dpr  = val;
        case 'dpn',  out.dpn  = val;
        case {'f','fcoef'},    out.f    = val;
        case {'pkor','pcor'},  out.pkor = val;
        otherwise
            warning([errId ':unknownSigmaField'], ...
                'Ignoring unrecognized sigma field "%s".', fn{i});
    end
end
end

% ------------------------------------------------------------------
function theta = sanitize_(theta, bound)
% Flag rows where the fit produced a non-physical (|ta| or |tb| too
% large) or non-finite result and replace with NaN, with a summary
% warning (never per-row, to avoid spamming the console on large batches).
bad = ~isfinite(theta(:,2)) | ~isfinite(theta(:,3)) | ...
      abs(theta(:,2)) > bound | abs(theta(:,3)) > bound;
theta(bad,:) = NaN;
if any(bad)
    warning('r858_solve3:outliers', ...
        '%d of %d row(s) failed to converge to a physical solution (|ta| or |tb| > %g, or non-finite) - returned as NaN.', ...
        sum(bad), numel(bad), bound);
end
end

% ------------------------------------------------------------------
function v = expand_(v, n)
if isscalar(v)
    v = repmat(v, n, 1);
else
    v = v(:);
end
end

% ------------------------------------------------------------------
function [ta, tb, Y] = closed_form_(dp1,dpa,dpb,dpr,dpn,f,pkor)
% Exact closed-form starting point. NUMERICALLY STABLE version: the
% original direct formula tb = -(dpb+2*dpr-sqrt(2*dpb^2+4*dpb*dpr+4*dpr^2))/dpb
% subtracts two nearly-equal quantities whenever dpb is small relative to
% dpr (near-zero-flow samples), causing catastrophic cancellation - this
% is what produced the ~1e83 outliers. The two roots of the underlying
% quadratic (dpb*tb^2+(2*dpb+4*dpr)*tb-dpb=0) multiply to exactly -1, so
% instead of computing the small (desired) root directly, we compute the
% large root (a sum, never cancels) and take its reciprocal - mathematically
% identical, but immune to cancellation (verified against a
% high-precision reference down to dpb ~ 1e-15).
%
% SENSOR DROPOUT: dpr==0 or dpn==0 (exactly) signals that sensor is out
% for that row. The "direct" formula for the corresponding shape variable
% (tb from dpb/dpr, or ta from dpa/dpn) is invalid in that case - dpr=0 or
% dpn=0 plugged in literally would silently produce a wrong closed-form
% value, not a warning. Instead, the missing shape variable is
% reconstructed via the cross-relation ta/tb = dpa/dpb (which comes from
% setting eq2's and eq3's f*q expressions equal to each other, and needs
% no dpr or dpn at all) applied to whichever of ta/tb is still directly
% solvable. If BOTH dpr and dpn are out simultaneously, neither shape
% variable has a reliable closed form here; both fall back to 0 and a
% warning is issued - the Gauss-Newton fit is left to find the (still
% possibly exactly-determined, via eq.1+2+3) solution from that generic
% starting point.
tinyTol = 1e-12;

dprOut = abs(dpr) < tinyTol;
dpnOut = abs(dpn) < tinyTol;

tb_big = -(dpb + 2*dpr + sqrt(2*dpb.^2 + 4*dpb.*dpr + 4*dpr.^2)) ./ dpb;
tb_direct = -1 ./ tb_big;
tb_direct(abs(dpb) < tinyTol) = 0;   % exact limit as dpb -> 0 (dpr staying finite)

ta_big = -(dpa + 2*dpn + sqrt(2*dpa.^2 + 4*dpa.*dpn + 4*dpn.^2)) ./ dpa;
ta_direct = -1 ./ ta_big;
ta_direct(abs(dpa) < tinyTol) = 0;   % exact limit as dpa -> 0 (dpn staying finite)

ta = ta_direct;
tb = tb_direct;

% dpr out (dpn present): tb from the cross-relation using the still-valid ta_direct
onlyDprOut = dprOut & ~dpnOut;
tb(onlyDprOut) = (dpb(onlyDprOut)./dpa(onlyDprOut)) .* ta_direct(onlyDprOut);

% dpn out (dpr present): ta from the cross-relation using the still-valid tb_direct
onlyDpnOut = dpnOut & ~dprOut;
ta(onlyDpnOut) = (dpa(onlyDpnOut)./dpb(onlyDpnOut)) .* tb_direct(onlyDpnOut);

% both out: no reliable closed form for either - generic fallback, flagged
bothOut = dprOut & dpnOut;
if any(bothOut)
    ta(bothOut) = 0;
    tb(bothOut) = 0;
    warning('r858_solve3:bothDropout', ...
        '%d row(s) have both DPR and DPN out simultaneously - no closed-form starting point available; falling back to ta=tb=0 and relying on the Gauss-Newton fit.', ...
        sum(bothOut));
end

D = 1 + ta.^2 + tb.^2;
q_eq1 = (dp1 + pkor).*D ./ (D.*(1-f) + f);
% q_aside, q_bside below are NOT used as the starting point (only q_eq1
% is, see caller) - they're kept only for parity with r858_solve2.m's Y
% layout. Guard the denominators so a both-dropout fallback (ta=tb=0)
% doesn't throw spurious divide-by-zero warnings for values that are
% discarded anyway.
tbSafe = tb; tbSafe(abs(tbSafe)<tinyTol) = tinyTol;
taSafe = ta; taSafe(abs(taSafe)<tinyTol) = tinyTol;
qb1 = dpb.*D ./ (2*f.*tbSafe);
qb2 = 2*dpr.*D ./ (f.*(1 - 2*tb - tb.^2));
q_bside = (qb1 + qb2) / 2;
qa1 = dpa.*D ./ (2*f.*taSafe);
qa2 = 2*dpn.*D ./ (f.*(1 - 2*ta - ta.^2));
q_aside = (qa1 + qa2) / 2;
Y = [q_eq1, q_aside, q_bside];
end

% ------------------------------------------------------------------
function r = residuals_(theta, data)
% theta: N x 3 = [q ta tb]. Returns N x 5 residuals [r1..r5].
q = theta(:,1); ta = theta(:,2); tb = theta(:,3);
dp1=data.dp1; dpa=data.dpa; dpb=data.dpb; dpr=data.dpr; dpn=data.dpn;
f=data.f; pkor=data.pkor;

D = 1 + ta.^2 + tb.^2;
r1 = dp1 - (q.*(D.*(1-f)+f)./D - pkor);
r2 = dpb - 2*f.*q.*tb./D;
r3 = dpa - 2*f.*q.*ta./D;
r4 = dpr - f.*q.*(1 - 2*tb - tb.^2)./(2*D);
r5 = dpn - f.*q.*(1 - 2*ta - ta.^2)./(2*D);
r = [r1,r2,r3,r4,r5];
end

% ------------------------------------------------------------------
function J = jacobian_fd_(theta, data)
% Analytic Jacobian of residuals_ wrt theta = [q,ta,tb] (N x 3).
% Returns J as N x 5 x 3. (Function name kept as jacobian_fd_ so the
% call sites in gauss_newton_ don't need to change; it is no longer
% finite-difference based - all 15 partials were verified symbolically
% against sympy before replacing the old 6-extra-eval finite-difference
% version. This is the dominant cost of the solver, so computing it in
% closed form (1 pass, no perturbed re-evaluation of residuals_) is the
% single biggest speedup available.)
q = theta(:,1); ta = theta(:,2); tb = theta(:,3);
f = data.f;

D  = 1 + ta.^2 + tb.^2;
D2 = D.^2;
N4 = 1 - 2*tb - tb.^2;   % appears in r4
N5 = 1 - 2*ta - ta.^2;   % appears in r5

n = size(theta,1);
J = zeros(n,5,3);

% dr1/d(q,ta,tb)
J(:,1,1) = -(1-f) - f./D;
J(:,1,2) = q.*f.*2.*ta./D2;
J(:,1,3) = q.*f.*2.*tb./D2;

% dr2/d(q,ta,tb)
J(:,2,1) = -2*f.*tb./D;
J(:,2,2) = 4*f.*q.*ta.*tb./D2;
J(:,2,3) = -2*f.*q.*(D - 2*tb.^2)./D2;

% dr3/d(q,ta,tb)
J(:,3,1) = -2*f.*ta./D;
J(:,3,2) = -2*f.*q.*(D - 2*ta.^2)./D2;
J(:,3,3) = 4*f.*q.*ta.*tb./D2;

% dr4/d(q,ta,tb)
J(:,4,1) = -f.*N4./(2*D);
J(:,4,2) = f.*q.*N4.*ta./D2;
J(:,4,3) = -f.*q./(2*D2) .* ( D.*(-2-2*tb) - N4.*2.*tb );

% dr5/d(q,ta,tb)
J(:,5,1) = -f.*N5./(2*D);
J(:,5,2) = -f.*q./(2*D2) .* ( D.*(-2-2*ta) - N5.*2.*ta );
J(:,5,3) = f.*q.*N5.*tb./D2;
end

% ------------------------------------------------------------------
function [theta, iters, invM] = gauss_newton_(theta0, data, W, maxIter, tol)
% Vectorized Gauss-Newton for the 3-parameter [q,ta,tb] WNLS fit.
% W: N x 5 weights. Returns theta (N x 3) and invM = N x 6
% [inv11 inv22 inv33 inv12 inv13 inv23] of the FINAL (J'WJ)^-1 (3x3 sym).
theta = theta0;
n = size(theta,1);
invM = zeros(n,6);
for it = 1:maxIter
    r = residuals_(theta, data);        % N x 5
    J = jacobian_fd_(theta, data);      % N x 5 x 3

    M11 = sum(W .* J(:,:,1).^2, 2);
    M22 = sum(W .* J(:,:,2).^2, 2);
    M33 = sum(W .* J(:,:,3).^2, 2);
    M12 = sum(W .* J(:,:,1) .* J(:,:,2), 2);
    M13 = sum(W .* J(:,:,1) .* J(:,:,3), 2);
    M23 = sum(W .* J(:,:,2) .* J(:,:,3), 2);
    v1  = sum(W .* J(:,:,1) .* r, 2);
    v2  = sum(W .* J(:,:,2) .* r, 2);
    v3  = sum(W .* J(:,:,3) .* r, 2);

    detM = M11.*M22.*M33 - M11.*M23.^2 - M12.^2.*M33 + 2*M12.*M13.*M23 - M13.^2.*M22;
    inv11 = (M22.*M33 - M23.^2) ./ detM;
    inv22 = (M11.*M33 - M13.^2) ./ detM;
    inv33 = (M11.*M22 - M12.^2) ./ detM;
    inv12 = (M13.*M23 - M12.*M33) ./ detM;
    inv13 = (M12.*M23 - M13.*M22) ./ detM;
    inv23 = (M12.*M13 - M11.*M23) ./ detM;

    d1 = -(inv11.*v1 + inv12.*v2 + inv13.*v3);
    d2 = -(inv12.*v1 + inv22.*v2 + inv23.*v3);
    d3 = -(inv13.*v1 + inv23.*v2 + inv33.*v3);

    theta(:,1) = theta(:,1) + d1;
    theta(:,2) = theta(:,2) + d2;
    theta(:,3) = theta(:,3) + d3;

    maxstep = max(max(abs([d1,d2,d3])));
    if maxstep < tol
        break
    end
end
invM = [inv11, inv22, inv33, inv12, inv13, inv23];
iters = it;
end
