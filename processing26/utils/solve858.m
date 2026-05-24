function [q, f, ta_out, tb_out, sigma_q, sigma_f, res, stats] = ...
         solve858(dp1a, dpa, dpb, varargin)
%SOLVE858  Vectorised closed-form solve for Rosemount 858 multi-hole probe.
%
%  Three required differential-pressure inputs, four unknowns (q,f,ta,tb).
%  At least one of 'dpr' or 'DPN' must be supplied to determine the angles.
%  Each measurement beyond the minimum adds one degree of freedom.
%
% -----------------------------------------------------------------------
%  EQUATION STRUCTURE  (P=1+ta²+tb², S=ta²+tb², fq≡f·q)
% -----------------------------------------------------------------------
%  [req]  eq1:  P·q − S·fq  = dp1a·P        centre-hole vs static
%  [req]  eq2:  2·ta·fq      = dpa·P          alpha pair
%  [req]  eq3:  2·tb·fq      = dpb·P          beta pair
%  [opt]  eq4:  cb·fq        = 2·dpr·P        p1 − right sideslip (cb=1−2tb−tb²)
%  [opt]  eq5:  ca·fq        = 2·DPN·P        p1 − bottom attack  (ca=1−2ta−ta²)
%  [opt]  eq6:  P·q − S·fq  = ptb_diff·P    ptb_abs − pstatic
%
%  ANGLE DETERMINATION STRATEGY
%  -----------------------------
%    dpr only  :  tb from dpb+dpr quadratic;  ta = tb·dpa/dpb  (ratio)
%    DPN only  :  ta from dpa+DPN quadratic;  tb = ta·dpb/dpa  (ratio)
%    dpr + DPN :  ta from dpa+DPN quadratic;  tb from dpb+dpr quadratic
%                 (best conditioned — avoids ratio instability at small dpb or dpa)
%
%  DEGREES OF FREEDOM
%  ------------------
%    dof = (use_dpr + use_DPN − 1) + use_ptb + use_fsim
%
%               | dpr only | DPN only | dpr + DPN |
%  no ptb,fsim  |    0     |    0     |     1     |
%  + ptb        |    1     |    1     |     2     |
%  + f_sim      |    1     |    1     |     2     |
%  + ptb+f_sim  |    2     |    2     |     3     |
%
%  dof=0: exact solution; sigma_resid=NaN; uncertainty from sigma_trans only.
%
%  UNCERTAINTY SOURCES
%  -------------------
%    sigma_resid  from LS residuals  (model / consistency error)
%    sigma_trans  propagated from transducer noise via central FD Jacobian
%    sigma_total  = sqrt(sigma_resid² + sigma_trans²)
%
% -----------------------------------------------------------------------
%  SYNTAX
%    % dpr only  (q_beta, dof=0)
%    [q,f,ta,tb] = solve858(dp1a, dpa, dpb, 'dpr', DPR)
%
%    % DPN only  (q_alpha, dof=0)
%    [q,f,ta,tb] = solve858(dp1a, dpa, dpb, 'DPN', DPN)
%
%    % both dpr+DPN  (dof=1)
%    [q,f,ta,tb,sq,sf,res,stats] = solve858(dp1a, dpa, dpb, ...
%                    'dpr',          DPR_vec,   ...
%                    'DPN',          DPN_vec,   ...
%                    'f_sim',        f_vec,     ...
%                    'ptb_abs',      PTB_vec,   'pstatic',       PS_vec,  ...
%                    'sigma_dp',     0.025,           ...
%                    'sigma_ptb',    1.0,       'sigma_pstatic', 1.0,     ...
%                    'ps_cor',       ps_cor_vec,                 ...
%                    'alpha',        0.05)
%
%  REQUIRED INPUTS  (all N×1)
%    dp1a   centre-hole minus static  differential pressure
%    dpa    alpha-pair  differential pressure
%    dpb    beta-pair   differential pressure
%
%  OPTIONAL NAME-VALUE INPUTS
%    'dpr'           N×1  p1 − right sideslip port  (at least one of dpr/DPN needed)
%    'DPN'           N×1  p1 − bottom attack port   (at least one of dpr/DPN needed)
%    'ptb_abs'       N×1  absolute pressure at port 1
%    'pstatic'       N×1  static pressure — both ptb_abs and pstatic required together
%    'sigma_ptb'     1-sigma for ptb_abs transducer            (default 0)
%    'sigma_pstatic' 1-sigma for pstatic transducer            (default 0)
%                    Combined: sigma_ptb_diff = sqrt(sigma_ptb²+sigma_pstatic²)
%    'f_sim'         N×1  predicted f — q solved by full LS, adds 1 dof
%    'sigma_dp'      1-sigma for [dp1a;dpa;dpb;(dpr);(DPN)]
%                    scalar = same for all; vector length 3, 4, or 5
%    'ps_cor'        N×1 or scalar  static pressure defect correction
%                    Convention: ps_cor = ps_measured − ps_true
%                    Applied:  dp1a_true     = dp1a     + ps_cor
%                              ptb_diff_true = ptb_diff + ps_cor  (if ptb supplied)
%                    Unaffected: dpr, DPN, dpa, dpb (probe-port references)
%    'alpha'         CI significance level                  (default 0.05)
%    'q_bounds'      [qmin qmax]                            (default [10 150])
%    'f_bounds'      [fmin fmax]                            (default [1 3])
%    'ta_bounds'     [tamin tamax]                          (default [-2 2])
%    'tb_bounds'     [tbmin tbmax]                          (default [-2 2])
%
%  OUTPUTS  (N×1 unless noted)
%    q        impact pressure         (0 where bounds violated)
%    f        probe sensitivity       (= f_sim when supplied; 0 if invalid)
%    ta_out   tan(angle of attack)    (NaN where invalid)
%    tb_out   tan(sideslip angle)     (NaN where invalid)
%    sigma_q  total 1-sigma on q      (NaN where invalid or dof=0 with no sigma_trans)
%    sigma_f  total 1-sigma on f      (NaN when f_sim supplied, dof=0, or invalid)
%    res      residuals N×n_eq
%             Columns (in order, when present):
%               [r_ptb, r_dp1a, r_dpa, r_dpb, r_dpr, r_DPN]
%               r_ptb / r_dp1a : only when ptb supplied (no f_sim) or f_sim supplied
%               r_dp1a alone   : only when f_sim (no ptb) — dp1a in LS
%               r_dpr          : only when dpr supplied
%               r_DPN          : only when DPN supplied
%    stats    struct (see STATS FIELDS)
%
%  STATS FIELDS
%    Per-point: .q_ci,.f_ci (N×2), .R2 (N×1), .valid (N×1 logical)
%               .sigma_q_resid, .sigma_q_trans
%               .sigma_f_resid, .sigma_f_trans
%               .sigma_ta, .sigma_tb, .sigma_alpha, .sigma_beta
%    Global:    .N_valid, .q_mean,.q_std,.q_ci_global
%               .f_mean, .f_std, .f_ci_global, .R2_median
%    Info:      .dof, .use_dpr, .use_DPN, .use_ptb, .use_fsim, .alpha, .sigma_dp
% -----------------------------------------------------------------------

% ---- parse inputs -------------------------------------------------------
p = inputParser;
addRequired(p,  'dp1a');
addRequired(p,  'dpa');
addRequired(p,  'dpb');
addParameter(p, 'dpr',           []);
addParameter(p, 'DPN',           []);
addParameter(p, 'ptb_abs',       []);
addParameter(p, 'pstatic',       []);
addParameter(p, 'sigma_ptb',     0);
addParameter(p, 'sigma_pstatic', 0);
addParameter(p, 'f_sim',         []);
addParameter(p, 'sigma_dp',      0);
addParameter(p, 'ps_cor',        0);
addParameter(p, 'q_bounds',      [5,  150]);
addParameter(p, 'f_bounds',      [1,   3  ]);
addParameter(p, 'ta_bounds',     [-2,  2  ]);
addParameter(p, 'tb_bounds',     [-2,  2  ]);
addParameter(p, 'alpha',         0.05);
parse(p, dp1a, dpa, dpb, varargin{:});

dpr      = p.Results.dpr;
DPN      = p.Results.DPN;
ptb_abs  = p.Results.ptb_abs;
pstatic  = p.Results.pstatic;
sig_ptb  = p.Results.sigma_ptb;
sig_ps   = p.Results.sigma_pstatic;
f_sim    = p.Results.f_sim;
sigma_dp = p.Results.sigma_dp(:);
ps_cor   = p.Results.ps_cor(:);
q_bounds = p.Results.q_bounds;
f_bounds = p.Results.f_bounds;
ta_bnd   = p.Results.ta_bounds;
tb_bnd   = p.Results.tb_bounds;
alpha    = p.Results.alpha;

use_dpr  = ~isempty(dpr);
use_DPN  = ~isempty(DPN);
use_ptb  = ~isempty(ptb_abs) && ~isempty(pstatic);
use_fsim = ~isempty(f_sim);

if ~use_dpr && ~use_DPN
    error('solve858:noAngleSensor', ...
        'Supply at least one of ''dpr'' or ''DPN''.');
end

% ---- ps_cor size check --------------------------------------------------
if ~isscalar(ps_cor) && numel(ps_cor) ~= numel(dp1a)
    error('solve858:pscorDim', ...
        'ps_cor must be scalar or length %d.', numel(dp1a));
end

% ---- ptb_diff -----------------------------------------------------------
if use_ptb
    ptb_diff       = ptb_abs - pstatic;
    sigma_ptb_diff = sqrt(sig_ptb.^2 + sig_ps.^2);
else
    ptb_diff       = zeros(size(dp1a));
    sigma_ptb_diff = 0;
end

% ---- static pressure defect correction ----------------------------------
%  ps_cor = ps_measured − ps_true.  Only dp1a and ptb_diff are referenced
%  to aircraft static; dpr, DPN, dpa, dpb are probe-port differentials.
if any(ps_cor ~= 0)
    dp1a = dp1a + ps_cor;
    if use_ptb
        ptb_diff = ptb_diff + ps_cor;
    end
end

% ---- degrees of freedom -------------------------------------------------
dof_val = (use_dpr + use_DPN - 1) + use_ptb + use_fsim;   % 0, 1, 2, or 3

% n_ctr: number of [P,−S] rows  (1 without ptb, 2 with ptb)
n_ctr = 1 + double(use_ptb);

% ---- sigma_dp -----------------------------------------------------------
%  Order: [dp1a; dpa; dpb; (dpr if present); (DPN if present)]
n_dp_main = 3 + use_dpr + use_DPN;
if isscalar(sigma_dp)
    sigma_dp = sigma_dp * ones(n_dp_main, 1);
elseif numel(sigma_dp) ~= n_dp_main
    lbl = 'dp1a;dpa;dpb';
    if use_dpr, lbl = [lbl ';dpr']; end
    if use_DPN, lbl = [lbl ';DPN']; end
    error('solve858:sigmaDim', ...
        'sigma_dp must be scalar or length %d [%s].', n_dp_main, lbl);
end

% ---- core solve ---------------------------------------------------------
[q_raw, f_raw, ta_raw, tb_raw, fq, S, P, ca, cb, D] = ...
    solve_core(ptb_diff, dp1a, dpa, dpb, dpr, DPN, ...
               use_dpr, use_DPN, use_ptb, n_ctr, ta_bnd, tb_bnd, f_sim);

% ---- bounds & validity --------------------------------------------------
if use_fsim
    valid = q_raw >= q_bounds(1) & q_raw <= q_bounds(2);
else
    valid = q_raw >= q_bounds(1) & q_raw <= q_bounds(2) & ...
            f_raw >= f_bounds(1) & f_raw <= f_bounds(2);
end

q = q_raw;  f = f_raw;
q(~valid) = 0;  f(~valid) = 0;
ta_out = ta_raw;  tb_out = tb_raw;
ta_out(~valid) = NaN;  tb_out(~valid) = NaN;

if nargout < 5, return; end

% ---- residuals ----------------------------------------------------------
%  With f_sim:  all equations contribute to the q LS; all have residuals.
%  Without f_sim:
%    dp1a (and ptb_diff when ptb supplied) give q exactly / averaged.
%    dpa, dpb, dpr, DPN give the fq LS residuals.
%    When ptb present: r_ptb = P*(dp1a−ptb_diff)/2 = −r_dp1a.
%    When ptb absent and no f_sim: dp1a is an exact equation — no residual.

if use_fsim
    r_dp1a = (P - f_raw.*S).*q_raw - dp1a.*P;
    r_dpa  = 2.*ta_raw.*f_raw.*q_raw - dpa.*P;
    r_dpb  = 2.*tb_raw.*f_raw.*q_raw - dpb.*P;
    cols   = {r_dp1a, r_dpa, r_dpb};
    if use_ptb
        r_ptb = (P - f_raw.*S).*q_raw - ptb_diff.*P;
        cols  = [{r_ptb}, cols];
    end
    if use_dpr, cols{end+1} = cb.*f_raw.*q_raw - 2.*dpr.*P;  end
    if use_DPN, cols{end+1} = ca.*f_raw.*q_raw - 2.*DPN.*P;  end
else
    r_dpa = 2.*ta_raw.*fq - dpa.*P;
    r_dpb = 2.*tb_raw.*fq - dpb.*P;
    if use_ptb
        r_ptb  =  P .* (dp1a - ptb_diff) ./ 2;
        cols   = {r_ptb, -r_ptb, r_dpa, r_dpb};
    else
        cols   = {r_dpa, r_dpb};
    end
    if use_dpr, cols{end+1} = cb.*fq - 2.*dpr.*P;  end
    if use_DPN, cols{end+1} = ca.*fq - 2.*DPN.*P;  end
end
res = [cols{:}];
res(~valid,:) = 0;

% ---- residual-based sigma -----------------------------------------------
if dof_val > 0
    sigma2_resid = sum(res.^2, 2) ./ dof_val;
else
    sigma2_resid = NaN(size(q_raw));
end

if use_fsim
    % Var(q) = sigma² / C'C  where C'C = n_ctr*(P−fS)² + f²*D_fq
    % (D returned from solve_core is C'C in the f_sim branch)
    var_q_resid   = sigma2_resid ./ D;
    sigma_f_resid = NaN(size(f_raw));
else
    % (A'A) normal equations with n_ctr rows of [P,−S]:
    %   [1,1] = (n_ctr·S²+D) / (n_ctr·P²·D)
    %   [2,2] = 1/D   [1,2] = S/(P·D)   (independent of n_ctr)
    var_q_resid = sigma2_resid .* (n_ctr.*S.^2 + D) ./ (n_ctr.*P.^2 .* D);
    var_f_resid = sigma2_resid ./ (q_raw.^2 .* n_ctr.*P.^2 .* D) .* ...
                  (f_raw.^2.*(n_ctr.*S.^2+D) - 2.*f_raw.*n_ctr.*S.*P + n_ctr.*P.^2);
    sigma_f_resid = sqrt(var_f_resid);
    sigma_f_resid(~valid) = NaN;
end
sigma_q_resid = sqrt(var_q_resid);
sigma_q_resid(~valid) = NaN;
if use_fsim, sigma_f_resid(:) = NaN; end

% ---- transducer uncertainty via central FD Jacobian ---------------------
%  Perturbed inputs: dp1a, dpa, dpb, [dpr], [DPN], [ptb_diff]
%  sigma_dp order:   1,    2,   3,   [4],   [4|5], combined for ptb_diff
wide_bnd = [-10, 10];

fd_names = {'dp1a','dpa','dpb'};
fd_vals  = {dp1a, dpa, dpb};
fd_sigma = sigma_dp(1:3);
if use_dpr
    fd_names{end+1} = 'dpr';
    fd_vals{end+1}  = dpr;
    fd_sigma(end+1) = sigma_dp(4);
end
if use_DPN
    fd_names{end+1} = 'DPN';
    fd_vals{end+1}  = DPN;
    fd_sigma(end+1) = sigma_dp(3 + use_dpr + 1);
end
if use_ptb
    fd_names{end+1} = 'ptb_diff';
    fd_vals{end+1}  = ptb_diff;
    fd_sigma(end+1) = sigma_ptb_diff;
end
n_fd = numel(fd_vals);

sq_q  = zeros(size(q_raw));
sq_f  = zeros(size(f_raw));
sq_ta = zeros(size(ta_raw));
sq_tb = zeros(size(tb_raw));

if any(fd_sigma > 0)
    h = 1e-4;

    for k = 1:n_fd
        if fd_sigma(k) == 0, continue; end

        fd_p = fd_vals;  fd_p{k} = fd_vals{k} + h;
        fd_m = fd_vals;  fd_m{k} = fd_vals{k} - h;

        % reconstruct named inputs
        dp1a_p = fd_p{1}; dpa_p = fd_p{2}; dpb_p = fd_p{3};
        dp1a_m = fd_m{1}; dpa_m = fd_m{2}; dpb_m = fd_m{3};
        dpr_p  = dpr;      dpr_m  = dpr;
        DPN_p  = DPN;      DPN_m  = DPN;
        ptd_p  = ptb_diff; ptd_m  = ptb_diff;

        for j = 4:n_fd
            switch fd_names{j}
                case 'dpr',      dpr_p = fd_p{j}; dpr_m = fd_m{j};
                case 'DPN',      DPN_p = fd_p{j}; DPN_m = fd_m{j};
                case 'ptb_diff', ptd_p = fd_p{j}; ptd_m = fd_m{j};
            end
        end

        [qp,fp,tap,tbp] = solve_core(ptd_p,dp1a_p,dpa_p,dpb_p,dpr_p,DPN_p,...
                                      use_dpr,use_DPN,use_ptb,n_ctr,...
                                      wide_bnd,wide_bnd,f_sim);
        [qm,fm,tam,tbm] = solve_core(ptd_m,dp1a_m,dpa_m,dpb_m,dpr_m,DPN_m,...
                                      use_dpr,use_DPN,use_ptb,n_ctr,...
                                      wide_bnd,wide_bnd,f_sim);

        sc    = (fd_sigma(k) / (2*h))^2;
        sq_q  = sq_q  + sc .* (qp - qm).^2;
        sq_f  = sq_f  + sc .* (fp - fm).^2;
        sq_ta = sq_ta + sc .* (tap - tam).^2;
        sq_tb = sq_tb + sc .* (tbp - tbm).^2;
    end
end

sigma_q_trans  = sqrt(sq_q);
sigma_f_trans  = sqrt(sq_f);
sigma_ta_trans = sqrt(sq_ta);
sigma_tb_trans = sqrt(sq_tb);

sigma_q_trans(~valid)  = NaN;
sigma_f_trans(~valid)  = NaN;
sigma_ta_trans(~valid) = NaN;
sigma_tb_trans(~valid) = NaN;
if use_fsim
    sigma_f_trans(:) = NaN;
    sigma_f_resid(:) = NaN;
end

% ---- angle uncertainties (delta method on FD Jacobian) ------------------
deg         = 180 / pi;
sigma_alpha = sigma_ta_trans .* deg ./ (1 + ta_raw.^2);
sigma_beta  = sigma_tb_trans .* deg ./ (1 + tb_raw.^2);
sigma_alpha(~valid) = NaN;
sigma_beta(~valid)  = NaN;

% ---- combined total uncertainty -----------------------------------------
sigma_q = sqrt(sigma_q_resid.^2 + sigma_q_trans.^2);
sigma_f = sqrt(sigma_f_resid.^2 + sigma_f_trans.^2);

if nargout < 8, return; end

% ---- confidence intervals -----------------------------------------------
if dof_val > 0
    try,   t_pt = tinv(1-alpha/2, dof_val);
    catch, t_pt = norminv(1-alpha/2); end
else
    t_pt = NaN;
end
q_ci = [q - t_pt.*sigma_q,  q + t_pt.*sigma_q];
f_ci = [f - t_pt.*sigma_f,  f + t_pt.*sigma_f];
q_ci(~valid,:) = NaN;
f_ci(~valid,:) = NaN;

% ---- R² (over equations represented in res) -----------------------------
SS_tot = (dpa.*P).^2 + (dpb.*P).^2;
if use_dpr,   SS_tot = SS_tot + (2.*dpr.*P).^2;              end
if use_DPN,   SS_tot = SS_tot + (2.*DPN.*P).^2;              end
if use_ptb
    SS_tot = SS_tot + (ptb_diff.*P).^2 + (dp1a.*P).^2;
elseif use_fsim
    SS_tot = SS_tot + (dp1a.*P).^2;
end
R2 = 1 - sum(res.^2, 2) ./ SS_tot;
R2(~valid) = NaN;

% ---- global statistics --------------------------------------------------
N_valid = sum(valid);
q_v    = q_raw(valid);  f_v = f_raw(valid);
q_mean = mean(q_v);     q_std = std(q_v);
f_mean = mean(f_v);     f_std = std(f_v);
q_sem  = q_std / sqrt(N_valid);
f_sem  = f_std / sqrt(N_valid);
try,   t_gl = tinv(1-alpha/2, N_valid-1);
catch, t_gl = norminv(1-alpha/2); end

% ---- pack stats ---------------------------------------------------------
stats.q_ci          = q_ci;
stats.f_ci          = f_ci;
stats.R2            = R2;
stats.valid         = valid;
stats.N_valid       = N_valid;
stats.q_mean        = q_mean;
stats.q_std         = q_std;
stats.q_ci_global   = [q_mean - t_gl*q_sem,  q_mean + t_gl*q_sem];
stats.f_mean        = f_mean;
stats.f_std         = f_std;
stats.f_ci_global   = [f_mean - t_gl*f_sem,  f_mean + t_gl*f_sem];
stats.R2_median     = median(R2(valid));
stats.sigma_q_resid = sigma_q_resid;
stats.sigma_q_trans = sigma_q_trans;
stats.sigma_f_resid = sigma_f_resid;
stats.sigma_f_trans = sigma_f_trans;
stats.sigma_ta      = sigma_ta_trans;
stats.sigma_tb      = sigma_tb_trans;
stats.sigma_alpha   = sigma_alpha;
stats.sigma_beta    = sigma_beta;
stats.dof           = dof_val;
stats.use_dpr       = use_dpr;
stats.use_DPN       = use_DPN;
stats.use_ptb       = use_ptb;
stats.use_fsim      = use_fsim;
stats.alpha         = alpha;
stats.sigma_dp      = sigma_dp;
stats.ps_cor        = ps_cor;

end   % main function

% =========================================================================
function [q, f, ta, tb, fq, S, P, ca, cb, D] = ...
         solve_core(ptb_diff, dp1a, dpa, dpb, dpr, DPN, ...
                    use_dpr, use_DPN, use_ptb, n_ctr, ta_bnd, tb_bnd, f_sim)
%SOLVE_CORE  Bare computation — called by solve858 and its FD Jacobian.
%
%  Angle determination:
%    DPN present : ta from dpa+DPN quadratic (eq2÷eq5)
%    dpr present : tb from dpb+dpr quadratic (eq3÷eq4)
%    dpr only    : ta from ratio  ta = tb·dpa/dpb
%    DPN only    : tb from ratio  tb = ta·dpb/dpa

use_fsim_here = nargin >= 13 && ~isempty(f_sim);

ta = zeros(size(dpa));
tb = zeros(size(dpb));

% ---- ta: from DPN quadratic or filled later via ratio -------------------
if use_DPN
    % dpa·ta² + (2·dpa + 4·DPN)·ta − dpa = 0
    b_a    = 2.*dpa + 4.*DPN;
    disc_a = b_a.^2 + 4.*dpa.^2;
    ta_r1  = (-b_a + sqrt(disc_a)) ./ (2.*dpa);
    ta_r2  = (-b_a - sqrt(disc_a)) ./ (2.*dpa);
    ta(dpa > 0) = max(ta_r1(dpa > 0), ta_r2(dpa > 0));
    ta(dpa < 0) = min(ta_r1(dpa < 0), ta_r2(dpa < 0));
end

% ---- tb: from dpr quadratic or filled later via ratio -------------------
if use_dpr
    % dpb·tb² + (2·dpb + 4·dpr)·tb − dpb = 0
    b_b    = 2.*dpb + 4.*dpr;
    disc_b = b_b.^2 + 4.*dpb.^2;
    tb_r1  = (-b_b + sqrt(disc_b)) ./ (2.*dpb);
    tb_r2  = (-b_b - sqrt(disc_b)) ./ (2.*dpb);
    tb(dpb > 0) = max(tb_r1(dpb > 0), tb_r2(dpb > 0));
    tb(dpb < 0) = min(tb_r1(dpb < 0), tb_r2(dpb < 0));
end

% ---- ratio fill-in for the angle not determined by quadratic ------------
if use_DPN && ~use_dpr
    % tb from ratio  tb = ta·dpb/dpa  (DPN-only mode)
    nz = abs(dpa) > eps * max(abs(dpa(:)));
    tb(nz) = ta(nz) .* dpb(nz) ./ dpa(nz);
elseif use_dpr && ~use_DPN
    % ta from ratio  ta = tb·dpa/dpb  (dpr-only mode)
    nz = abs(dpb) > eps * max(abs(dpb(:)));
    ta(nz) = tb(nz) .* dpa(nz) ./ dpb(nz);
end
% Both present: ta from DPN quad, tb from dpr quad — no ratio needed.

ta = max(ta_bnd(1), min(ta_bnd(2), ta));
tb = max(tb_bnd(1), min(tb_bnd(2), tb));

% ---- geometry -----------------------------------------------------------
S  = ta.^2 + tb.^2;
P  = 1 + S;
ca = 1 - 2.*ta - ta.^2;
cb = 1 - 2.*tb - tb.^2;

% ---- D and Nv: include only available angle-sensor equations ------------
D  = 4.*S;
Nv = 2.*tb.*dpb + 2.*ta.*dpa;
if use_dpr
    D  = D  + cb.^2;
    Nv = Nv + 2.*cb.*dpr;
end
if use_DPN
    D  = D  + ca.^2;
    Nv = Nv + 2.*ca.*DPN;
end

% ---- dp_sum: n_ctr centre-hole measurements -----------------------------
if use_ptb
    dp_sum = ptb_diff + dp1a;   % n_ctr = 2
else
    dp_sum = dp1a;               % n_ctr = 1
end

% ---- Stage 3: q and f ---------------------------------------------------
if use_fsim_here
    % f given — scalar LS for q over all equations
    %   C'C = n_ctr·(P−f·S)² + f²·D_fq
    D_fq  = 4.*S;
    if use_dpr, D_fq = D_fq + cb.^2; end
    if use_DPN, D_fq = D_fq + ca.^2; end
    c_ctr = P - f_sim.*S;
    num   = c_ctr .* dp_sum .* P  +  f_sim .* Nv .* P;
    den   = n_ctr .* c_ctr.^2    +  f_sim.^2 .* D_fq;
    q     = num ./ den;
    f     = f_sim;
    fq    = f .* q;
    D     = den;   % C'C — used as denominator in variance formula
else
    % f unknown — fq from LS; q from dp_ctr + S·Nv/D
    fq = P .* Nv ./ D;
    q  = dp_sum ./ n_ctr  +  S .* Nv ./ D;
    f  = fq ./ q;
end

end   % solve_core