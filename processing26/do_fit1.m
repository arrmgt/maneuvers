function [x,fu,fv,fw,mag,dir,resid,jacobian,CI,beta_samp]=do_fit1(method,fun,X)

switch method
    case 'lsqnonlin'
        options=optimset('lsqnonlin');
    case 'fsolve'
        options=optimset('fsolve');
end

%fitting function
switch method
    case 'lsqnonlin'
        %Upper and lower limits of parameters
        LB(1) = -1*pi./180;   UB(1) = +1*pi./180;   % roll offset
        LB(2) = -1*pi./180;   UB(2) = +1*pi./180;   % pitch offset
        LB(3) = -2*pi./180;   UB(3) = +2*pi./180;   % heading offset
        LB(4) = 0.0;          UB(4) = +1;           % attack factor
        LB(5) = -0.2;         UB(5) = +0.2;         % attack offset
        LB(6) = 0.0;          UB(6) = +1;           % sideslip factor
        LB(7) = -0.2;         UB(7) = +0.2;         % sideslip offset
        LB(8) = -20;          UB(8) = +20;          % ps offset
        options=optimset('lsqnonlin');
        options=optimset(options,'Display','iter');
        %options=optimset(options,'MaxFunEvals',5000);
        %options=optimset(options,'TolX',5.e-3,'TolFun',5.e-3);
        x0=zeros(1,8);
        try
            [x,resnorm,resid,exitflag,output,lambda,jacobian]=lsqnonlin(fun,x0,LB,UB,options,X);
        catch ME
            catchME(ME)
        end
    case 'fsolve'
        options=optimset('fsolve');
        options=optimset(options,'Display','iter');
        %options=optimset(options,'MaxFunEvals',5000);
        %options=optimset(options,'TolX',5.e-3,'TolFun',5.e-3);
        [x,fval,exitflag,output,jacobian]=fsolve(fun,x0,options,X);
        resid=[fu,fv,fw];
end
%
J = jacobian;                  % 105033×8 sparse
r = resid(:);                  % stack residuals to match J rows
m = size(J,1);
p = numel(x);

mse = sum(r.^2) / (m - p);

[Q,R] = qr(full(J),0);         % R is p×p (8×8)
if rank(R) == p
    CovB = mse * inv(R'*R);
else
    CovB = mse * pinv(R'*R);   % fallback for rank-deficient case
end
CI = nlparci(x, r, 'Covar', CovB);

[f,fu,fv,fw,va,vg,mmm]=wcal8(x,X);

mag=sqrt(fu.^2+fv.^2);
dir=(unwrap(atan2(fu,fv),pi)).*180./pi;

% Do some statistics on result

% Inputs (from lsqnonlin)
beta = x(:);        % 8x1
J = jacobian;              % N x 8 (stacked scalar outputs)
r = resid(:);         % N x 1 (same stacking as J)
N = numel(r);
p = numel(beta);
df = N - p;

% Covariance
Jfull = full(J);                 % Jfull is 17070x8 dense
[U,S,V] = svd(Jfull, 'econ');    % works now
s = diag(S);
tol = max(size(Jfull)) * eps(max(s));
s_inv = zeros(size(s));
s_inv(s > tol) = 1 ./ s(s > tol);
mse = sum(r(:).^2) / (numel(r) - numel(beta));
CovB = mse * (V * diag(s_inv.^2) * V');   % p x p


% Draw samples (multivariate normal)
nSamp = 10000;
beta_samp = mvnrnd(beta.', CovB, nSamp);   % nSamp x 8

% Summaries
beta_mean = mean(beta_samp, 1).';
beta_med  = median(beta_samp, 1).';
beta_ci   = prctile(beta_samp, [2.5, 97.5])';  % 8x2 [lower upper]

end

