function [yhat, lm, BETA, T, ci] = calc_lm (X,variableNames)

[lm,beta,yhat,R2, T]=do_regress(X,variableNames);
R2 = lm.Rsquared.Ordinary;
BETA = lm.Coefficients.Estimate;   % [Intercept; ...]
yhat = lm.Fitted;
ci = coefCI(lm);        % Nx2 matrix, lower and upper bounds
%C = [1 1 1 1];
%d = 0
%p = coefTest(lm, C, d); % general linear hypothesis C*beta = d
%ANOVA=anova(lm,'summary')

[yhat_ship, lm_ship, beta_ship, Tship] =do_lm(X,variableNames);
