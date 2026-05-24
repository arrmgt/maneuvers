%%%%%%%%%%%%%%%%%%%%%%% Helper
% Outliers are points are outside indices (kk)
% Leading outliers are forced to first good point
% Trailing outliers are force to last good point
function x = setOutPoints(kk,x)
v1 = kk(1)-1; v2 = kk(end)+1;
x(1:v1) = x(1); x(v2:end) = x(end);
