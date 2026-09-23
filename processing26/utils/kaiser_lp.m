function b = kaiser_lp(fs, fc, transWidth, Ap, Ast)
% KAISER_LP  Kaiser-windowed FIR lowpass filter.
 
if Ast >= 50
    beta = 0.1102 * (Ast - 8.7);
elseif Ast >= 21
    beta = 0.5842 * (Ast - 21)^0.4 + 0.07886 * (Ast - 21);
else
    beta = 0;
end
 
deltaP = (10^(Ap/20) - 1) / (10^(Ap/20) + 1);
deltaS = 10^(-Ast/20);
delta  = min(deltaP, deltaS);
Aeff   = -20 * log10(delta);
 
deltaW = 2 * pi * transWidth / fs;
N      = ceil((Aeff - 8) / (2.285 * deltaW));
if mod(N, 2) ~= 0, N = N + 1; end
 
fc_mid = fc + transWidth / 2;    % sinc cutoff at midpoint of transition band
fc_n   = fc_mid / (fs / 2);
 
n       = 0 : N;
h_ideal = fc_n * sinc(fc_n * (n - N/2));
w       = kaiser(N+1, beta)';
b       = h_ideal .* w;
b       = b / sum(b);
 
end