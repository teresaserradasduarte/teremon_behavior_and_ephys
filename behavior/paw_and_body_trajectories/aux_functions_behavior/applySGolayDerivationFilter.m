function [ y ] = applySGolayDerivationFilter( x, N, K, F)
%SGOLAY_DERIVATION_FILTER returns the smooth Nth-derivative of signal x
%
% Inputs:
%   x = row vector signal
%   N = derivation order
%   K = order of polinomial fit
%   F = window length

% make sure x is a row vector
[r, ~] = size(x);
if(r ~= 1)
    x = x';
end

if(N < 2)
    factor = 1;
elseif(N == 2)
    factor = 2;
else
    error('Derivatives above the 2nd order are not handled. Check the factor required');
end

[b,g]=sgolay(K,F);

for n = (F+1)/2: length(x)-(F+1)/2,
    y(n) = factor * g(:,N+1)' * x(n - (F+1)/2 + 1: n + (F+1)/2 - 1)';
end

y(length(x)) = 0;

end