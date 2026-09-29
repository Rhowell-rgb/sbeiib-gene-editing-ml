function Z = clr_transform(X, floorVal)
%CLR_TRANSFORM  Centred log-ratio transform for compositional data (rows sum to ~100).
%   Zeros/negatives are replaced by floorVal (default 0.01) before taking logs.
if nargin < 2, floorVal = 0.01; end
X(X <= 0) = floorVal;
L = log(X);
Z = L - mean(L, 2);             % rows with any NaN stay NaN
end
