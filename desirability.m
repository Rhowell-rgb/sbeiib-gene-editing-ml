function d = desirability(y, type, lims, s)
%DESIRABILITY  Derringer-Suich desirability in [0,1].
%   'min'   : lims = [L U]  -> 1 if y<=L, 0 if y>=U        (e.g. GI)
%   'max'   : lims = [L U]  -> 0 if y<=L, 1 if y>=U        (e.g. RS)
%   'target': lims = tol, y = relative change vs control -> 1 at 0, 0 at |y|>=tol
if nargin < 4, s = 1; end
switch type
    case 'min'
        d = ((lims(2) - y) ./ (lims(2) - lims(1))) .^ s;
    case 'max'
        d = ((y - lims(1)) ./ (lims(2) - lims(1))) .^ s;
    case 'target'
        d = (1 - abs(y) ./ lims) .^ s;
end
d = min(max(real(d), 0), 1);
d(isnan(y)) = NaN;
end
