function [D, names] = one_hot(c, prefix)
%ONE_HOT  Dummy-code a categorical/string column. Keeps all levels (trees and
%   penalised models cope with the redundancy).
c = categorical(c);
lv = categories(c);
D = double(dummyvar(c));
D(isundefined(c), :) = NaN;
names = strcat(prefix, '_', matlab.lang.makeValidName(lv))';
end
