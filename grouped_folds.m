function G2F = grouped_folds(groupIdx, nFolds, nRepeats)
%GROUPED_FOLDS  Assign whole groups (T0 events) to folds, balanced by size.
%   G2F = grouped_folds(groupIdx, nFolds, nRepeats) returns an nGroups-by-nRepeats
%   matrix; the fold of observation i in repeat r is G2F(groupIdx(i), r).
%   All progeny of one T0 event therefore stay in the same fold (no sibling leakage).
nG   = max(groupIdx);
sz   = accumarray(groupIdx(:), 1, [nG 1]);
G2F  = zeros(nG, nRepeats);
for r = 1:nRepeats
    ord = randperm(nG);                        % random tie-breaking
    [~, k] = sort(sz(ord), 'descend');
    ord = ord(k);
    load_ = zeros(nFolds, 1);
    for g = ord
        cand = find(load_ == min(load_));
        f = cand(randi(numel(cand)));
        G2F(g, r) = f;
        load_(f) = load_(f) + sz(g);
    end
end
end
