function [oof, M] = run_grouped_cv(model, X, y, groupIdx, G2F)
%RUN_GROUPED_CV  Repeated grouped K-fold CV. Returns out-of-fold predictions
%   (nObs-by-nRepeats) and a struct array of metrics per repeat.
nR  = size(G2F, 2);
oof = nan(numel(y), nR);
for r = 1:nR
    fold = G2F(groupIdx, r);
    for f = unique(fold)'
        te = fold == f;  tr = ~te;
        if nnz(tr) < 10 || nnz(te) == 0, continue; end
        oof(te, r) = fit_predict(model, X(tr, :), y(tr), X(te, :));
    end
    M(r) = reg_metrics(y, oof(:, r)); %#ok<AGROW>
end
end
