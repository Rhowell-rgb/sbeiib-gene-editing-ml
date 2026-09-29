%% S04_shap_interpretation.m
% Fit a random forest on all Layer-2 rows (feature set C: genotype + SEC + TPA) and explain it
% with SHAP (MATLAB shapley, R2021a+) and out-of-bag permutation importance.
% Question addressed: which edit features / structural traits push GI down and RS up?
% NOTE: SHAP explains the MODEL, not causality; report it alongside the CV R^2 from S03
%       (interpret only targets whose CV R^2 is clearly > 0).
clear; clc; S00_config;
load(fullfile(cfg.outDir, 'prepared.mat'));

X = [XG, XP1];  names = [namesG, namesP1];
maxQuery = cfg.shapMaxQuery;                  % set in S00_config (shapley is slow)

for t = 1:numel(cfg.L2_targets)
    tgt = cfg.L2_targets{t};
    if cfg.useDelta, yName = ['d_' tgt]; else, yName = tgt; end
    y = T.(yName);
    rows = find(~isnan(y) & T.GenotypeKnown == 1 & T.has_TPA == 1 & T.has_SEC == 1);
    Xr = X(rows, :);
    keep = std(Xr, 0, 1, 'omitnan') > 0;
    Xr = Xr(:, keep);  nm = names(keep);

    tpl = templateTree('MinLeafSize', 5, 'Surrogate', 'on', 'NumVariablesToSample', max(1, round(numel(nm) / 3)));
    mdl = fitrensemble(Xr, y(rows), 'Method', 'Bag', 'NumLearningCycles', 500, ...
                       'Learners', tpl, 'PredictorNames', nm);

    % --- permutation importance (out-of-bag) ------------------------------
    imp = oobPermutedPredictorImportance(mdl);
    [~, o] = sort(imp, 'descend');
    fig = figure('Visible', 'off', 'Position', [100 100 520 520]);
    k = min(15, numel(o));
    barh(imp(o(k:-1:1))); yticks(1:k); yticklabels(strrep(nm(o(k:-1:1)), '_', '\_'));
    xlabel('OOB permutation importance'); title(strrep(yName, '_', '\_'));
    exportgraphics(fig, fullfile(cfg.outDir, ['L2_permImp_' yName '.png']), 'Resolution', 200); close(fig);

    % --- SHAP values, one query point at a time (works R2021a onwards) ----
    Xfill = Xr;                                  % shapley needs complete query rows
    med = median(Xr, 1, 'omitnan');
    for j = 1:size(Xfill, 2), Xfill(isnan(Xfill(:, j)), j) = med(j); end
    q = rows; if numel(q) > maxQuery, q = rows(randperm(numel(rows), maxQuery)); end
    [~, qi] = ismember(q, rows);
    S = nan(numel(qi), numel(nm));
    explainer = shapley(mdl, Xfill);
    for i = 1:numel(qi)
        explainer = fit(explainer, Xfill(qi(i), :));
        S(i, :) = explainer.ShapleyValues.ShapleyValue';
        if mod(i, 25) == 0, fprintf('  SHAP %s: %d/%d\n', yName, i, numel(qi)); end
    end
    shapTbl = array2table(S, 'VariableNames', matlab.lang.makeValidName(nm));
    shapTbl = [T(q, {'SrNo', 'EventID', 'Construct', 'Background', 'Generation', 'Zygosity'}), shapTbl];
    writetable(shapTbl, fullfile(cfg.outDir, ['L2_SHAP_' yName '.csv']));

    % --- beeswarm-style summary plot --------------------------------------
    [~, o] = sort(mean(abs(S), 1), 'descend');  k = min(15, numel(o));
    fig = figure('Visible', 'off', 'Position', [100 100 620 560]); hold on
    for r = 1:k
        j = o(r);  v = Xfill(qi, j);
        c = (v - min(v)) ./ max(eps, max(v) - min(v));       % feature value scaled 0-1
        swarmchart(S(:, j), (k - r + 1) * ones(numel(qi), 1), 12, c, 'filled', ...
                   'YJitter', 'density', 'YJitterWidth', 0.7, 'MarkerFaceAlpha', 0.7);
    end
    colormap(parula); cb = colorbar; cb.Label.String = 'feature value (low \rightarrow high)';
    yticks(1:k); yticklabels(strrep(nm(o(k:-1:1)), '_', '\_'));
    xline(0, 'k-'); xlabel(sprintf('SHAP value (effect on %s)', strrep(yName, '_', '\_'))); grid on
    exportgraphics(fig, fullfile(cfg.outDir, ['L2_SHAP_beeswarm_' yName '.png']), 'Resolution', 200); close(fig);

    % --- SHAP aggregated by genotype group (answers "which mutation group") --
    grp = string(T.Construct(q)) + " | " + string(T.Zygosity(q));
    tot = sum(S(:, startsWith(nm, {'dose_', 'AAlen', 'A1_', 'A2_', 'Dom_', 'DOF', 'Prom_'})), 2);
    G = groupsummary(table(grp, tot, 'VariableNames', {'Group', 'SHAP_editFeatures'}), 'Group', {'mean', 'std'});
    writetable(G, fullfile(cfg.outDir, ['L2_SHAP_byGenotypeGroup_' yName '.csv']));
    disp(G);
end
