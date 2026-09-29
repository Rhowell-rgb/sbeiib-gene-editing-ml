%% S02_layer1_genotype_to_phenotype.m
% LAYER 1 (T0/T1 seed): can the edit (allele class, dose, protein length, promoter DOF loss)
% plus background predict starch structure (SEC) and texture (TPA)?
% Targets are control-relative (d_*) when cfg.useDelta = true.
% Output: results/L1_cv_results.csv, results/L1_oof.mat (out-of-fold predictions for stacking), plots.
clear; clc; S00_config;
load(fullfile(cfg.outDir, 'prepared.mat'));

baseCols = startsWith(namesG, {'BG_', 'CON_'});   % design-only baseline: background + construct
res = table();  L1oof = struct();

for t = 1:numel(cfg.L1_targets)
    tgt = cfg.L1_targets{t};
    if cfg.useDelta, yName = ['d_' tgt]; else, yName = tgt; end
    y = T.(yName);
    rows = ~isnan(y) & T.GenotypeKnown == 1;
    nGrp = numel(unique(groupIdx(rows)));
    fprintf('\n== %s : n = %d, T0 groups = %d\n', yName, nnz(rows), nGrp);

    % baseline: construct x background means (elastic net on one-hot only)
    [~, Mb] = run_grouped_cv('EN', XG(rows, baseCols), y(rows), groupIdx(rows), G2F);
    res = [res; summarize_cv(yName, 'Baseline(BG+CON)', 'EN', Mb, nGrp)]; %#ok<AGROW>

    bestR2 = -Inf; best = cfg.models{1}; bestPred = nan(height(T), 1);
    for m = 1:numel(cfg.models)
        mdlName = cfg.models{m};
        [oof, M] = run_grouped_cv(mdlName, XG(rows, :), y(rows), groupIdx(rows), G2F);
        res = [res; summarize_cv(yName, 'G', mdlName, M, nGrp)]; %#ok<AGROW>
        fprintf('   %-4s R2 = %5.2f (sd %4.2f)\n', mdlName, mean([M.R2]), std([M.R2]));
        full = nan(height(T), 1);  full(rows) = mean(oof, 2, 'omitnan');
        L1oof.(matlab.lang.makeValidName(tgt)).(mdlName) = full;
        if mean([M.R2]) > bestR2, bestR2 = mean([M.R2]); best = mdlName; bestPred = full; end
    end
    L1oof.(matlab.lang.makeValidName(tgt)).best = best;

    % observed vs predicted for the best model
    fig = figure('Visible', 'off', 'Position', [100 100 420 400]);
    gscatter(bestPred(rows), y(rows), T.Construct(rows)); hold on
    lim = [min([bestPred(rows); y(rows)]) max([bestPred(rows); y(rows)])];
    plot(lim, lim, 'k--'); axis equal tight; grid on
    xlabel('Out-of-fold prediction'); ylabel('Observed');
    title(sprintf('%s  (%s, R^2 = %.2f)', strrep(yName, '_', '\_'), best, bestR2));
    exportgraphics(fig, fullfile(cfg.outDir, ['L1_obs_pred_' yName '.png']), 'Resolution', 200);
    close(fig);
end

writetable(res, fullfile(cfg.outDir, 'L1_cv_results.csv'));
save(fullfile(cfg.outDir, 'L1_oof.mat'), 'L1oof');
disp(res);
% Read-out: if G beats Baseline(BG+CON), allele-level information (dose, class,
% protein length) predicts grain quality beyond "which construct / which background".
