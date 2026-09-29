%% S03_layer2_nutrition.m
% LAYER 2: predict nutritional traits (GI, RS; later AC, RVA, PC, amino acids)
% from (A) genotype only, (B) early phenotypes only (SEC + TPA), (C) G + P1,
% (D) G + P1 + Layer-1 out-of-fold predictions (stacking).
% All feature sets are compared on the SAME rows and the SAME grouped folds.
clear; clc; S00_config;
load(fullfile(cfg.outDir, 'prepared.mat'));
load(fullfile(cfg.outDir, 'L1_oof.mat'));

% stacked features: RF out-of-fold Layer-1 predictions (never saw the row's own label)
stk = []; stkN = {};
f = fieldnames(L1oof);
for k = 1:numel(f)
    stk = [stk, L1oof.(f{k}).RF]; %#ok<AGROW>
    stkN{end + 1} = ['L1hat_' f{k}]; %#ok<AGROW>
end
bgCols = startsWith(namesG, 'BG_');

sets = struct( ...
  'name', {'A_G', 'B_P1', 'C_G+P1', 'D_G+P1+L1hat'}, ...
  'X',    {XG, [XP1, XG(:, bgCols)], [XG, XP1], [XG, XP1, stk]}, ...
  'names',{namesG, [namesP1, namesG(bgCols)], [namesG, namesP1], [namesG, namesP1, stkN]});

res = table();  L2oof = struct();
for t = 1:numel(cfg.L2_targets)
    tgt = cfg.L2_targets{t};
    if cfg.useDelta, yName = ['d_' tgt]; else, yName = tgt; end
    y = T.(yName);
    rows = ~isnan(y) & T.GenotypeKnown == 1 & T.has_TPA == 1 & T.has_SEC == 1;
    nGrp = numel(unique(groupIdx(rows)));
    fprintf('\n== %s : n = %d, T0 groups = %d\n', yName, nnz(rows), nGrp);
    for s = 1:numel(sets)
        for m = 1:numel(cfg.models)
            [oof, M] = run_grouped_cv(cfg.models{m}, sets(s).X(rows, :), y(rows), groupIdx(rows), G2F);
            res = [res; summarize_cv(yName, sets(s).name, cfg.models{m}, M, nGrp)]; %#ok<AGROW>
            fprintf('   %-14s %-4s R2 = %5.2f (sd %4.2f)\n', sets(s).name, cfg.models{m}, mean([M.R2]), std([M.R2]));
            full = nan(height(T), 1); full(rows) = mean(oof, 2, 'omitnan');
            L2oof.(matlab.lang.makeValidName(tgt)).(matlab.lang.makeValidName([sets(s).name '_' cfg.models{m}])) = full;
        end
    end
end
writetable(res, fullfile(cfg.outDir, 'L2_cv_results.csv'));
save(fullfile(cfg.outDir, 'L2_oof.mat'), 'L2oof', 'sets');
disp(sortrows(res, {'Target', 'R2_mean'}, {'ascend', 'descend'}));
% Read-out for the paper: if B or C clearly beats A, the cheap early assays
% (SEC/TPA on ~20 grains) carry information for GI/RS -> supports early selection.
