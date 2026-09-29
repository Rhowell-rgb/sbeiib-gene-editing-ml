%% S01_prepare_data.m  —  load master table, normalise to controls, build feature matrices, fix CV folds
clear; clc; S00_config;

if ~isfile(cfg.dataFile)
    error('Data file not found: %s\nSee data/README.md (build it with python/build_master.py).', cfg.dataFile);
end
T = readtable(cfg.dataFile, 'TextType', 'string');
T.Background = categorical(T.Background);
T.Construct  = categorical(T.Construct);
fprintf('Loaded %d rows x %d columns\n', height(T), width(T));

%% 1) Control-relative traits (reference = MEDIAN of WT + TC controls of the same background;
%     median, because several controls are flagged as outliers in QC_issues)
traits = unique([cfg.SEC_parts, cfg.SEC_summary, {'SEC_SCAP_total','SEC_MCAP_total'}, ...
                 cfg.TPA, cfg.L2_targets], 'stable');
traits = traits(ismember(traits, T.Properties.VariableNames));
bgs = categories(T.Background);
ctrlMean = array2table(nan(numel(bgs), numel(traits)), 'VariableNames', traits, 'RowNames', bgs);
for j = 1:numel(traits)
    tr = traits{j};
    T.(['d_' tr])   = nan(height(T), 1);    % absolute difference from control
    T.(['rel_' tr]) = nan(height(T), 1);    % relative change (fraction) from control
    for b = 1:numel(bgs)
        inBg = T.Background == bgs{b};
        mu = median(T.(tr)(inBg & T.IsControl == 1), 'omitnan');
        ctrlMean{b, j} = mu;
        T.(['d_' tr])(inBg)   = T.(tr)(inBg) - mu;
        T.(['rel_' tr])(inBg) = (T.(tr)(inBg) - mu) ./ abs(mu);
    end
end
disp('Background control medians:'); disp(ctrlMean(:, intersect(traits, [cfg.TPA, cfg.L2_targets, {'SEC_AM_total'}], 'stable')));

%% 2) Genotype feature matrix XG
[bgD, bgN] = one_hot(T.Background, 'BG');
[csD, csN] = one_hot(T.Construct,  'CON');
XG = [T{:, cfg.G_numeric}, bgD, csD];
namesG = [cfg.G_numeric, bgN, csN];

%% 3) Layer-1 phenotype matrix XP1 (SEC + TPA)
if cfg.useCLR
    Z = clr_transform(T{:, cfg.SEC_parts});
    secN = strrep(cfg.SEC_parts, 'SEC_', 'clr_');
else
    Z = T{:, cfg.SEC_summary};  secN = cfg.SEC_summary;
end
if cfg.useDelta
    tpaX = T{:, strcat('rel_', cfg.TPA)};  tpaN = strcat('rel_', cfg.TPA);
else
    tpaX = T{:, cfg.TPA};                  tpaN = cfg.TPA;
end
XP1 = [Z, tpaX];
namesP1 = [secN, tpaN];

%% 4) Grouped CV folds (one assignment, reused by every script -> comparable & stackable)
[grpNames, ~, groupIdx] = unique(T.CVGroup);
G2F = grouped_folds(groupIdx, cfg.nFolds, cfg.nRepeats);
fprintf('%d CV groups (T0 events + controls) -> %d folds x %d repeats\n', numel(grpNames), cfg.nFolds, cfg.nRepeats);

%% 5) Quick overview
okL1 = T.GenotypeKnown == 1 & T.has_TPA == 1 & T.has_SEC == 1;
fprintf('Rows usable for Layer 1 (G + TPA + SEC): %d\n', nnz(okL1));
fprintf('Rows usable for Layer 2 (Layer 1 + GI/RS): %d\n', nnz(okL1 & T.has_GIRS == 1));
disp(groupsummary(T(T.IsControl == 0, :), {'Construct', 'Background'}, 'mean', {'GI', 'RS', 'TPA_Hardness', 'SEC_AM_total'}));

save(fullfile(cfg.outDir, 'prepared.mat'), 'T', 'XG', 'namesG', 'XP1', 'namesP1', ...
     'groupIdx', 'grpNames', 'G2F', 'ctrlMean', 'cfg');
fprintf('Saved %s\n', fullfile(cfg.outDir, 'prepared.mat'));
