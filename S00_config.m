%% S00_config.m  —  shared settings for the SBEIIb ML pipeline
% Called at the top of every Sxx script. Edit here, not in the scripts.
% Requires: Statistics and Machine Learning Toolbox (R2022a or newer recommended;
%           shapley() needs R2021a+, fitrnet needs R2021a+).

cfg = struct();
cfg.root    = fileparts(mfilename('fullpath'));
if isempty(cfg.root), cfg.root = pwd; end
cfg.repoRoot = fileparts(cfg.root);                       % repository root (parent of matlab/)
cfg.dataFile = fullfile(cfg.repoRoot, 'data', 'SBEIIb_ML_master.csv');   % built by python/build_master.py
cfg.outDir   = fullfile(cfg.repoRoot, 'results');
if ~exist(cfg.outDir, 'dir'), mkdir(cfg.outDir); end
addpath(fullfile(cfg.root, 'helpers'));

cfg.seed     = 2026;
cfg.nFolds   = 5;      % grouped K-fold (groups = T0 events, column CVGroup)
cfg.nRepeats = 5;      % repeat CV with different group-to-fold assignments
cfg.useDelta = true;   % model traits as difference from the background control mean
cfg.models   = {'EN', 'RF', 'GPR', 'NN'};   % elastic net, random forest, Gaussian process, shallow ANN

% ---------------- genotype (G) features --------------------------------
% Background and Construct are one-hot encoded in S01. GBSS SNPs are NOT used:
% they are identical within each background (perfectly confounded).
cfg.G_numeric = {'GenOrder', 'IsBackcross', ...
    'dose_KO', 'dose_AS', 'dose_CEXT', 'dose_INF', 'dose_edited_CDS', 'dose_functional', ...
    'AAlen_min', 'AAlen_dev_mean', 'A1_indel_bp', 'A2_indel_bp', ...
    'Dom_CentralCatalytic_w', 'Dom_ActiveSiteI_w', 'Dom_CDbindV_w', 'Dom_loss_total_w', ...
    'DOF_lost_total', 'Prom_TS1_edited_alleles', 'Prom_TS2_edited_alleles', 'Prom_spanning_del'};
cfg.G_categorical = {'Background', 'Construct'};

% Generation / backcross status is confounded with assay batch or season for some traits
% (e.g. SEC amylose and RS). true  = use it as a design covariate (best prediction);
% false = edit-only model (cleaner attribution of effects to the edit). Run both and report both.
cfg.includeGeneration = true;
if ~cfg.includeGeneration
    cfg.G_numeric = setdiff(cfg.G_numeric, {'GenOrder', 'IsBackcross'}, 'stable');
end

% ---------------- Layer-1 phenotypes (cheap, early: T0/T1 seed) -------
cfg.SEC_parts = {'SEC_AM1','SEC_AM2','SEC_MCAP1','SEC_MCAP2','SEC_MCAP3','SEC_MCAP4','SEC_MCAP5', ...
                 'SEC_SCAP1','SEC_SCAP2','SEC_SCAP3','SEC_SCAP4','SEC_SCAP5'};
cfg.SEC_summary = {'SEC_AM_total','SEC_SCAP_MCAP_ratio','SEC_DP6to12_of_AP'};
cfg.TPA = {'TPA_Hardness','TPA_Adhesiveness','TPA_Cohesiveness','TPA_Springiness'};
cfg.useCLR = true;     % use centred log-ratio SEC fractions as features (compositional data)

% Layer-1 targets: can the edit predict starch structure and texture?
cfg.L1_targets = {'SEC_AM_total','SEC_SCAP_MCAP_ratio','SEC_DP6to12_of_AP', cfg.TPA{:}};

% Layer-2 targets (nutritional). Add AC, RVA_*, PC, AA_* here once those columns exist.
cfg.L2_targets = {'GI','RS'};

% ---------------- desirability (S05) -----------------------------------
cfg.des.GIrelative = true;  % true: score the GI CHANGE from the background control median (fairer across backgrounds)
cfg.des.GIrel = [-10 0];    %   d = 1 at <= -10 GI units vs control, 0 at >= 0
cfg.des.GI  = [55 65];      % used only if GIrelative = false: d = 1 at <= 55, 0 at >= 65 (absolute GI)
cfg.des.RS  = [1 6];        % d = 0 at <= 1 %, 1 at >= 6 % (larger is better)
cfg.des.texTol = 0.40;      % texture: d = 1 at control median, 0 at +/-40 % change (TPA replicate CV is 15-40 %)
cfg.des.texTraits = {'TPA_Hardness','TPA_Adhesiveness'};
cfg.des.PC  = [7 9];        % placeholder once protein content is available

% ---------------- SHAP (S04) -------------------------------------------
cfg.shapMaxQuery = 50;      % number of lines explained per target (shapley is slow; 200 for final figures)

rng(cfg.seed, 'twister');
