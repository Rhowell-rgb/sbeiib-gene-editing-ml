%% S06_mixed_model_baseline.m
% Interpretable statistical baseline that reviewers will expect next to the ML results:
% linear mixed model with allele dose effects, background and generation as fixed effects
% and T0 event as a random effect (siblings are not independent).
% If the ML models in S02/S03 do not beat this, report the LMM as the main result
% and the ML as a screening tool.
clear; clc; S00_config;
load(fullfile(cfg.outDir, 'prepared.mat'));

D = T(T.GenotypeKnown == 1 & T.Construct ~= 'IRS1625', :);   % coding-region edits + controls
D.CVGroup = categorical(D.CVGroup);
traits = [cfg.L2_targets, cfg.TPA, {'SEC_AM_total'}];
out = table();
for k = 1:numel(traits)
    y = ['d_' traits{k}];
    if ~ismember(y, D.Properties.VariableNames), continue; end
    f = sprintf('%s ~ dose_KO + dose_AS + dose_CEXT + dose_INF + Background + GenOrder + (1|CVGroup)', y);
    lme = fitlme(D, f);
    c = lme.Coefficients;
    tb = table(repmat(string(y), numel(c.Name), 1), string(c.Name), c.Estimate, c.SE, c.pValue, ...
               'VariableNames', {'Trait', 'Term', 'Estimate', 'SE', 'p'});
    out = [out; tb]; %#ok<AGROW>
    fprintf('\n%s  (n = %d, AIC = %.1f)\n', y, lme.NumObservations, lme.ModelCriterion.AIC);
    disp(tb(startsWith(tb.Term, 'dose_'), :));
end
writetable(out, fullfile(cfg.outDir, 'S06_LMM_dose_effects.csv'));
% Interpretation: Estimate for dose_KO = change per knock-out allele relative to control
% (in trait units), adjusted for background and generation.
