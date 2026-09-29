function row = summarize_cv(target, featSet, model, M, nGroups)
%SUMMARIZE_CV  One results-table row (mean and SD over CV repeats).
R2 = [M.R2]; RM = [M.RMSE]; rho = [M.rho];
row = table(string(target), string(featSet), string(model), M(1).n, nGroups, ...
    mean(R2), std(R2), mean(RM), mean(rho), ...
    'VariableNames', {'Target','FeatureSet','Model','n','nT0groups','R2_mean','R2_sd','RMSE_mean','Spearman_mean'});
end
