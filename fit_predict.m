function [yhat, mdl, keep] = fit_predict(model, Xtr, ytr, Xte)
%FIT_PREDICT  Train one model on the training fold and predict the test fold.
%   model: 'EN' (elastic net), 'RF' (bagged trees), 'GPR', 'NN' (shallow ANN).
%   Missing predictors: trees use surrogate splits; other models get
%   training-fold median imputation (never uses test-fold information).
keep = std(Xtr, 0, 1, 'omitnan') > 0 & mean(isnan(Xtr), 1) < 0.5;   % drop constant / mostly-missing columns
Xtr = Xtr(:, keep);  Xte = Xte(:, keep);
if ~strcmp(model, 'RF')
    med = median(Xtr, 1, 'omitnan');
    for j = 1:size(Xtr, 2)
        Xtr(isnan(Xtr(:, j)), j) = med(j);
        Xte(isnan(Xte(:, j)), j) = med(j);
    end
end
switch model
    case 'EN'
        [B, S] = lasso(Xtr, ytr, 'Alpha', 0.5, 'CV', 5, 'Standardize', true);
        k = S.Index1SE;                              % more parsimonious than IndexMinMSE
        yhat = Xte * B(:, k) + S.Intercept(k);
        mdl = struct('B', B(:, k), 'Intercept', S.Intercept(k), 'Lambda', S.Lambda(k));
    case 'RF'
        p = size(Xtr, 2);
        t = templateTree('MinLeafSize', 5, 'Surrogate', 'on', ...
                         'NumVariablesToSample', max(1, round(p / 3)));
        mdl = fitrensemble(Xtr, ytr, 'Method', 'Bag', 'NumLearningCycles', 300, 'Learners', t);
        yhat = predict(mdl, Xte);
    case 'GPR'
        mdl = fitrgp(Xtr, ytr, 'KernelFunction', 'ardsquaredexponential', ...
                     'BasisFunction', 'constant', 'Standardize', true);
        yhat = predict(mdl, Xte);
    case 'NN'
        nRep = 5; yhat = zeros(size(Xte, 1), 1); mdl = cell(nRep, 1);
        for r = 1:nRep                                % average 5 nets: single small nets are unstable
            mdl{r} = fitrnet(Xtr, ytr, 'LayerSizes', 8, 'Activations', 'relu', ...
                             'Standardize', true, 'Lambda', 1e-2, 'IterationLimit', 1000);
            yhat = yhat + predict(mdl{r}, Xte) / nRep;
        end
    otherwise
        error('Unknown model %s', model);
end
end
