function m = reg_metrics(y, yhat)
%REG_METRICS  Out-of-fold R^2 (vs. the overall mean), RMSE, MAE and Spearman rho.
ok = ~isnan(y) & ~isnan(yhat);
y = y(ok); yhat = yhat(ok);
m.n    = numel(y);
m.R2   = 1 - sum((y - yhat).^2) / sum((y - mean(y)).^2);
m.RMSE = sqrt(mean((y - yhat).^2));
m.MAE  = mean(abs(y - yhat));
m.rho  = corr(y, yhat, 'type', 'Spearman');
end
