%% S05_desirability_ranking.m
% Turn "balanced nutrition + texture" into an explicit, tunable score instead of ad-hoc classes.
% Derringer desirability: each trait -> d in [0,1]; overall D = geometric mean (a line that fails
% one criterion completely gets D = 0). Scenarios follow the outline:
%   S1 low GI + texture | S2 high RS + texture | S3 low GI + high RS + texture | S4 S3 + protein
% Texture desirability = closeness to the background control median (relative change within +/- texTol).
% Thresholds in S00_config are PLACEHOLDERS - agree them with your team before reporting rankings.
% Also reports the Pareto (non-dominated) set, which needs no weights at all.
clear; clc; S00_config;
load(fullfile(cfg.outDir, 'prepared.mat'));

if cfg.des.GIrelative
    dGI = desirability(T.d_GI, 'min', cfg.des.GIrel);   % change from background control median
else
    dGI = desirability(T.GI, 'min', cfg.des.GI);        % absolute GI
end
dRS  = desirability(T.RS, 'max', cfg.des.RS);
dT   = ones(height(T), 1);
for k = 1:numel(cfg.des.texTraits)
    dT = dT .* desirability(T.(['rel_' cfg.des.texTraits{k}]), 'target', cfg.des.texTol);
end
dTex = dT .^ (1 / numel(cfg.des.texTraits));
gm = @(varargin) prod(cat(2, varargin{:}), 2) .^ (1 / numel(varargin));

T.dGI = dGI;  T.dRS = dRS;  T.dTex = dTex;
T.D_S1_lowGI_tex       = gm(dGI, dTex);
T.D_S2_highRS_tex      = gm(dRS, dTex);
T.D_S3_lowGI_highRS_tex = gm(dGI, dRS, dTex);
if ismember('PC', T.Properties.VariableNames)
    T.dPC = desirability(T.PC, 'max', cfg.des.PC);
    T.D_S4_lowGI_highRS_tex_PC = gm(dGI, dRS, dTex, T.dPC);
end
Dcols = T.Properties.VariableNames(startsWith(T.Properties.VariableNames, 'D_S'));

%% Pareto front (minimise GI, maximise RS, minimise texture deviation)
ok  = ~isnan(T.GI) & ~isnan(T.RS) & ~isnan(dTex);
obj = [T.GI, -T.RS, 1 - dTex];                 % all "smaller is better"
idx = find(ok);  isPareto = false(height(T), 1);
for a = idx'
    dom = all(obj(idx, :) <= obj(a, :), 2) & any(obj(idx, :) < obj(a, :), 2);
    isPareto(a) = ~any(dom);
end
T.Pareto = isPareto;

%% Line-level ranking (edited lines only)
cols = [{'SrNo','EventID','Construct','Background','Generation','Zygosity','A1_class','A2_class', ...
         'GI','RS','TPA_Hardness','TPA_Adhesiveness','rel_TPA_Hardness','rel_TPA_Adhesiveness', ...
         'dGI','dRS','dTex'}, Dcols, {'Pareto','T0Event'}];
L = T(ok & T.IsControl == 0, cols);
L = sortrows(L, 'D_S3_lowGI_highRS_tex', 'descend');
writetable(L, fullfile(cfg.outDir, 'S05_line_ranking.csv'));
fprintf('Top 15 lines (scenario S3):\n'); disp(L(1:min(15, height(L)), {'EventID','Construct','Background','Generation','Zygosity','GI','RS','dTex','D_S3_lowGI_highRS_tex'}));

%% Event-level summary (selection unit for advancement)
E = groupsummary(L, 'T0Event', {'mean', 'max'}, Dcols);
E = sortrows(E, ['mean_' 'D_S3_lowGI_highRS_tex'], 'descend');
writetable(E, fullfile(cfg.outDir, 'S05_event_ranking.csv'));

%% Genotype-group summary ("which fine-tuned mutation class gives the balance?")
L.GenoGroup = string(L.Construct) + " | " + string(L.A1_class) + "/" + string(L.A2_class);
Gs = groupsummary(L, 'GenoGroup', {'mean', 'std'}, [{'GI', 'RS', 'rel_TPA_Hardness', 'rel_TPA_Adhesiveness'}, Dcols]);
Gs = sortrows(Gs, 'mean_D_S3_lowGI_highRS_tex', 'descend');
writetable(Gs, fullfile(cfg.outDir, 'S05_genotype_group_summary.csv'));
disp(Gs(:, {'GenoGroup', 'GroupCount', 'mean_GI', 'mean_RS', 'mean_rel_TPA_Hardness', 'mean_D_S3_lowGI_highRS_tex'}));

%% Plot: GI vs RS, colour = texture desirability, circles = Pareto set
fig = figure('Visible', 'off', 'Position', [100 100 640 500]);
e = ok & T.IsControl == 0;  c = ok & T.IsControl == 1;
scatter(T.RS(e), T.GI(e), 30, T.dTex(e), 'filled', 'MarkerFaceAlpha', 0.8); hold on
scatter(T.RS(c), T.GI(c), 45, 'k', 'x', 'LineWidth', 1.2);
p = T.Pareto & e;
scatter(T.RS(p), T.GI(p), 90, 'r', 'LineWidth', 1.2);
colormap(parula); cb = colorbar; cb.Label.String = 'texture desirability (1 = control-like)';
xlabel('Resistant starch (%)'); ylabel('Glycaemic index (in vitro)'); grid on
legend({'edited lines', 'WT/TC controls', 'Pareto set'}, 'Location', 'northeast');
title('Nutrition-texture trade-off');
exportgraphics(fig, fullfile(cfg.outDir, 'S05_tradeoff_pareto.png'), 'Resolution', 200); close(fig);
