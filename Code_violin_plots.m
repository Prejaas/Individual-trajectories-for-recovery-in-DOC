% ========================================================================
% Project: Individual Trajectories for Recovery of Neocortical Activity in DoC
% Paper: "Individual trajectories for recovery of neocortical activity in disorders of consciousness" (2025)
%
% Authors:
%   Prejaas K.B. Tewarie^a,b,c,* , Romesh Abeysuriya^d,e , Rajanikant Panda^f,g ,
%   Pablo Núñez^f,g , Marie M. Vitello^f,g , Glenn van der Lande^f,g ,
%   Olivia Gosseries^f,g , Aurore Thibaut^f,g , Steven Laureys^a,f,g ,
%   Gustavo Deco^c,i , Jitka Annen^f,g
%
% * Correspondence: prbat@ulaval.ca
%
% Repository: https://github.com/Prejaas/Individual-trajectories-for-recovery-in-DOC
%
% Purpose: PLOT VIOLIN PLOTS OF ESTIMATED PARAMETERS
%
% Inputs:
%   - Fitted parameter .mat files (HC & patients): S.fit_data.fitted_params
%   - Demographics (one of):
%       • Anonym_clean_joined_demodata_and_lookup_v2.xlsx (private, not committed)
%       • demodata_minimal.mat (names, diag_cell)  <-- public fallback
%
% Outputs:
%   - Figure (HC–MCS–UWS violins) saved as .fig and .png
%   - CSV of p-values per parameter (optional)
%
% Dependencies:
%   - MATLAB R2023b+ (tested)
%   - FieldTrip (for I/O if used elsewhere)
%   - violinplot.m (https://github.com/bastibe/Violinplot-Matlab or equivalent)
%   - sigstar.m (optional significance bars)
%
%
% How to run:
%   1) Ensure paths to fit folders and demographics are set below.
%   2) Run this script; outputs will be written under SAVE_FIG_PATH*.
%
% ========================================================================


load("parameters_all.mat")

opts = struct( ...
  'thr_low',2.5,'thr_high',97.5,'alpha',0.05,'p_adjust','none', ...
  'param_names', {'G_{ee}','G_{ei}','G_{ese}','G_{esre}','G_{srs}','\alpha','\beta','t_0'}, ...
  'treat_zero_as_nan', true, ...
  'figure_id', 82, ...
  'save_path', 'figures/violin_params_3groups');
out = plot_violin_params_3groups(params_pat, params_hc, uws_ind, mcs_ind, opts);

% ========================================================================
% File: scripts/plot_violin_params_3groups.m
% Purpose: Publish-ready violin plots for HC vs MCS vs UWS (no EMCS, no MCS*)
% Methods: Wilcoxon rank-sum; percentile trimming (HC/MCS only); optional FDR
% Deps: violinplot.m, sigstar.m
% ========================================================================
function out = plot_violin_params_3groups(params_pat, params_hc, uws_ind, mcs_ind, opts)

% ---- options (defaults keep your legacy behavior) ----
if nargin < 5, opts = struct; end
thr_low   = g(opts,'thr_low',  2.5);       % percentile low
thr_high  = g(opts,'thr_high', 97.5);      % percentile high
alpha     = g(opts,'alpha',     0.05);
p_adjust  = lower(g(opts,'p_adjust','none'));  % 'none'|'fdr'
labels    = g(opts,'labels',   {'HC','MCS','UWS'});
param_names = g(opts,'param_names',[]);    % cellstr; default below
treat_zero_as_nan = g(opts,'treat_zero_as_nan', true);
fig_id    = g(opts,'figure_id', 82);
save_path = g(opts,'save_path','');
indg      = g(opts,'indg', []);            % optional patient reindex

% colors (match your palette)
colHC  = [1.0000 0.3673 0.4132];
colMCS = [0.0000 0.6637 1.0000];
colUWS = [0.9000 0.5637 0.4000];

% optional reorder of patients to mirror your 'indg'
if ~isempty(indg)
    params_pat = params_pat(:, indg);
    inv_idx = zeros(1,numel(indg)); inv_idx(indg) = 1:numel(indg);
    remap = @(idx) inv_idx(idx(idx>=1 & idx<=numel(inv_idx)));
    mcs_ind = remap(mcs_ind);
    uws_ind = remap(uws_ind);
end

P = size(params_pat,1);
if isempty(param_names)
    param_names = arrayfun(@(k)sprintf('Param %d',k), 1:P, 'uni',0);
    if P >= 8
        param_names{6} = '\alpha'; param_names{7} = '\beta'; param_names{8} = 't_0';
    end
end

pairs = {[1 2],[1 3],[2 3]}; % HC-MCS, HC-UWS, MCS-UWS

% outputs
out.pairs      = pairs;
out.labels     = labels;
out.pvals      = nan(P,3);
out.pvals_adj  = nan(P,3);

% figure
figure(fig_id); clf;
tiledlayout('flow','TileSpacing','compact','Padding','compact');

for p = 1:P
    % --- extract groups ---
    hc  = params_hc(p,:);
    mcs = take(params_pat(p,:), mcs_ind);
    uws = take(params_pat(p,:), uws_ind);

    % --- trim (HC/MCS only), keep UWS as-is (to match legacy) ---
    hc_rm  = rmoutliers(hc,  'percentiles', [thr_low thr_high]);
    mcs_rm = rmoutliers(mcs, 'percentiles', [thr_low thr_high]);
    uws_rm = uws;

    % --- build padded matrix: [HC MCS UWS] ---
    T = padcols({hc_rm, mcs_rm, uws_rm});
    if treat_zero_as_nan, T(T==0) = NaN; end

    % --- stats ----
    pv = nan(1,3);
    pv(1) = ranksum(hc_rm,  mcs_rm);   % 1-2
    pv(2) = ranksum(hc_rm,  uws_rm);   % 1-3
    pv(3) = ranksum(mcs_rm, uws_rm);   % 2-3

    if strcmp(p_adjust,'fdr')
        pv_adj = fdr_bh(pv);
    else
        pv_adj = pv;
    end
    out.pvals(p,:)     = pv;
    out.pvals_adj(p,:) = pv_adj;

    rw = {}; st = [];
    for k = 1:numel(pairs)
        if pv_adj(k) < alpha
            rw{end+1} = pairs{k}; %#ok<AGROW>
            st(end+1) = pv_adj(k); %#ok<AGROW>
        end
    end

    % --- plot ---
    nexttile;
    assert(exist('violinplot','file')==2, 'violinplot.m not found on path.');
    h = violinplot(T, 'ShowData', true, 'ShowMean', false);
    cols = {colHC, colMCS, colUWS};
    for ii = 1:3
        h(ii).ViolinPlot.LineWidth = 1.5;
        h(ii).ViolinColor          = cols{ii};
        h(ii).ViolinAlpha          = 0.45;
        if isprop(h(ii),'ShowData') && h(ii).ShowData
            set(h(ii).ScatterPlot, 'Marker','.', 'SizeData', 60);
        end
    end

    set(gca,'XTickLabel',labels,'FontSize',16,'Box','off','TickDir','in', ...
        'XColor',[.1 .1 .1],'YColor',[.1 .1 .1],'LineWidth',1,'FontName','Arial');
    xtickangle(35);
    title(param_names{p},'FontWeight','normal');
    xlim([0.6 3.4]);

    if exist('sigstar','file')==2 && ~isempty(rw)
        try, sigstar(rw, st); catch, end
    end
end

% optional save
if ~isempty(save_path)
    [d,~,~] = fileparts(save_path); if ~isempty(d) && ~exist(d,'dir'), mkdir(d); end
    savefig(gcf, [save_path '.fig']);
    exportgraphics(gcf, [save_path '.png'], 'Resolution', 300);
end

end

% ============================== helpers ==============================
function v = g(s, k, d)
if isfield(s,k) && ~isempty(s.(k)), v = s.(k); else, v = d; end
end

function arr = take(vec, idx)
if isempty(idx), arr = []; return; end
idx = idx(idx>=1 & idx<=numel(vec));
arr = vec(idx);
end

function T = padcols(cols)
n = max(cellfun(@numel, cols)); T = nan(n, numel(cols));
for i = 1:numel(cols), ci = cols{i}(:); T(1:numel(ci), i) = ci; end
end

function p = fdr_bh(p_raw)
p = p_raw; [ps, ix] = sort(p_raw); m = numel(p_raw);
q = ps .* (m ./ (1:m)); q = cummin(q(end:-1:1)); q = q(end:-1:1);
p(ix) = min(q, 1);
end
