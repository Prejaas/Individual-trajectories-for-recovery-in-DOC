% ========================================================================
% Project: Individual Trajectories for Recovery of Neocortical Activity in DoC
% Paper: "Individual trajectories for recovery of neocortical activity in disorders of consciousness" (2025)
%
% Authors:
%   Prejaas K.B. Tewarie^a,b,c,* , Romesh Abeysuriya^d,e , Rajanikant Panda^f,g ,
%   Pablo Núñez^f,g , Marie M. Vitello^f,g , Glenn van der Lande^f,g ,
%   Olivia Gosseries^f,g , Aurore Thibaut^f,g , Steven Laureys^a,f,g ,
%   Gustavo Deco^c,i , Jitka Annen^f,g
% * Correspondence: prbat@ulaval.ca
%
%
% Repository: https://github.com/Prejaas/Individual-trajectories-for-recovery-in-DOC
% 
% GOAL: create Figure 4 in the paper 
% ======================================================================
% File: scripts/02_analyze_plasticity_spectra.m
% ======================================================================

clear; clc; warning off;

BASE_RESULTS = 'D:\MATLAB\Fellowship\';
FREEZE_IN    = fullfile(BASE_RESULTS, 'plasticity_spectra.mat');
SPEC_MAT     = 'results_spec.mat';

L = load(FREEZE_IN); plasticity = L.plasticity; % outcome of plasticity simulations
met_idx = L.met_idx; % metabolic index
nSub  = numel(met_idx);

% Group membership (reuse exact subjects)
uws_idx_all  = plasticity.subjectIdx.uws(:);
mcs_idx_all  = plasticity.subjectIdx.mcs(:);
uws_mask_all = false(nSub,1); uws_mask_all(uws_idx_all) = true;
mcs_mask_all = false(nSub,1); mcs_mask_all(mcs_idx_all) = true;

uws_keep = plasticity.keepMask.uws(:);
mcs_keep = plasticity.keepMask.mcs(:);

% Colors
cUWS = [1.0000 0.3673 0.4132];
cMCS = [0.0282 0.7141 0.2944];
cFIT = [0.25   0.25   0.25];

% ========================= Figure (re-FOOOF plasticity spectra) =========================
freq     = plasticity.freq(:)';          
settings = struct();
settings.peak_width_limits = [1 14];
settings.max_n_peaks       = 2;
settings.min_peak_height   = 0.3;
settings.aperiodic_mode    = 'knee';

ap_uws = NaN(numel(uws_idx_all),1);
for k = 1:numel(uws_idx_all)
    sub_id  = uws_idx_all(k);
    psd_col = plasticity.psd.uws(:,k);
    f_range = [0.5, 30]; 
    foo = fooof(freq, psd_col, f_range, settings, true);
    ap_uws(k) = foo.aperiodic_params(3);
end
ap_mcs = NaN(numel(mcs_idx_all),1);
for k = 1:numel(mcs_idx_all)
    sub_id  = mcs_idx_all(k);
    psd_col = plasticity.psd.mcs(:,k);
    f_range = [0.5, 30]; 
    foo = fooof(freq, psd_col, f_range, settings, true);
    ap_mcs(k) = foo.aperiodic_params(3);
end

valval = NaN(nSub,1);
valval(uws_idx_all(uws_keep)) = ap_uws(uws_keep);
valval(mcs_idx_all(mcs_keep)) = ap_mcs(mcs_keep);

Xk = met_idx(:); Yk = valval(:);
in_groups = (uws_mask_all | mcs_mask_all);
keep = ~isnan(Xk) & ~isnan(Yk) & in_groups;

Xk = Xk(keep); Yk = Yk(keep);
uws_keep_global = uws_mask_all(keep);
mcs_keep_global = mcs_mask_all(keep);

figure(303); clf; hold on; box off;
set(gca,'FontSize',14,'TickDir','out','LineWidth',1);
scatter(Xk(uws_keep_global), Yk(uws_keep_global), 55, cUWS, 'filled', 'MarkerFaceAlpha',0.8);
scatter(Xk(mcs_keep_global), Yk(mcs_keep_global), 55, cMCS, 'filled', 'MarkerFaceAlpha',0.8);
p = polyfit(Xk, Yk, 1);
xx = linspace(min(Xk), max(Xk), 200); yy = polyval(p, xx);
plot(xx, yy, '-', 'Color', cFIT, 'LineWidth', 2);
[R303, P303] = corr(Xk, Yk, 'type','Pearson', 'rows','complete');
xlabel('Metabolic index'); ylabel('Aperiodic exponent (k)');
legend({'UWS','MCS'},'Location','southeast');
xl = xlim; yl = ylim;
text(xl(1)+0.02*range(xl), yl(2)-0.08*range(yl), sprintf('R = %.2f, p = %.3g', R303, P303), ...
    'Color', cFIT, 'FontSize', 14);
grid off; xlim([2 10]);


% Re-FOOOF empirical spectra for audit only (does not affect plot)
S2 = load(SPEC_MAT);
if ~isfield(S2,'freq') || ~isfield(S2,'mean_psdx'), error('results_spec.mat must have freq, mean_psdx'); end
freq_emp  = S2.freq(:)'; 
mean_psdx = S2.mean_psdx;
if isfield(S2,'erski') && ~isempty(S2.erski)
    try, mean_psdx(S2.erski(1:end-1),:) = []; catch, end
end
nSub_emp = size(mean_psdx,1);
nLoop    = min(nSub_emp, nSub);

settings_emp = struct('peak_width_limits',[1 14], 'max_n_peaks',2, ...
                      'min_peak_height',0.3, 'aperiodic_mode','knee');
f_range_emp  = [0.5, 45];
k_emp2 = NaN(nSub,1);
for s = 1:nLoop
    foo = fooof(freq_emp, mean_psdx(s,:), f_range_emp, settings_emp, true);
    k_emp2(s) = foo.aperiodic_params(3);
end

Y_viz = plasticity.fig202.Y_viz(:);

uws_mask = uws_mask_all;
mcs_mask = mcs_mask_all;
valid    = ~isnan(met_idx(:)) & ~isnan(Y_viz(:));
in_groups = (uws_mask | mcs_mask);
keep202  = valid & in_groups;

X2  = met_idx(keep202);
Y2  = Y_viz(keep202);
uw2 = uws_mask(keep202);
mc2 = mcs_mask(keep202);

figure(202); clf; hold on; box off;
set(gca,'FontSize',14,'TickDir','out','LineWidth',1);
scatter(X2(uw2), Y2(uw2),   55, cUWS, 'filled', 'MarkerFaceAlpha',0.8);
scatter(X2(mc2), Y2(mc2),   55, cMCS, 'filled', 'MarkerFaceAlpha',0.8);
p2 = polyfit(X2, Y2, 1);
xx2 = linspace(min(X2), max(X2), 200); yy2 = polyval(p2, xx2);
plot(xx2, yy2, '-', 'Color', cFIT, 'LineWidth', 2);
[R202, P202] = corr(X2, Y2, 'type','Pearson', 'rows','complete');
xlabel('Metabolic index'); ylabel('Aperiodic exponent (k)');  % label preserved
legend({'UWS','MCS'},'Location','southeast');
xl = xlim; yl = ylim;
text(xl(1)+0.02*range(xl), yl(2)-0.065*range(yl), sprintf('R = %.2f, p = %.3g', R202, P202), ...
     'Color', cFIT, 'FontSize', 14);
grid off; xlim([2 10]);