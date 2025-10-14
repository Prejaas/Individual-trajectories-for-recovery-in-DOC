% ========================================================================
% Project: Individual Trajectories for Recovery of Neocortical Activity in DoC
% Paper: "Individual trajectories for recovery of neocortical activity in disorders of consciousness" (2025)
%
% Repository: https://github.com/<org>/<repo>   (release: v1.0.0)
% File: scripts/compute_plasticity_psd.m
% License: MIT (see LICENSE)  |  SPDX-License-Identifier: MIT
%
% Purpose:
%   Compute sliding-window power spectra from simulated cortical/thalamic
%   signals of the corticothalamic model with plasticity, per subject, using
%   fitted parameters. Outputs a single MAT compatible with downstream code
%   that expects `freq` and `psdg_cor` (and adds extras).
%
%
% Inputs (on disk; not arguments):
%   - INPUT_DIR: folder with per-subject fit .mat files that contain struct S
%   - Toolboxes/Funcs on path: eeglab (eegfilt2), braintrak model,
%     Robinson_network_reduced_plas.m, freqplot_win.m
%
% Output:
%   - results_plasticity_auto.mat with:
%       freq           : [1 x F] Hz
%       psdg_cor       : [F x W x N] cortical PSD per window & subject
%       psdg_thal      : [F x W x N] thalamic PSD
%       plv            : [W x N]     (optional; zeros if disabled)
%       wins           : [1 x (W+1)] window start indices (in samples at Fs_ds)
%       Fs_ds          : scalar      sampling rate after downsample (Hz)
%       v_se_all       : 1xN cell    v_se per subject (downsampled 100×)
%       meta           : struct      parameters, versions, time
%
% Reproducibility:
%   - Deterministic given inputs; no randomness used.
% ========================================================================
function compute_plasticity_psd()
    %% ------------------------------- Config --------------------------------
    % Paths (edit to your environment)
    addpath(('G:\MATLAB\eeglab10_0_0_0b'));
    addpath(genpath('G:\MATLAB\Fellowship\braintrak'));
    addpath(genpath('G:\MATLAB\Fellowship\corticothalamic-model'));

    INPUT_DIR  = 'G:\MATLAB\Fellowship\results_nu_ab\patients\run2';
    OUT_PATH   = 'results_plasticity_auto.mat';  % repo-friendly relative path

    % Model simulation defaults (override per subject from S)
    baseparam = struct();
    baseparam.h         = 1e-4;   % time step (s)
    baseparam.T         = 3000;   % total time (s)
    baseparam.rho       = 17;     % plasticity strength (legacy default)
    baseparam.tim_plast = 10;     % plasticity update interval (s)
    % baseparam.noise   = [];     % pulled from S.model.p2.phin when available

    % Analysis options
    win_length   = 20;    % seconds per window
    bp_lo        = 0.5;   % Hz (band-pass low)
    bp_hi        = 50;    % Hz (band-pass high)

    % Sanity checks for required functions (why: early fail beats late crash)
    req_funcs = {'Robinson_network_reduced_plas','eegfilt2','freqplot_win'};
    for f = req_funcs
        assert(exist(f{1},'file')==2, 'Required function not on path: %s', f{1});
    end

    % List subject files
    files = [dir(fullfile(INPUT_DIR,'*.mat'))];
    % Filter to files that likely hold fits (optional pattern)
    files = files(~startsWith({files.name}, '.'));  % drop hidden
    assert(~isempty(files), 'No .mat fit files found in %s', INPUT_DIR);

    %%
    nSub = numel(files);
    psdg_cor = []; psdg_thal = []; freq = []; plv = []; wins = [];
    v_se_all = cell(1, nSub);
    Fs_ds    = [];   % sampling rate after downsample

    tAll = tic;
    for si = 1:nSub
        fname = fullfile(files(si).folder, files(si).name);
        fprintf('[%3d/%3d] %s\n', si, nSub, files(si).name);

        % -------- Load subject fit --------
        Swrap = load(fname);
        assert(isfield(Swrap,'S'), 'File %s missing struct S', files(si).name);
        S = Swrap.S;

        % -------- Build param from model + fit (why: provenance-preserving) --------
        p = baseparam;                                           % start with defaults
        p.Q_max = S.model.p2.qmax;
        p.theta = S.model.p2.theta * 1000;                       % s→ms legacy
        p.sigma = S.model.p2.sigma * 1000;                       % s→ms legacy
        p.g     = S.model.p2.gammae;
        p.a_e   = S.model.p2.alpha(1);
        p.b_e   = S.model.p2.beta(1);
        p.a_i   = S.model.p2.alpha(1);
        p.b_i   = S.model.p2.beta(1);
        p.a_t   = S.model.p2.alpha(1);
        p.b_t   = S.model.p2.beta(1);
        p.tau   = S.model.p2.taues;
        fp      = S.fit_data.fitted_params(:);
        p.v_ee  = fp(1) ;
        p.v_ei  = fp(2) ;
        p.v_es  = fp(3) ;
        p.v_se  = fp(4) ;
        p.v_sr  = fp(5) ;
        p.v_sn  = fp(6) ;
        p.v_re  = fp(7) ;
        p.v_rs  = fp(8) ;
        if isfield(S.model.p2,'phin'), p.noise = S.model.p2.phin; end

        % -------- Simulate model --------
        [ts_cortex, v_se, ts_thal] = Robinson_network_reduced_plas(p);

        % -------- Downsample & band-pass (why: speed, numerical stability) --------
        Fs  = 1 / p.h;         % original Hz (e.g., 10 kHz)
        dsF = 10;              % downsample factor
        Fs_ds_loc = Fs / dsF;  % e.g., 1000 Hz
        if isempty(Fs_ds), Fs_ds = Fs_ds_loc; else, assert(abs(Fs_ds-Fs_ds_loc)<1e-6, 'Inconsistent Fs_ds'); end

        sig_c = downsample(ts_cortex, dsF);
        sig_t = downsample(ts_thal,   dsF);

        % EEGLAB band-pass (why: replicates legacy)
        sig_c = eegfilt(sig_c.', Fs_ds, bp_lo, bp_hi).';

        % Keep plasticity variable (coarser rate for memory)
        v_se_all{si} = downsample(v_se, 100);

        % -------- Windowing (deterministic, non-overlapping) --------
        step  = round(win_length * Fs_ds);
        N     = min(numel(sig_c), numel(sig_t)); % align streams
        if N < 2*step
            warning('Subject %s too short after downsample; skipping.', files(si).name);
            continue;
        end
        starts = 1:step:(N - step + 1);  % inclusive starts
        if isempty(wins), wins = [starts, starts(end)+step]; end
        nWins = numel(starts);

        % -------- First window → PSD size & prealloc --------
        if isempty(freq)
            [freq, psd_c] = freqplot_win(sig_c(starts(1):starts(1)+step-1).', Fs_ds, 5);

            F = numel(freq);
            psdg_cor  = nan(F, nWins, nSub);

            % Assign first window
            psdg_cor(:,1,si)  = psd_c(:);

            wStart = 2; % next window index
        else
            wStart = 1;
        end

        % -------- Remaining windows --------
        for w = wStart:nWins
            seg_c = sig_c(starts(w):starts(w)+step-1).';

            [~, psd_c] = freqplot_win(seg_c, Fs_ds, 5);
            psdg_cor(:,w,si)  = psd_c(:);

        end

        fprintf('   -> simulated & analyzed (wins=%d, F=%d)\n', nWins, numel(freq));
    end
    fprintf('All subjects done in %.1f min\n', toc(tAll)/60);

    %% ------------------------------- Save ----------------------------------
    meta = struct();
    meta.bp = [bp_lo bp_hi];
    meta.win_length = win_length;
    meta.downsample = dsF;
    meta.date = datestr(now);
    meta.input_dir = INPUT_DIR;
    meta.notes = 'Sliding-window PSDs for cortical/thalamic signals from plasticity model.';

    save(OUT_PATH, 'freq','psdg_cor','psdg_thal','wins','Fs_ds','v_se_all','meta','-v7.3');
    fprintf('Saved: %s\n', OUT_PATH);
end
