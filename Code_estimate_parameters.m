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
% Purpose:
%   Compute channel-averaged PSDs from multi-epoch EEG with FieldTrip,
%   then fit BrainTrak’s full model per subject on the mean PSD (DC bin removed).
%   Writes per-subject fit files and an aggregated PSD file for downstream use.
%
% GOAL: run on raw data (not provided due to ethical restrictions)
% Methods:
%   - Channel selection: HydroCel 257 scalp regions (C/T/P/O), electrode 257 excluded.
%   - PSD: Welch (window = win_sec * Fs, 50% overlap), averaged across selected channels.
%   - Optional downsample: 500→250 Hz for uniformity.
%   - Model: bt.model.full fit on (freq(2:end), mean_psd(2:end)).
%
% ========================================================================

%% ------------------------------- Config --------------------------------
clearvars; clc;

% Toolboxes (edit to your local clones)
addpath('D:\MATLAB\fieldtrip-master');
addpath(genpath('D:\MATLAB\Fellowship\braintrak'));
addpath(genpath('D:\MATLAB\Fellowship\corticothalamic-model'));
addpath(genpath('D:\MATLAB\Fellowship\HydroCel256-master'));  % for HCGSN257_scalp_regions.m

% I/O
RAW_DIR        = 'D:\Liege_EEG\patients';                    % input EEG files
OUTPUT_DIR     = 'D:\MATLAB\Fellowship\results\patients'; % per-subject fit .mat
AGGR_OUT_MAT   = fullfile('D:\MATLAB\Fellowship', 'results_spec.mat'); % aggregated PSD

% File filter (set to {'*.fif'} if you know your extension)
VALID_EXTS     = {'*.fif','*.edf','*.set','*.mat','*.*'};   % last is catch-all

% Processing options
WIN_SEC        = 5;     % Welch window length in seconds
OVERLAP_FRAC   = 0.5;   % 50% overlap
SKIP_SUB_IDX   = [27 29 44 88 102 106 141 168];  % skip by index (kept for compatibility)
SUPPRESS_WARN  = true;  % hide noisy warnings (why: some file readers are verbose)
WRITE_AGGREGATE = true; % write results_spec.mat at the end

% Sanity
if ~exist(RAW_DIR,'dir'); error('RAW_DIR not found: %s', RAW_DIR); end
if ~exist(OUTPUT_DIR,'dir'); mkdir(OUTPUT_DIR); end

if SUPPRESS_WARN, warning('off','all'); end

%% -------------------------- Channel selection --------------------------
% restrict to scalp C/T/P/O regions and drop ref 257 for stable PSDs
run('HCGSN257_scalp_regions.m');  % defines CChans, TChans, PChans, OChans
sel = [CChans TChans PChans OChans];
sel = sel(sel ~= 257);

%% ----------------------------- File listing ----------------------------
files = [];
for k = 1:numel(VALID_EXTS)
    files = [files; dir(fullfile(RAW_DIR, VALID_EXTS{k}))]; %#ok<AGROW>
end
% de-duplicate by name
[~, iu] = unique(string({files.name}),'stable');
files = files(iu);

nFiles = numel(files);
if nFiles == 0
    error('No input files found in %s matching %s', RAW_DIR, strjoin(VALID_EXTS,', '));
end
fprintf('Found %d files.\n', nFiles);

%% --------------------------- Prealloc holders --------------------------
mean_psdx = [];  % will size after first PSD
all_freq  = [];

params_pat_v2 = []; % will size after first fit (why: param count unknown)

%% --------------------------- Main processing ---------------------------
tStart = tic;
for sub = 1:nFiles
    if ismember(sub, SKIP_SUB_IDX)
        fprintf('[%3d/%3d] SKIP (in exclusion list): %s\n', sub, nFiles, files(sub).name);
        continue;
    end

    fname = fullfile(RAW_DIR, files(sub).name);
    fprintf('[%3d/%3d] Reading: %s\n', sub, nFiles, files(sub).name);

    % -------- Read raw EEG via FieldTrip --------
    try
        dataFT = ft_read_data(fname);
        hdrFT  = ft_read_header(fname);
    catch ME
        warning('FieldTrip read failed for %s: %s', files(sub).name, ME.message);
        continue;
    end

    % Flatten trials to continuous [nChan x nSamples] (why: uniform PSD handling)
    eeg = reshape(dataFT, size(dataFT,1), []);     % [channels x samples]
    if max(sel) > size(eeg,1)
        warning('Channel selection exceeds data channels for %s. Skipping.', files(sub).name);
        continue;
    end
    eeg = eeg(sel, :);

    Fs = hdrFT.Fs;
    % Downsample 500→250 Hz (why: match earlier pipeline; avoids window drift)
    if abs(Fs - 500) < 1e-6
        eeg = downsample(eeg.', 2).';  % keep memory friendly
        Fs  = Fs / 2;
    end

    % -------- PSD via Welch (per channel), then mean across channels -----
    winSamples   = round(WIN_SEC * Fs);
    overlap      = round(OVERLAP_FRAC * winSamples);

    % First channel to establish freq grid
    [px, freq] = pwelch(eeg(1,:), winSamples, overlap, [], Fs);
    nFreq = numel(freq);
    psdx_temp = zeros(nFreq, size(eeg,1), 'like', px);
    psdx_temp(:,1) = px;

    for ch = 2:size(eeg,1)
        psdx_temp(:,ch) = pwelch(eeg(ch,:), winSamples, overlap, [], Fs);
    end

    psd_mean = mean(psdx_temp, 2).';  % row vector

    % Grow holders
    if isempty(mean_psdx)
        mean_psdx = nan(nFiles, nFreq);
        all_freq  = freq(:).';
    end
    mean_psdx(sub, :) = psd_mean;

    % -------- BrainTrak fit on mean PSD (DC removed) --------------------
    freqs_fit = all_freq; freqs_fit(1) = [];             % drop DC
    psd_fit   = psd_mean;  psd_fit(1)   = [];

    try
        a = bt.fit(bt.model.full, freqs_fit, psd_fit.'); % column PSD per bt.fit
    catch ME
        warning('BrainTrak fit failed for %s: %s', files(sub).name, ME.message);
        continue;
    end

    % Capture fitted params (why: downstream QC)
    this_params = a.fit_data.fitted_params;
    if isempty(params_pat_v2)
        params_pat_v2 = nan(numel(this_params), nFiles);
    end
    params_pat_v2(:, sub) = this_params;

    % -------- Persist subject-level fit ---------------------------------
    [~, base, ~] = fileparts(files(sub).name);
    outName = fullfile(OUTPUT_DIR, sprintf('full_track_%s.mat', base));
    S = struct(a);             % compact serialization of fit object
    A = struct(a.model);       % model settings snapshot
    save(outName, 'S', 'A');

    fprintf('   -> Saved: %s\n', outName);
end
elapsedMin = toc(tStart)/60;
fprintf('Done. Elapsed: %.1f min\n', elapsedMin);

%% ------------------------ Aggregate outputs ---------------------------
if WRITE_AGGREGATE
    % Remove fully empty rows (why: if some subs were skipped/failed)
    keepRows = any(isfinite(mean_psdx), 2);
    mean_psdx_out = mean_psdx(keepRows, :);
    erski = SKIP_SUB_IDX; %#ok<NASGU> % kept for compatibility with prior code

    save(AGGR_OUT_MAT, 'mean_psdx_out', 'all_freq', 'erski', '-v7.3');
    fprintf('Wrote aggregate PSD: %s  (subjects=%d, nFreq=%d)\n', AGGR_OUT_MAT, size(mean_psdx_out,1), numel(all_freq));
end
