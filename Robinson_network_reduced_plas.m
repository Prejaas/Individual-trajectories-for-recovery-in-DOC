function [timeseries, v_se, thalamus] = Robinson_network_reduced_plas(param)
% ========================================================================
% Project: Individual Trajectories for Recovery of Neocortical Activity in DoC
% File: models/Robinson_network_reduced_plas.m
% License: MIT (see LICENSE) | SPDX-License-Identifier: MIT
%
% Purpose:
%   Simulate the reduced corticothalamic model with activity-dependent
%   plasticity on v_se and additive white noise (Euler–Maruyama).
%
% Inputs (param fields; units in parentheses):
%   Required:
%     Q_max (1/s), theta (mV), sigma (mV), g (1/s), a_e (1/s), b_e (1/s),
%     a_t (1/s), b_t (1/s), v_ee..v_rs (mV·s), v_se (mV·s), q_std (1) noise SD,
%     tau (s) cortico-thalamic delay, h (s) step, T (s) duration,
%     rho (1) plasticity target, tim_plast (s) plasticity onset.
%   Optional:
%     rng_seed (int) for reproducible noise; default [] (non-deterministic)
%
% Outputs (column vectors, transients removed):
%   timeseries : cortical field φ_e(t)   (state 9)
%   v_se       : plastic synapse v_se(t)
%   thalamus   : thalamic state V_s^1(t) (state 6)
%
% Notes:
%   - Keeps your original equations and noise placement.
%   - Fixes delay indexing to use an integer number of steps consistently.
% ========================================================================

% --- validate & pull parameters ---
req = @(f) assert(isfield(param,f), 'Missing param.%s', f);
ff  = @(f) param.(f);

cellfun(req, {'Q_max','theta','sigma','g','a_e','b_e','a_t','b_t', ...
              'v_ee','v_ei','v_es','v_sr','v_se','v_rs','v_re', ...
              'q_std','tau','h','T','rho','tim_plast'});

Q_max = ff('Q_max');   theta = ff('theta'); sigma = ff('sigma');
g     = ff('g');       a_e   = ff('a_e');    b_e    = ff('b_e');
a_t   = ff('a_t');     b_t   = ff('b_t');
v_ee  = ff('v_ee');    v_ei  = ff('v_ei');   v_es   = ff('v_es');
v_sr  = ff('v_sr');    v_se0 = ff('v_se');   v_rs   = ff('v_rs');
v_re  = ff('v_re');    q_std = ff('q_std');  tau    = ff('tau');
h     = ff('h');       T     = ff('T');      rho    = ff('rho');
tim_plast = ff('tim_plast');

if isfield(param,'rng_seed') && ~isempty(param.rng_seed)
    rng(param.rng_seed, 'twister'); % why: reproducible noise if needed
end

% --- derived sizes & init ---
delay_steps = max(1, round(tau / h));         % integer delay in steps
N           = floor(T / h) + 2 + delay_steps; % robust length
v_se        = ones(1, N) * v_se0;             % plasticity trajectory
tau_isp_steps = 20 * (1 / h);                 % keep legacy scaling in steps
plast_onset   = max(1, round(tim_plast / h));

X = randn(10, N);            % state vector init
NOISE_STATE_IDX = 6;         % why: match original placement

% --- main loop (Euler–Maruyama) ---
for n = delay_steps+1 : N-1
    % noise (why: additive white noise on chosen state)
    noise = zeros(10,1);
    noise(NOISE_STATE_IDX) = a_e * b_e * q_std * randn(1);

    % plasticity update after onset (why: avoid early transients)
    if n > plast_onset
        % consistent delay indexing
        pre  = X(5, n);                       % V_s^1(t)
        post = X(9, n - delay_steps);         % φ_e(t - τ)
        dv   = synaptic_plas(pre, post, rho, tau_isp_steps);
        v_se(n+1) = v_se(n) + h * dv;
    end

    % one step with delayed state
    X(:,n+1) = X(:,n) + h * dynamics(X(:,n), X(:, n - delay_steps), v_se(n)) ...
                     + sqrt(h) * noise;
end

% --- outputs (drop 5 s transient; keep vectors as columns) ---
trim = max(1, round(5 / h));
timeseries = X(9, trim+1:end).';
thalamus   = X(6, trim+1:end).';
v_se       = v_se(:, trim+1:end).';

% ======================= nested local functions =========================
    function dX = dynamics(x, x_delay, v_se_now)
        % why: compact ODE/SDE RHS; inhibitory block omitted as in original
        dX       = zeros(10,1);
        dX(1)    = x(2);
        dX(2)    = a_e*b_e*(-x(1) + v_ee*x(9) + v_ei*(x(9)) + v_es*S(x_delay(5))) ...
                   - (a_e + b_e)*x(2);

        dX(5)    = x(6);
        dX(6)    = a_t*b_t*(-x(5) + v_se_now*x_delay(9) + v_sr*S(x(7))) ...
                   - (a_t + b_t)*x(6);

        dX(7)    = x(8);
        dX(8)    = a_t*b_t*(-x(7) + v_re*x_delay(9) + v_rs*S(x(5))) ...
                   - (a_t + b_t)*x(8);

        dX(9)    = x(10);
        dX(10)   = g^2*(-x(9) + S(x(1))) - 2*g*x(10);
    end

    function out = S(v)
        % why: standard sigmoid transfer
        out = Q_max ./ (1 + exp(-(v - theta) ./ sigma));
    end

    function d_vse = synaptic_plas(pre_v, post_phi_delay, rho_loc, tau_steps)
        % why: keep legacy scaling in 'steps' to match original code
        d_vse = -S(pre_v) .* (post_phi_delay - rho_loc) ./ tau_steps;
    end
end
