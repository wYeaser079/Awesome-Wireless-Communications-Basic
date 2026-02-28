%% main_MU_OMP_hybrid_beamforming.m
%
%  Multi-User OMP-Based Hybrid Beamforming Simulation
%
%  This script simulates a multi-user mmWave MIMO downlink system with:
%    1. OMP-based analog RF precoder design (per-user)
%    2. Zero-forcing (ZF) digital baseband precoding on the effective channel
%    3. Comparison with fully-digital ZF precoding upper bound
%    4. Comparison with analog-only beam steering
%
%  System: BS with Nt antennas serves K single-antenna users
%  Channel: Saleh-Valenzuela clustered mmWave channel with ULA
%
%  References:
%    [1] El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%        MIMO Systems," IEEE TWC, 2014.
%    [2] Le Liang, "Low-complexity hybrid precoding in massive multiuser
%        MIMO systems," IEEE Wireless Communications Letters, 2014.
%    [3] le-liang/HybridPrecodingMassiveMIMO
%
% =========================================================================

clear; clc; close all;

%% ===== System Parameters =====
Nt       = 128;       % Number of BS transmit antennas (ULA)
K        = 4;         % Number of single-antenna users
Nr       = 1;         % Antennas per user (single-antenna)
N_RF     = K;         % Number of RF chains at BS (>= K for ZF)
Ns       = 1;         % Data streams per user

% Channel parameters
Ncl      = 10;        % Number of clusters per user
Nray     = 1;         % Rays per cluster (simplified: 1 ray = 1 path per cluster)
sigma_AS = deg2rad(7.5);  % Angular spread (radians)
d        = 0.5;       % Antenna spacing (wavelengths)

% Simulation parameters
SNR_dB   = -30:5:0;
SNR_lin  = 10.^(SNR_dB / 10);
N_iter   = 500;       % Monte Carlo iterations

% Dictionary for OMP
G_dict   = 4 * Nt;    % Oversampled dictionary size
[A_dict, ~] = build_dictionary_ULA(Nt, G_dict, d);

fprintf('=== Multi-User OMP Hybrid Beamforming ===\n');
fprintf('Nt = %d, K = %d, N_RF = %d\n', Nt, K, N_RF);
fprintf('Ncl = %d, Nray = %d\n', Ncl, Nray);
fprintf('Dictionary size: %d\n', G_dict);
fprintf('Monte Carlo iterations: %d\n\n', N_iter);

%% ===== Preallocate Results =====
SE_fulldig_ZF    = zeros(length(SNR_dB), 1);
SE_hybrid_OMP_ZF = zeros(length(SNR_dB), 1);
SE_hybrid_PR_ZF  = zeros(length(SNR_dB), 1);
SE_analog_only   = zeros(length(SNR_dB), 1);

%% ===== Main Simulation Loop =====
for s_idx = 1:length(SNR_dB)
    SNR = SNR_lin(s_idx);
    P = SNR;  % Total power (noise variance = 1)

    temp_fulldig  = 0;
    temp_omp_zf   = 0;
    temp_pr_zf    = 0;
    temp_analog   = 0;

    for iter = 1:N_iter
        % ================================================================
        % Generate multi-user mmWave channel (K x Nt)
        % Each row is a single-antenna user's channel
        % ================================================================
        [H, ~, At_all] = generate_MU_channel_ULA(Nt, K, Ncl, Nray, sigma_AS, d);

        % ================================================================
        % 1. Fully-Digital ZF Precoding (Upper Bound)
        % ================================================================
        W_ZF_full = H' / (H * H');                            % Nt x K
        % Normalize per-user
        for k = 1:K
            W_ZF_full(:, k) = W_ZF_full(:, k) / norm(W_ZF_full(:, k));
        end
        temp_fulldig = temp_fulldig + compute_MU_sum_rate(H, W_ZF_full, P/K, K);

        % ================================================================
        % 2. OMP Hybrid: Analog via OMP + Digital via ZF
        % ================================================================
        F_RF_omp = zeros(Nt, 0);
        rf_per_user = floor(N_RF / K);

        for k = 1:K
            hk = H(k, :);  % 1 x Nt
            % "Optimal" direction for single-antenna user = matched filter
            fopt_k = hk' / norm(hk);  % Nt x 1

            % OMP: select rf_per_user atoms from dictionary
            f_rf_k = [];
            f_res = fopt_k;
            for r = 1:rf_per_user
                psi = A_dict' * f_res;               % G x 1
                [~, idx] = max(abs(psi).^2);
                f_rf_k = [f_rf_k, A_dict(:, idx)];
                f_bb_temp = (f_rf_k' * f_rf_k) \ (f_rf_k' * fopt_k);
                f_res = fopt_k - f_rf_k * f_bb_temp;
                if norm(f_res) > 1e-10
                    f_res = f_res / norm(f_res);
                end
            end
            F_RF_omp = [F_RF_omp, f_rf_k];
        end

        % Effective channel through RF precoder
        H_eff_omp = H * F_RF_omp;                     % K x N_RF

        % ZF on effective channel
        F_BB_omp = H_eff_omp' / (H_eff_omp * H_eff_omp');  % N_RF x K

        % Total precoder and normalization
        W_omp = F_RF_omp * F_BB_omp;                  % Nt x K
        for k = 1:K
            W_omp(:, k) = W_omp(:, k) / norm(W_omp(:, k));
        end
        temp_omp_zf = temp_omp_zf + compute_MU_sum_rate(H, W_omp, P/K, K);

        % ================================================================
        % 3. Phase-Reversal Hybrid: Analog via conjugate phase + ZF
        %    (Le Liang's low-complexity method)
        % ================================================================
        F_RF_pr = zeros(Nt, K);
        for k = 1:K
            ph = -angle(H(k, :));
            F_RF_pr(:, k) = (1/sqrt(Nt)) * exp(1j * ph(:));
        end

        H_eff_pr = H * F_RF_pr;                        % K x K
        F_BB_pr = H_eff_pr' / (H_eff_pr * H_eff_pr');  % K x K

        W_pr = F_RF_pr * F_BB_pr;
        for k = 1:K
            W_pr(:, k) = W_pr(:, k) / norm(W_pr(:, k));
        end
        temp_pr_zf = temp_pr_zf + compute_MU_sum_rate(H, W_pr, P/K, K);

        % ================================================================
        % 4. Analog-only: DFT Beam Selection (Beam-Space MIMO)
        % ================================================================
        D = dftmtx(Nt);
        Hf = H * D;                                    % K x Nt in beam space
        beam_power = sum(abs(Hf).^2, 1);               % 1 x Nt
        [~, sorted_idx] = sort(beam_power, 'descend');
        F_beams = D(:, sorted_idx(1:K));                % Nt x K

        H_eff_beam = H * F_beams;
        F_BB_beam = H_eff_beam' / (H_eff_beam * H_eff_beam');
        W_beam = F_beams * F_BB_beam;
        for k = 1:K
            W_beam(:, k) = W_beam(:, k) / norm(W_beam(:, k));
        end
        temp_analog = temp_analog + compute_MU_sum_rate(H, W_beam, P/K, K);
    end

    SE_fulldig_ZF(s_idx)    = temp_fulldig / N_iter;
    SE_hybrid_OMP_ZF(s_idx) = temp_omp_zf / N_iter;
    SE_hybrid_PR_ZF(s_idx)  = temp_pr_zf / N_iter;
    SE_analog_only(s_idx)   = temp_analog / N_iter;

    fprintf('SNR = %3d dB: FD-ZF=%.2f, OMP-ZF=%.2f, PR-ZF=%.2f, Beam=%.2f bits/s/Hz\n', ...
            SNR_dB(s_idx), SE_fulldig_ZF(s_idx), SE_hybrid_OMP_ZF(s_idx), ...
            SE_hybrid_PR_ZF(s_idx), SE_analog_only(s_idx));
end

%% ===== Plotting =====
figure; hold on; grid on;
plot(SNR_dB, abs(SE_fulldig_ZF),    'k-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'Fully-Digital ZF');
plot(SNR_dB, abs(SE_hybrid_OMP_ZF), 'r-*', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'OMP Hybrid + ZF');
plot(SNR_dB, abs(SE_hybrid_PR_ZF),  'b-^', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'Phase-Reversal Hybrid + ZF');
plot(SNR_dB, abs(SE_analog_only),   'm-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'Beam-Space MIMO + ZF');
hold off;

xlabel('SNR (dB)', 'FontSize', 14);
ylabel('Sum Spectral Efficiency (bits/s/Hz)', 'FontSize', 14);
title(sprintf('MU Hybrid Beamforming: Nt=%d, K=%d, N_{RF}=%d', Nt, K, N_RF), 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12);

fprintf('\nSimulation complete.\n');

%% ===== Local Helper Functions =====

function [H, Gain, At] = generate_MU_channel_ULA(Nt, K, Ncl, Nray, sigma_AS, d)
% Generate simplified multi-user mmWave channel with ULA.
% Each user has a single antenna. K users, Nt BS antennas.
%
% Reference: le-liang/HybridPrecodingMassiveMIMO (GenChannelSimp.m)

if nargin < 6, d = 0.5; end

spread = pi;  % angular spread for cluster means
ang = (2 * rand(Ncl, K) - 1) * spread;  % cluster mean angles
Gain = (randn(Ncl, K) + 1j * randn(Ncl, K)) / sqrt(2);

n_vec = (0:Nt-1).';
At = zeros(Nt, Ncl, K);
H = zeros(K, Nt);

for k = 1:K
    for c = 1:Ncl
        % Generate ray angles with truncated Laplacian
        if Nray == 1
            ray_angles = ang(c, k);
        else
            u = rand(Nray, 1) - 0.5;
            b = sigma_AS / sqrt(2);
            ray_angles = ang(c, k) - b * sign(u) .* log(1 - 2 * abs(u));
        end

        for l = 1:Nray
            a_t = (1/sqrt(Nt)) * exp(1j * 2 * pi * d * n_vec * sin(ray_angles(l)));
            H(k, :) = H(k, :) + Gain(c, k) * a_t.';
        end
    end
    At(:, :, k) = (1/sqrt(Nt)) * exp(1j * 2 * pi * d * n_vec * sin(ang(:, k)).');
end

H = sqrt(Nt / (Ncl * Nray)) * H;

end

function R = compute_MU_sum_rate(H, W, P_per_user, K)
% Compute multi-user sum rate for single-antenna users.
%
%   R_k = log2(1 + P * |h_k * w_k|^2 / (sum_{j!=k} P * |h_k * w_j|^2 + 1))
%
% Reference: le-liang/HybridPrecodingMassiveMIMO (CalRate.m)

R = 0;
for k = 1:K
    hk = H(k, :);
    wk = W(:, k);

    sig_power  = P_per_user * abs(hk * wk)^2;
    int_power  = 0;
    for j = 1:K
        if j ~= k
            wj = W(:, j);
            int_power = int_power + P_per_user * abs(hk * wj)^2;
        end
    end

    R = R + log2(1 + sig_power / (int_power + 1));
end
R = real(R);

end
