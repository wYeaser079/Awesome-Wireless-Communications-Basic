%% main_SU_OMP_hybrid_beamforming.m
%
%  Single-User OMP-Based Hybrid Beamforming Simulation
%  (Reproduces Fig. 3 from El Ayach et al. 2014)
%
%  This script compares:
%    1. Optimal unconstrained precoding (full-digital SVD)
%    2. OMP-based hybrid precoding and combining
%    3. Beam steering (single-beam analog-only)
%
%  System: Narrowband mmWave MIMO with UPA antennas
%  Channel: Saleh-Valenzuela clustered channel model
%
%  References:
%    [1] O. El Ayach, S. Rajagopal, S. Abu-Surra, Z. Pi, and R. W. Heath,
%        "Spatially Sparse Precoding in Millimeter Wave MIMO Systems,"
%        IEEE Trans. Wireless Commun., vol. 13, no. 3, 2014.
%    [2] Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%    [3] metegenez/spatially-sparse-precoding
%
% =========================================================================

clear; clc; close all;

%% ===== System Parameters =====
Nt       = 64;        % Number of transmit antennas (8x8 UPA)
Nr       = 16;        % Number of receive antennas  (4x4 UPA)
Ns_vec   = [1, 2];    % Number of data streams to simulate
N_RF     = 4;         % Number of RF chains at both TX and RX

% Channel parameters
Ncl      = 8;         % Number of scattering clusters
Nray     = 10;        % Number of rays per cluster
AS       = 7.5;       % Angular spread (degrees)

% Simulation parameters
SNR_dB   = -40:5:0;   % SNR range in dB
SNR_lin  = 10.^(SNR_dB / 10);
N_iter   = 100;        % Number of Monte Carlo channel realizations

fprintf('=== Single-User OMP Hybrid Beamforming ===\n');
fprintf('Nt = %d, Nr = %d, N_RF = %d\n', Nt, Nr, N_RF);
fprintf('Ncl = %d, Nray = %d, AS = %.1f deg\n', Ncl, Nray, AS);
fprintf('Monte Carlo iterations: %d\n\n', N_iter);

%% ===== Main Simulation Loop =====
for ns_idx = 1:length(Ns_vec)
    Ns = Ns_vec(ns_idx);
    fprintf('--- Simulating Ns = %d ---\n', Ns);

    SE_optimal     = zeros(N_iter, length(SNR_lin));
    SE_hybrid      = zeros(N_iter, length(SNR_lin));
    SE_beamsteer   = zeros(N_iter, length(SNR_lin));

    for s = 1:length(SNR_lin)
        SNR = SNR_lin(s);

        for iter = 1:N_iter
            % ----- Generate mmWave channel (UPA, Saleh-Valenzuela) -----
            [H, At, Ar, Alpha] = generate_mmWave_channel_UPA(Nt, Nr, Ncl, Nray, AS);

            % ----- SVD of channel -----
            [U, S, V] = svd(H);
            F_opt = V(:, 1:Ns);           % Optimal precoder
            W_opt_mmse = ((1/sqrt(SNR)) * (F_opt' * H' * H * F_opt + ...
                          (Ns/SNR) * eye(Ns)) \ (F_opt' * H'))';
            % W_opt_mmse is Nr x Ns

            % ============================================================
            % A. Optimal Unconstrained Spectral Efficiency
            % ============================================================
            Rn_opt = W_opt_mmse' * W_opt_mmse;
            SE_optimal(iter, s) = abs(log2(det(eye(Ns) + (SNR/Ns) * ...
                (Rn_opt \ (W_opt_mmse' * H * F_opt * F_opt' * H' * W_opt_mmse)))));

            % ============================================================
            % B. OMP Hybrid Precoding
            % ============================================================
            [F_RF, F_BB] = OMP_hybrid_precoder(F_opt, At, Ns, N_RF);

            % ============================================================
            % C. OMP Hybrid Combining
            % ============================================================
            [W_RF, W_BB] = OMP_hybrid_combiner(H, F_RF, F_BB, Ar, Ns, N_RF, SNR);

            % Hybrid spectral efficiency
            Rn_hyb = W_BB' * W_RF' * W_RF * W_BB;
            SE_hybrid(iter, s) = abs(log2(det(eye(Ns) + (SNR/Ns) * ...
                (Rn_hyb \ (W_BB' * W_RF' * H * F_RF * F_BB * F_BB' * ...
                F_RF' * H' * W_RF * W_BB)))));

            % ============================================================
            % D. Beam Steering (only for Ns = 1)
            % ============================================================
            if Ns == 1
                [idx_t, idx_r] = find_beam_steering_vectors(H, At, Ar);
                F_BS = At(:, idx_t);
                W_BS = Ar(:, idx_r);
                Rn_bs = W_BS' * W_BS;
                SE_beamsteer(iter, s) = abs(log2(det(eye(Ns) + (SNR/Ns) * ...
                    (Rn_bs \ (W_BS' * H * F_BS * F_BS' * H' * W_BS)))));
            end
        end
    end

    % Average over Monte Carlo iterations
    SE_opt_avg  = mean(SE_optimal, 1);
    SE_hyb_avg  = mean(SE_hybrid, 1);
    SE_bs_avg   = mean(SE_beamsteer, 1);

    % ===== Plotting =====
    figure; hold on; grid on;
    plot(SNR_dB, SE_opt_avg, '-s', 'Color', [0 0.5 0], 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', sprintf('Optimal unconstrained (Ns=%d)', Ns));
    plot(SNR_dB, SE_hyb_avg, '-o', 'Color', [0 0.45 0.74], 'LineWidth', 2, 'MarkerSize', 8, ...
         'DisplayName', sprintf('Hybrid precoding/combining (Ns=%d)', Ns));
    if Ns == 1
        plot(SNR_dB, SE_bs_avg, '-d', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'MarkerSize', 8, ...
             'DisplayName', 'Beam steering (Ns=1)');
    end
    xlabel('SNR (dB)', 'FontSize', 14);
    ylabel('Spectral Efficiency (bits/s/Hz)', 'FontSize', 14);
    title(sprintf('SU Hybrid Beamforming: %dx%d, N_{RF}=%d, Ns=%d', Nt, Nr, N_RF, Ns), 'FontSize', 14);
    legend('Location', 'northwest', 'FontSize', 12);
    hold off;
end

fprintf('\nSimulation complete.\n');

%% ===== Local Helper Function: Beam Steering =====
function [idx_t, idx_r] = find_beam_steering_vectors(H, At, Ar)
% Exhaustive search for the best transmit/receive steering vector pair.
% Only applicable for single data stream (Ns = 1).
%
% Reference: Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm (findSteeringVector.m)

Npath = size(At, 2);
EffGain = zeros(Npath, Npath);

for t = 1:Npath
    for r = 1:Npath
        EffGain(t, r) = abs(Ar(:, r)' * H * At(:, t));
    end
end

[idx_t, idx_r] = find(EffGain == max(EffGain(:)), 1, 'first');
end
