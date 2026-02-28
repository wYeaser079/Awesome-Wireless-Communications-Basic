%% main_SU_OMP_hybrid_beamforming_ULA.m
%
%  Single-User OMP-Based Hybrid Beamforming with ULA and Dictionary Matrix
%
%  This script demonstrates the complete OMP hybrid beamforming pipeline
%  using ULA steering vectors and an oversampled dictionary matrix.
%
%  Components demonstrated:
%    - Saleh-Valenzuela channel model with ULA
%    - Dictionary matrix construction (oversampled angular grid)
%    - OMP hybrid precoder using dictionary
%    - OMP hybrid combiner using dictionary
%    - Spectral efficiency comparison
%
%  References:
%    [1] El Ayach et al., IEEE TWC, 2014.
%    [2] ge99210/Hybrid-Precoding-Combining-
%    [3] MathWorks omphybweights function
%
% =========================================================================

clear; clc; close all;

%% ===== System Parameters =====
Nt      = 64;         % Transmit antennas (ULA)
Nr      = 16;         % Receive antennas  (ULA)
Ns      = 2;          % Number of data streams
NtRF    = 4;          % TX RF chains
NrRF    = 4;          % RX RF chains
d       = 0.5;        % Antenna spacing (wavelengths)

% Channel parameters
Ncl     = 6;          % Number of clusters
Nray    = 3;          % Rays per cluster
sigma_AS = deg2rad(7.5);

% Dictionary parameters (oversampled by factor of 4)
Gt      = 4 * Nt;     % TX dictionary size
Gr      = 4 * Nr;     % RX dictionary size

% Simulation parameters
SNR_dB  = -30:5:0;
SNR_lin = 10.^(SNR_dB / 10);
N_iter  = 50;

fprintf('=== SU OMP Hybrid Beamforming with ULA & Dictionary ===\n');
fprintf('Nt=%d, Nr=%d, NtRF=%d, NrRF=%d, Ns=%d\n', Nt, Nr, NtRF, NrRF, Ns);
fprintf('TX dictionary: %d atoms, RX dictionary: %d atoms\n\n', Gt, Gr);

%% ===== Build Dictionary Matrices =====
[At_dict, ~] = build_dictionary_ULA(Nt, Gt, d);  % Nt x Gt
[Ar_dict, ~] = build_dictionary_ULA(Nr, Gr, d);  % Nr x Gr

%% ===== Preallocate =====
SE_optimal = zeros(length(SNR_dB), 1);
SE_hybrid  = zeros(length(SNR_dB), 1);
SE_fulldig = zeros(length(SNR_dB), 1);

%% ===== Main Loop =====
for s_idx = 1:length(SNR_dB)
    SNR = SNR_lin(s_idx);
    noisevar = 1 / SNR;

    temp_opt = 0;
    temp_hyb = 0;
    temp_fd  = 0;

    for iter = 1:N_iter
        % --- Generate ULA mmWave channel ---
        [H, At_ch, Ar_ch, ~] = generate_mmWave_channel_ULA(Nt, Nr, Ncl, Nray, sigma_AS, d);

        % --- SVD ---
        [U, S_mat, V] = svd(H);
        F_opt = V(:, 1:Ns);
        W_opt = U(:, 1:Ns);

        % --- Full-digital spectral efficiency ---
        temp_fd = temp_fd + compute_spectral_efficiency(H, F_opt, W_opt, Ns, SNR);

        % --- Optimal with MMSE combiner ---
        W_mmse = ((F_opt' * H' * H * F_opt + noisevar * Ns * eye(Ns)) \ (F_opt' * H'))';
        Rn = W_mmse' * W_mmse;
        temp_opt = temp_opt + abs(log2(det(eye(Ns) + (SNR/Ns) * ...
                   (Rn \ (W_mmse' * H * F_opt * F_opt' * H' * W_mmse)))));

        % --- OMP Hybrid Precoder (using DICTIONARY, not channel steering vectors) ---
        [F_RF, F_BB] = OMP_hybrid_precoder(F_opt, At_dict, Ns, NtRF);

        % --- OMP Hybrid Combiner (using DICTIONARY) ---
        [W_RF, W_BB] = OMP_hybrid_combiner(H, F_RF, F_BB, Ar_dict, Ns, NrRF, SNR);

        % --- Hybrid spectral efficiency ---
        F_total = F_RF * F_BB;
        W_total = W_RF * W_BB;
        temp_hyb = temp_hyb + compute_spectral_efficiency(H, F_total, W_total, Ns, SNR);
    end

    SE_fulldig(s_idx) = temp_fd / N_iter;
    SE_optimal(s_idx) = temp_opt / N_iter;
    SE_hybrid(s_idx)  = temp_hyb / N_iter;

    fprintf('SNR=%3d dB: Full-Digital=%.2f, Opt-MMSE=%.2f, OMP-Hybrid=%.2f bits/s/Hz\n', ...
            SNR_dB(s_idx), SE_fulldig(s_idx), SE_optimal(s_idx), SE_hybrid(s_idx));
end

%% ===== Plot Results =====
figure; hold on; grid on;
plot(SNR_dB, SE_fulldig, 'k--s', 'LineWidth', 2, 'MarkerSize', 8, ...
     'DisplayName', 'Full-Digital SVD');
plot(SNR_dB, SE_optimal, '-s', 'Color', [0 0.5 0], 'LineWidth', 2, 'MarkerSize', 8, ...
     'DisplayName', 'Optimal + MMSE Combiner');
plot(SNR_dB, SE_hybrid,  '-o', 'Color', [0 0.45 0.74], 'LineWidth', 2, 'MarkerSize', 8, ...
     'DisplayName', sprintf('OMP Hybrid (Dict=%dx)', Gt/Nt));
hold off;

xlabel('SNR (dB)', 'FontSize', 14);
ylabel('Spectral Efficiency (bits/s/Hz)', 'FontSize', 14);
title(sprintf('SU OMP Hybrid (ULA): %dx%d, N_{RF}=%d, Ns=%d', Nt, Nr, NtRF, Ns), 'FontSize', 14);
legend('Location', 'northwest', 'FontSize', 12);

fprintf('\nSimulation complete.\n');
