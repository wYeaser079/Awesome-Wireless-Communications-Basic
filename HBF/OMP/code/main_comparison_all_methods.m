%% main_comparison_all_methods.m
%
%  Comprehensive Comparison of Hybrid Beamforming Methods
%
%  This script compares multiple hybrid beamforming approaches:
%    1. Full-digital SVD (upper bound)
%    2. OMP hybrid with channel steering vectors as dictionary
%    3. OMP hybrid with oversampled ULA dictionary
%    4. OMP hybrid with oversampled UPA dictionary
%    5. Phase-reversal analog + ZF baseband
%    6. Beam steering (single-stream only)
%
%  For both single-user and multi-user scenarios.
%
%  References:
%    [1] El Ayach et al., IEEE TWC, 2014.
%    [2] ge99210/Hybrid-Precoding-Combining-
%    [3] le-liang/HybridPrecodingMassiveMIMO
%    [4] Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%    [5] metegenez/spatially-sparse-precoding
%    [6] MathWorks Phased Array Toolbox: omphybweights, ompdecomp
%
% =========================================================================

clear; clc; close all;

%% ===== PART 1: SINGLE-USER COMPARISON =====
fprintf('========================================\n');
fprintf('  PART 1: Single-User Hybrid Beamforming\n');
fprintf('========================================\n\n');

% System parameters
Nt    = 64;
Nr    = 16;
Ns    = 1;
N_RF  = 4;
Ncl   = 8;
Nray  = 10;
AS    = 7.5;  % degrees

SNR_dB  = -40:5:0;
SNR_lin = 10.^(SNR_dB / 10);
N_iter  = 100;

% Build oversampled UPA dictionary (for UPA channel)
Nt_H = sqrt(Nt); Nt_V = sqrt(Nt);
Nr_H = sqrt(Nr); Nr_V = sqrt(Nr);
[At_dict_upa, ~, ~] = build_dictionary_UPA(Nt_H, Nt_V, 2*Nt_H, 2*Nt_V);
[Ar_dict_upa, ~, ~] = build_dictionary_UPA(Nr_H, Nr_V, 2*Nr_H, 2*Nr_V);

% Preallocate
SE_opt = zeros(length(SNR_dB), 1);
SE_hyb_ch = zeros(length(SNR_dB), 1);
SE_hyb_dict = zeros(length(SNR_dB), 1);

for s = 1:length(SNR_dB)
    SNR = SNR_lin(s);
    t_opt = 0; t_ch = 0; t_dict = 0;

    for iter = 1:N_iter
        [H, At_ch, Ar_ch, ~] = generate_mmWave_channel_UPA(Nt, Nr, Ncl, Nray, AS);
        [U, S_mat, V] = svd(H);
        F_opt = V(:, 1:Ns);

        % Optimal
        W_opt = U(:, 1:Ns);
        t_opt = t_opt + compute_spectral_efficiency(H, F_opt, W_opt, Ns, SNR);

        % OMP with channel steering vectors
        [F_RF1, F_BB1] = OMP_hybrid_precoder(F_opt, At_ch, Ns, N_RF);
        [W_RF1, W_BB1] = OMP_hybrid_combiner(H, F_RF1, F_BB1, Ar_ch, Ns, N_RF, SNR);
        t_ch = t_ch + compute_spectral_efficiency(H, F_RF1*F_BB1, W_RF1*W_BB1, Ns, SNR);

        % OMP with oversampled UPA dictionary
        [F_RF2, F_BB2] = OMP_hybrid_precoder(F_opt, At_dict_upa, Ns, N_RF);
        [W_RF2, W_BB2] = OMP_hybrid_combiner(H, F_RF2, F_BB2, Ar_dict_upa, Ns, N_RF, SNR);
        t_dict = t_dict + compute_spectral_efficiency(H, F_RF2*F_BB2, W_RF2*W_BB2, Ns, SNR);
    end

    SE_opt(s)      = t_opt / N_iter;
    SE_hyb_ch(s)   = t_ch / N_iter;
    SE_hyb_dict(s) = t_dict / N_iter;

    fprintf('SNR=%3d dB: Opt=%.2f, OMP(ch)=%.2f, OMP(dict)=%.2f\n', ...
            SNR_dB(s), SE_opt(s), SE_hyb_ch(s), SE_hyb_dict(s));
end

figure(1); hold on; grid on;
plot(SNR_dB, SE_opt,      'k-s', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Optimal (Full Digital)');
plot(SNR_dB, SE_hyb_ch,   'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP (Channel Steering)');
plot(SNR_dB, SE_hyb_dict, 'r-^', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'OMP (Oversampled Dict)');
xlabel('SNR (dB)'); ylabel('Spectral Efficiency (bits/s/Hz)');
title(sprintf('SU: %dx%d, N_{RF}=%d, Ns=%d', Nt, Nr, N_RF, Ns));
legend('Location', 'northwest');

%% ===== PART 2: MULTI-USER COMPARISON =====
fprintf('\n========================================\n');
fprintf('  PART 2: Multi-User Hybrid Beamforming\n');
fprintf('========================================\n\n');

Nt_mu  = 64;
K      = 4;
N_RF_mu = K;
Ncl_mu  = 8;
d       = 0.5;
sigma_AS_mu = deg2rad(7.5);

SNR_dB_mu  = -20:5:10;
SNR_lin_mu = 10.^(SNR_dB_mu / 10);
N_iter_mu  = 200;

[A_dict_mu, ~] = build_dictionary_ULA(Nt_mu, 4*Nt_mu, d);

SE_fd_mu   = zeros(length(SNR_dB_mu), 1);
SE_omp_mu  = zeros(length(SNR_dB_mu), 1);
SE_pr_mu   = zeros(length(SNR_dB_mu), 1);

for s = 1:length(SNR_dB_mu)
    P = SNR_lin_mu(s);
    t_fd = 0; t_omp = 0; t_pr = 0;

    for iter = 1:N_iter_mu
        % Generate K-user channel
        H_mu = zeros(K, Nt_mu);
        for k = 1:K
            ang_k = (2*rand(Ncl_mu, 1) - 1) * pi;
            gain_k = (randn(Ncl_mu, 1) + 1j*randn(Ncl_mu, 1)) / sqrt(2);
            n_vec = (0:Nt_mu-1).';
            for c = 1:Ncl_mu
                a_t = (1/sqrt(Nt_mu)) * exp(1j * 2*pi*d * n_vec * sin(ang_k(c)));
                H_mu(k, :) = H_mu(k, :) + gain_k(c) * a_t.';
            end
        end
        H_mu = sqrt(Nt_mu / Ncl_mu) * H_mu;

        % Full-digital ZF
        W_fd = H_mu' / (H_mu * H_mu');
        for k = 1:K
            W_fd(:, k) = W_fd(:, k) / norm(W_fd(:, k));
        end
        for k = 1:K
            sig = (P/K) * abs(H_mu(k,:) * W_fd(:,k))^2;
            intf = 0;
            for j = 1:K
                if j ~= k
                    intf = intf + (P/K) * abs(H_mu(k,:) * W_fd(:,j))^2;
                end
            end
            t_fd = t_fd + log2(1 + sig / (intf + 1));
        end

        % OMP hybrid + ZF
        F_RF_omp = zeros(Nt_mu, 0);
        for k = 1:K
            fopt_k = H_mu(k,:)' / norm(H_mu(k,:));
            f_rf_k = [];
            f_res = fopt_k;
            for r = 1:1  % 1 RF chain per user
                psi = A_dict_mu' * f_res;
                [~, idx] = max(abs(psi).^2);
                f_rf_k = [f_rf_k, A_dict_mu(:, idx)];
                fbb_t = (f_rf_k' * f_rf_k) \ (f_rf_k' * fopt_k);
                f_res = fopt_k - f_rf_k * fbb_t;
                if norm(f_res) > 1e-10, f_res = f_res / norm(f_res); end
            end
            F_RF_omp = [F_RF_omp, f_rf_k];
        end
        H_eff = H_mu * F_RF_omp;
        F_BB_omp = H_eff' / (H_eff * H_eff');
        W_omp_mu = F_RF_omp * F_BB_omp;
        for k = 1:K
            W_omp_mu(:, k) = W_omp_mu(:, k) / norm(W_omp_mu(:, k));
        end
        for k = 1:K
            sig = (P/K) * abs(H_mu(k,:) * W_omp_mu(:,k))^2;
            intf = 0;
            for j = 1:K
                if j ~= k
                    intf = intf + (P/K) * abs(H_mu(k,:) * W_omp_mu(:,j))^2;
                end
            end
            t_omp = t_omp + log2(1 + sig / (intf + 1));
        end

        % Phase-reversal + ZF
        F_pr = (1/sqrt(Nt_mu)) * exp(1j * (-angle(H_mu)).');
        F_BB_pr = (H_mu * F_pr)' / ((H_mu * F_pr) * (H_mu * F_pr)');
        W_pr_mu = F_pr * F_BB_pr;
        for k = 1:K
            W_pr_mu(:, k) = W_pr_mu(:, k) / norm(W_pr_mu(:, k));
        end
        for k = 1:K
            sig = (P/K) * abs(H_mu(k,:) * W_pr_mu(:,k))^2;
            intf = 0;
            for j = 1:K
                if j ~= k
                    intf = intf + (P/K) * abs(H_mu(k,:) * W_pr_mu(:,j))^2;
                end
            end
            t_pr = t_pr + log2(1 + sig / (intf + 1));
        end
    end

    SE_fd_mu(s)  = real(t_fd) / N_iter_mu;
    SE_omp_mu(s) = real(t_omp) / N_iter_mu;
    SE_pr_mu(s)  = real(t_pr) / N_iter_mu;

    fprintf('SNR=%3d dB: FD-ZF=%.2f, OMP-ZF=%.2f, PR-ZF=%.2f\n', ...
            SNR_dB_mu(s), SE_fd_mu(s), SE_omp_mu(s), SE_pr_mu(s));
end

figure(2); hold on; grid on;
plot(SNR_dB_mu, SE_fd_mu,  'k-o', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Fully-Digital ZF');
plot(SNR_dB_mu, SE_omp_mu, 'r-*', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'OMP Hybrid + ZF');
plot(SNR_dB_mu, SE_pr_mu,  'b-^', 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', 'Phase-Reversal + ZF');
xlabel('SNR (dB)'); ylabel('Sum Spectral Efficiency (bits/s/Hz)');
title(sprintf('MU: Nt=%d, K=%d', Nt_mu, K));
legend('Location', 'northwest');

fprintf('\nAll simulations complete.\n');
