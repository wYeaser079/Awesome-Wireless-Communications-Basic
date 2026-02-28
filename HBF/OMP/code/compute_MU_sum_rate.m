function [R_sum, R_per_user, SINR_per_user] = compute_MU_sum_rate(H_all, FRF, FBB, K, Ns_per_user, SNR_lin)
% COMPUTE_MU_SUM_RATE  Compute the sum-rate spectral efficiency and
%   per-user SINR for a multi-user MIMO downlink with hybrid precoding.
%
%   [R_sum, R_per_user, SINR_per_user] = compute_MU_sum_rate(H_all, FRF, FBB, K, Ns_per_user, SNR_lin)
%
%   Inputs:
%       H_all       - Cell array of K channel matrices, H_all{k} is Nr_k x Nt
%       FRF         - Nt x NtRF analog RF precoding matrix
%       FBB         - NtRF x (K*Ns_per_user) digital baseband precoding matrix
%       K           - Number of users
%       Ns_per_user - Number of data streams per user
%       SNR_lin     - SNR in linear scale (rho = P_total / sigma^2)
%
%   Outputs:
%       R_sum          - Sum-rate spectral efficiency (bits/s/Hz)
%       R_per_user     - K x 1 vector of per-user rates
%       SINR_per_user  - K x 1 vector of per-user SINR values (linear)
%                        (only for Ns_per_user = 1)
%
%   For Ns_per_user = 1 (single stream per user, common in MU-MIMO):
%   -----------------------------------------------------------------
%   The received signal at user k is:
%       y_k = H_k * FRF * fbb_k * s_k + sum_{j~=k} H_k * FRF * fbb_j * s_j + n_k
%
%   SINR for user k:
%       SINR_k = (rho/Ns) * |H_k * FRF * fbb_k|^2
%                / ( (rho/Ns) * sum_{j~=k} |H_k * FRF * fbb_j|^2 + 1 )
%
%   where fbb_k is the k-th column of FBB.
%
%   Rate for user k:
%       R_k = log2(1 + SINR_k)
%
%   Sum rate:
%       R_sum = sum_{k=1}^K R_k
%
%   For Ns_per_user > 1 (multiple streams per user):
%   -------------------------------------------------
%       R_k = log2(det(I + (rho/Ns) * Rn_k^{-1} * H_k*FRF*FBB_k*FBB_k'*FRF'*H_k'))
%   where Rn_k = (rho/Ns) * sum_{j~=k} H_k*FRF*FBB_j*FBB_j'*FRF'*H_k' + I
%
%   References:
%       [1] C. B. Peel et al., "A Vector-Perturbation Technique for
%           Near-Capacity Multiantenna Multiuser Communication," IEEE
%           Trans. Commun., 2005.
%       [2] A. Alkhateeb et al., "Limited Feedback Hybrid Precoding for
%           Multi-Layer Millimeter Wave MIMO," IEEE JSAC, 2015.

    Ns_total = K * Ns_per_user;
    F_total = FRF * FBB;    % Nt x Ns_total  (effective total precoder)

    R_per_user    = zeros(K, 1);
    SINR_per_user = zeros(K, 1);

    if Ns_per_user == 1
        % ==================================================================
        % CASE 1: Single stream per user
        % ==================================================================
        for k = 1:K
            H_k = H_all{k};    % Nr_k x Nt
            Nr_k = size(H_k, 1);

            % Effective channel for user k through each precoding column
            % h_eff_kj = H_k * f_j  (Nr_k x 1 vector for each column j)

            if Nr_k == 1
                % Single-antenna user: h_eff_kj is scalar
                h_eff = H_k * F_total;   % 1 x K row vector

                % Signal power (desired)
                signal_power = (SNR_lin / Ns_total) * abs(h_eff(k))^2;

                % Interference power (from other users)
                interf_idx = setdiff(1:K, k);
                interf_power = (SNR_lin / Ns_total) * sum(abs(h_eff(interf_idx)).^2);

                % SINR
                SINR_per_user(k) = signal_power / (interf_power + 1);

            else
                % Multi-antenna user: use strongest direction
                [U_k, ~, ~] = svd(H_k);
                w_k = U_k(:, 1);  % Receive combining vector (MRC on strongest mode)

                h_eff = w_k' * H_k * F_total;   % 1 x K effective channel

                signal_power = (SNR_lin / Ns_total) * abs(h_eff(k))^2;
                interf_idx = setdiff(1:K, k);
                interf_power = (SNR_lin / Ns_total) * sum(abs(h_eff(interf_idx)).^2);

                SINR_per_user(k) = signal_power / (interf_power + 1);
            end

            R_per_user(k) = log2(1 + SINR_per_user(k));
        end

    else
        % ==================================================================
        % CASE 2: Multiple streams per user
        % ==================================================================
        for k = 1:K
            H_k = H_all{k};
            Nr_k = size(H_k, 1);

            % Columns of FBB assigned to user k
            col_idx_k = (k-1)*Ns_per_user + (1:Ns_per_user);
            F_k = F_total(:, col_idx_k);   % Nt x Ns_per_user

            % Desired signal covariance at user k
            S_k = (SNR_lin / Ns_total) * (H_k * F_k) * (H_k * F_k)';

            % Interference-plus-noise covariance at user k
            Rn_k = eye(Nr_k);
            for j = 1:K
                if j ~= k
                    col_idx_j = (j-1)*Ns_per_user + (1:Ns_per_user);
                    F_j = F_total(:, col_idx_j);
                    Rn_k = Rn_k + (SNR_lin / Ns_total) * (H_k * F_j) * (H_k * F_j)';
                end
            end

            % Rate for user k
            R_per_user(k) = real(log2(det(eye(Nr_k) + Rn_k \ S_k)));
            R_per_user(k) = abs(R_per_user(k));
        end

        % SINR per user is not well-defined for multi-stream; set to NaN
        SINR_per_user = NaN(K, 1);
    end

    % Sum rate
    R_sum = sum(R_per_user);

end
