function [F_RF, F_BB, W_ZF] = OMP_hybrid_precoder_MU(H, At, K, Ns_per_user, N_RF, SNR)
% OMP_HYBRID_PRECODER_MU  Multi-user OMP-based hybrid precoding with
%   zero-forcing (ZF) digital precoding on the effective channel.
%
%   [F_RF, F_BB, W_ZF] = OMP_hybrid_precoder_MU(H, At, K, Ns_per_user, N_RF, SNR)
%
%   Inputs:
%       H            - K*Nr x Nt  stacked channel matrix (or K x Nt for
%                      single-antenna users)
%       At           - Nt x G     transmit dictionary matrix
%       K            - Number of users
%       Ns_per_user  - Scalar or K x 1 vector of data streams per user
%       N_RF         - Number of RF chains at the BS
%       SNR          - Signal-to-noise ratio (linear scale, for normalization)
%
%   Outputs:
%       F_RF    - Nt x N_RF   analog/RF precoding matrix (from OMP)
%       F_BB    - N_RF x Ns_total  digital/baseband precoding matrix (from ZF)
%       W_ZF    - Nt x Ns_total   total hybrid precoder F_RF * F_BB
%
%   Algorithm:
%   -------------------------------------------------------
%   Step 1: Compute the analog RF precoder F_RF using OMP.
%           For each user, use OMP to select the best RF beamforming
%           directions from the dictionary At.
%
%   Step 2: Compute the effective channel: H_eff = H * F_RF  (K x N_RF)
%
%   Step 3: Apply ZF digital precoding on the effective channel:
%           F_BB = H_eff^H * (H_eff * H_eff^H)^{-1}
%           Then normalize per-user columns.
%
%   There are two main approaches for multi-user RF precoder design:
%
%   Approach A (Per-user OMP): Run OMP independently per user's channel
%       to find each user's best steering directions, then concatenate.
%
%   Approach B (Joint OMP): Run OMP on the aggregate channel to find
%       N_RF columns of At that span the multi-user channel space.
%
%   This function implements Approach A (per-user) by default, which is
%   more standard in the literature.
%
%   References:
%       [1] El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%           MIMO Systems," IEEE TWC, 2014.
%       [2] Le Liang, HybridPrecodingMassiveMIMO (mainCompareScheme_mmWave.m)
%       [3] A. Alkhateeb et al., "Hybrid Precoding for Millimeter Wave
%           Cellular Systems with Partial Channel Knowledge," ITA, 2013.
%
% =========================================================================

if isscalar(Ns_per_user)
    Ns_per_user = Ns_per_user * ones(K, 1);
end
Ns_total = sum(Ns_per_user);

[Nr_total, Nt] = size(H);
Nr_per_user = Nr_total / K;   % antennas per user (assumed equal)

% --- Determine RF chains allocated per user ---
% Simple equal allocation: N_RF_per_user(k) = floor(N_RF / K)
% with remainder distributed to first users
N_RF_per_user = floor(N_RF / K) * ones(K, 1);
remainder = N_RF - sum(N_RF_per_user);
N_RF_per_user(1:remainder) = N_RF_per_user(1:remainder) + 1;

% ====================================================================
%  APPROACH A: Per-user OMP for analog precoder
% ====================================================================

F_RF = zeros(Nt, 0);

for k = 1:K
    % Extract user k's channel
    row_idx = (k-1)*Nr_per_user + (1:Nr_per_user);
    Hk = H(row_idx, :);     % Nr_per_user x Nt

    % Optimal unconstrained precoder for user k (from SVD)
    [~, ~, Vk] = svd(Hk);
    Fopt_k = Vk(:, 1:Ns_per_user(k));   % Nt x Ns_k

    % OMP to find the RF precoding columns for user k
    Frf_k = [];
    F_res = Fopt_k;
    for r = 1:N_RF_per_user(k)
        Psi = At' * F_res;
        [~, idx] = max(diag(Psi * Psi'));
        Frf_k = [Frf_k, At(:, idx)];
        F_bb_temp = (Frf_k' * Frf_k) \ (Frf_k' * Fopt_k);
        F_res = Fopt_k - Frf_k * F_bb_temp;
        if norm(F_res, 'fro') > 1e-10
            F_res = F_res / norm(F_res, 'fro');
        end
    end

    F_RF = [F_RF, Frf_k];   % Nt x (accumulated RF chains)
end

% ====================================================================
%  ZF digital precoding on effective channel
% ====================================================================

% Effective channel: H_eff = H * F_RF    (K*Nr_per_user x N_RF)
H_eff = H * F_RF;

% Zero-forcing (pseudo-inverse)
F_BB = H_eff' / (H_eff * H_eff');   % N_RF x K*Nr_per_user

% Select only the required Ns_total streams (for single-antenna users, K = Ns_total)
if size(F_BB, 2) > Ns_total
    % If users have multiple antennas, use block diagonalization instead
    F_BB = block_diag_precoder(H_eff, Ns_per_user);
end

% --- Power normalization (per-user equal power) ---
W_ZF = F_RF * F_BB;       % Nt x Ns_total total precoder

% Normalize so that each user's columns have unit norm
offset = 0;
for k = 1:K
    cols = offset + (1:Ns_per_user(k));
    W_ZF(:, cols) = W_ZF(:, cols) / norm(W_ZF(:, cols), 'fro') * sqrt(Ns_per_user(k));
    offset = offset + Ns_per_user(k);
end

% Recompute F_BB from normalized total precoder
F_BB = (F_RF' * F_RF) \ (F_RF' * W_ZF);

end

% =========================================================================
function F_BD = block_diag_precoder(H, Ns_per_user)
% Block diagonalization precoder for multi-antenna users.
%
% Reference: Le Liang, HybridPrecodingMassiveMIMO (CalBDPrecoder.m)
%            Spencer et al., "Zero-Forcing Methods for Downlink Spatial
%            Multiplexing in Multiuser MIMO Channels," IEEE TSP, 2004.

[NR, Nt] = size(H);
K = length(Ns_per_user);
Nr = NR / K;
F_BD = [];

for k = 1:K
    Gk = H((k-1)*Nr+1:k*Nr, :);
    G_tilde = H;
    G_tilde((k-1)*Nr+1:k*Nr, :) = [];

    [~, ~, V] = svd(G_tilde);
    rv = rank(G_tilde);

    Htmp = Gk * V(:, rv+1:Nt);
    [~, ~, Vt] = svd(Htmp);
    Ftmp = V(:, rv+1:Nt) * Vt(:, 1:Ns_per_user(k));
    F_BD = [F_BD, Ftmp];
end

end
