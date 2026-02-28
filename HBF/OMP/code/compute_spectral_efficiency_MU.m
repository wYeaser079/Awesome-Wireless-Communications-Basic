function [R_sum, R_per_user] = compute_spectral_efficiency_MU(H_cell, W_total, P_per_user, Ns_per_user, sigma2)
% COMPUTE_SPECTRAL_EFFICIENCY_MU  Compute the achievable sum spectral
%   efficiency for a multi-user MIMO downlink system.
%
%   [R_sum, R_per_user] = compute_spectral_efficiency_MU(H_cell, W_total, P_per_user, Ns_per_user, sigma2)
%
%   Inputs:
%       H_cell       - Cell array of K channel matrices: H_cell{k} is Nr_k x Nt
%       W_total      - Nt x sum(Ns_per_user) total precoding matrix
%       P_per_user   - K x 1 vector of transmit powers per user
%       Ns_per_user  - K x 1 vector of data streams per user
%       sigma2       - Noise variance
%
%   Outputs:
%       R_sum        - Sum spectral efficiency (bits/s/Hz)
%       R_per_user   - K x 1 vector of per-user spectral efficiencies
%
%   For single-antenna users with ZF precoding on the effective channel,
%   inter-user interference is eliminated and each user sees:
%       R_k = log2(1 + P_k * ||h_k * w_k||^2 / sigma2)
%
%   For multi-antenna users, the SINR accounting considers interference:
%       R_k = log2|I + P_k * H_k * W_k * W_k^H * H_k^H * (sigma2*I + sum_{j!=k} P_j * H_k * W_j * W_j^H * H_k^H)^{-1}|
%
%   Reference:
%       Le Liang, HybridPrecodingMassiveMIMO (CalRate.m)
%
% =========================================================================

K = length(H_cell);
R_per_user = zeros(K, 1);

% Build stream index mapping
stream_idx = cell(K, 1);
offset = 0;
for k = 1:K
    stream_idx{k} = offset + (1:Ns_per_user(k));
    offset = offset + Ns_per_user(k);
end

for k = 1:K
    Hk = H_cell{k};                        % Nr_k x Nt
    Wk = W_total(:, stream_idx{k});         % Nt x Ns_k
    Nr_k = size(Hk, 1);

    % Signal power for user k
    P_sig = P_per_user(k) * Hk * Wk * Wk' * Hk';

    % Interference from all other users
    P_int = zeros(Nr_k);
    for j = 1:K
        if j ~= k
            Wj = W_total(:, stream_idx{j});
            P_int = P_int + P_per_user(j) * Hk * Wj * Wj' * Hk';
        end
    end

    % Rate for user k
    R_per_user(k) = real(log2(det(eye(Nr_k) + (sigma2 * eye(Nr_k) + P_int) \ P_sig)));
end

R_sum = sum(R_per_user);

end
