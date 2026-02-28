function [W_RF, W_BB] = OMP_hybrid_combiner(H, F_RF, F_BB, Ar, Ns, N_RF, SNR)
% OMP_HYBRID_COMBINER  OMP-based hybrid combiner design (El Ayach et al. 2014).
%
%   [W_RF, W_BB] = OMP_hybrid_combiner(H, F_RF, F_BB, Ar, Ns, N_RF, SNR)
%
%   Inputs:
%       H      - Nr x Nt channel matrix
%       F_RF   - Nt x N_RF  analog/RF precoding matrix (from OMP precoder)
%       F_BB   - N_RF x Ns  digital/baseband precoding matrix
%       Ar     - Nr x G     receive dictionary matrix
%       Ns     - Number of data streams
%       N_RF   - Number of RF chains at the receiver
%       SNR    - Signal-to-noise ratio (linear scale)
%
%   Outputs:
%       W_RF   - Nr x N_RF  analog/RF combining matrix
%       W_BB   - N_RF x Ns  digital/baseband combining matrix
%
%   Algorithm:
%   -------------------------------------------------------
%   First compute the MMSE optimal combiner:
%       W_MMSE = (1/sqrt(SNR)) * (F_BB^H * F_RF^H * H^H * H * F_RF * F_BB
%                + (Ns/SNR) * I)^{-1} * F_BB^H * F_RF^H * H^H
%
%   Then compute the receive covariance:
%       C_yy = (SNR/Ns) * H * F_RF * F_BB * F_BB^H * F_RF^H * H^H + I
%
%   Then use OMP to decompose W_MMSE using Ar as dictionary:
%   1. W_res = W_MMSE, W_RF = []
%   2. For r = 1, ..., N_RF:
%      a. Psi = Ar^H * C_yy * W_res
%      b. k = argmax_i (Psi * Psi^H)_{i,i}
%      c. W_RF = [W_RF, Ar(:,k)]
%      d. W_BB = (W_RF^H * C_yy * W_RF)^{-1} * W_RF^H * C_yy * W_MMSE
%      e. W_res = (W_MMSE - W_RF * W_BB) / ||...||_F
%
%   References:
%       [1] El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%           MIMO Systems," IEEE TWC, 2014.
%       [2] meuseabe/deepHybridBeamforming (helperOMPHybridWeights.m)
%       [3] Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%
% =========================================================================

Nr = size(H, 1);

% --- Compute MMSE optimal (unconstrained) combiner ---
F_total = F_RF * F_BB;   % Nt x Ns effective precoder
W_MMSE = ((F_BB' * F_RF' * (H' * H) * F_RF * F_BB + ...
           (Ns / SNR) * eye(Ns)) \ (F_BB' * F_RF' * H'))';
% W_MMSE is Nr x Ns

% --- Receive covariance matrix ---
Ess = (1/Ns) * eye(Ns);  % signal covariance per stream
C_yy = H * F_RF * F_BB * Ess * F_BB' * F_RF' * H' + (1/SNR) * eye(Nr);
% Alternatively: C_yy = (SNR/Ns) * H * F_total * F_total' * H' + eye(Nr);

% --- OMP for hybrid combiner ---
W_RF  = [];
W_res = W_MMSE;

for r = 1:N_RF
    % Project weighted residual onto dictionary
    Psi = Ar' * C_yy * W_res;              % G x Ns

    % Select best atom
    [~, k] = max(diag(Psi * Psi'));

    % Augment RF combiner
    W_RF = [W_RF, Ar(:, k)];               % Nr x r

    % Weighted least-squares baseband combiner
    W_BB = (W_RF' * C_yy * W_RF) \ (W_RF' * C_yy * W_MMSE);  % r x Ns

    % Update residual
    W_res = W_MMSE - W_RF * W_BB;
    W_res = W_res / norm(W_res, 'fro');
end

end
