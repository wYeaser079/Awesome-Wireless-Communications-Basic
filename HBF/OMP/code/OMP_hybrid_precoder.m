function [F_RF, F_BB] = OMP_hybrid_precoder(F_opt, At, Ns, N_RF)
% OMP_HYBRID_PRECODER  OMP-based hybrid precoder design (El Ayach et al. 2014).
%
%   [F_RF, F_BB] = OMP_hybrid_precoder(F_opt, At, Ns, N_RF)
%
%   Inputs:
%       F_opt  - Nt x Ns optimal (unconstrained) digital precoder
%                (typically the first Ns columns of V from SVD of H)
%       At     - Nt x G transmit dictionary matrix (columns are candidate
%                steering vectors for the RF precoder)
%       Ns     - Number of data streams
%       N_RF   - Number of RF chains at the transmitter
%
%   Outputs:
%       F_RF   - Nt x N_RF analog/RF precoding matrix
%                (each column is a steering vector from At)
%       F_BB   - N_RF x Ns digital/baseband precoding matrix
%                (normalized so that ||F_RF * F_BB||_F^2 = Ns)
%
%   Algorithm (Algorithm 1 in El Ayach et al. 2014):
%   -------------------------------------------------------
%   1. Initialize: F_res = F_opt, F_RF = []
%   2. For each RF chain r = 1, ..., N_RF:
%      a. Compute Psi = At^H * F_res
%      b. Find k = argmax_i (Psi * Psi^H)_{i,i}   (best matching atom)
%      c. Append: F_RF = [F_RF, At(:,k)]
%      d. Compute: F_BB = (F_RF^H * F_RF)^{-1} * F_RF^H * F_opt   (LS)
%      e. Update residual: F_res = (F_opt - F_RF * F_BB) / ||...||_F
%   3. Normalize: F_BB = sqrt(Ns) * F_BB / ||F_RF * F_BB||_F
%
%   References:
%       [1] O. El Ayach, S. Rajagopal, S. Abu-Surra, Z. Pi, and R. W. Heath,
%           "Spatially Sparse Precoding in Millimeter Wave MIMO Systems,"
%           IEEE Trans. Wireless Commun., vol. 13, no. 3, pp. 1499-1513, 2014.
%       [2] ge99210/Hybrid-Precoding-Combining- (hp_omp.m)
%       [3] Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%       [4] metegenez/spatially-sparse-precoding
%
% =========================================================================

% Validate inputs
[Nt, ~] = size(At);
assert(size(F_opt, 1) == Nt, 'F_opt and At must have the same number of rows (Nt)');
assert(N_RF >= Ns, 'Number of RF chains must be >= number of data streams');

F_RF  = [];
F_res = F_opt;

for r = 1:N_RF
    % Step 2a: Project residual onto dictionary
    Psi = At' * F_res;                     % G x Ns

    % Step 2b: Select the dictionary atom with largest projection energy
    [~, k] = max(diag(Psi * Psi'));        % scalar index

    % Step 2c: Augment the RF precoder
    F_RF = [F_RF, At(:, k)];              % Nt x r

    % Step 2d: Least-squares baseband precoder
    F_BB = (F_RF' * F_RF) \ (F_RF' * F_opt);   % r x Ns

    % Step 2e: Update residual
    F_res = F_opt - F_RF * F_BB;
    F_res = F_res / norm(F_res, 'fro');
end

% Step 3: Power normalization so that ||F_RF * F_BB||_F^2 = Ns
F_BB = sqrt(Ns) * F_BB / norm(F_RF * F_BB, 'fro');

end
