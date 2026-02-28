function R = compute_spectral_efficiency(H, F, W, Ns, SNR)
% COMPUTE_SPECTRAL_EFFICIENCY  Compute the achievable spectral efficiency
%   (bits/s/Hz) for a single-user MIMO system with given precoder and combiner.
%
%   R = compute_spectral_efficiency(H, F, W, Ns, SNR)
%
%   Inputs:
%       H    - Nr x Nt channel matrix
%       F    - Nt x Ns precoding matrix (can be F_RF*F_BB for hybrid)
%       W    - Nr x Ns combining matrix (can be W_RF*W_BB for hybrid)
%       Ns   - Number of data streams
%       SNR  - Signal-to-noise ratio (linear scale)
%
%   Output:
%       R    - Spectral efficiency in bits/s/Hz
%
%   Formula (from El Ayach et al. 2014, Eq. 6):
%
%       R = log2 |I_Ns + (SNR/Ns) * R_n^{-1} * W^H * H * F * F^H * H^H * W|
%
%   where:
%       R_n = W^H * W  (noise covariance after combining)
%
%   For the optimal unconstrained case (F = V(:,1:Ns), W = U(:,1:Ns)):
%       R = sum_{i=1}^{Ns} log2(1 + (SNR/Ns) * sigma_i^2)
%
%   References:
%       [1] El Ayach et al., IEEE TWC, 2014.
%       [2] Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%       [3] meuseabe/deepHybridBeamforming (helperComputeSpectralEfficiency.m)
%
% =========================================================================

% Noise covariance after combining (noise coloring correction)
% For hybrid combining where W = W_RF * W_BB, the noise after combining
% is n_eff = W' * n, so its covariance is E[n_eff * n_eff'] = sigma^2 * W' * W.
% The term Rn = W' * W accounts for the fact that the analog combiner W_RF
% is not unitary (its columns are unit-norm steering vectors but NOT orthogonal),
% so the noise is "colored" after combining. For the optimal unconstrained case
% where W = U(:,1:Ns), W'*W = I and this reduces to the standard formula.
Rn = W' * W;

% Signal term after combining
HF = H * F;       % Nr x Ns
WHF = W' * HF;    % Ns x Ns

% Spectral efficiency
R = real(log2(det(eye(Ns) + (SNR / Ns) * (Rn \ (WHF * WHF')))));

% Take absolute value to handle small numerical issues
R = abs(R);

end
