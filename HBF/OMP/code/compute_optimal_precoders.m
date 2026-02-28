function [Fopt, Wopt] = compute_optimal_precoders(H, Ns, SNR_lin)
% COMPUTE_OPTIMAL_PRECODERS  Compute the optimal unconstrained digital
%   precoder and MMSE combiner for a MIMO channel.
%
%   [Fopt, Wopt] = compute_optimal_precoders(H, Ns, SNR_lin)
%
%   Inputs:
%       H       - Nr x Nt channel matrix
%       Ns      - Number of data streams
%       SNR_lin - SNR in linear scale (rho = P / sigma^2)
%
%   Outputs:
%       Fopt    - Nt x Ns optimal digital precoder
%                 (the Ns dominant right singular vectors of H)
%       Wopt    - Nr x Ns optimal MMSE combiner
%
%   The optimal precoder is obtained from the SVD of H:
%       [U, S, V] = svd(H)
%       Fopt = V(:, 1:Ns)
%
%   This maximizes the channel capacity by directing energy along the
%   strongest singular modes of the channel.
%
%   The MMSE combiner is:
%       Wopt = ((1/rho)*Fopt'*H'*H*Fopt + (Ns/rho)*I_Ns)^{-1} * Fopt'*H'
%       Wopt = (1/sqrt(rho)) * Wopt.'
%
%   where the transpose converts from the row-form solution to the
%   standard Nr x Ns combining matrix format.
%
%   Reference:
%       O. El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., 2014.

    % ======================================================================
    % SVD of channel matrix
    % ======================================================================
    [U, S, V] = svd(H);

    % ======================================================================
    % Optimal precoder: Ns dominant right singular vectors
    % ======================================================================
    Fopt = V(:, 1:Ns);

    % ======================================================================
    % Optimal MMSE combiner
    % ======================================================================
    % The MMSE combiner minimizes E[||s - W'*y||^2]
    % where y = H*F*s + n
    %
    % Solution:
    %   W_MMSE = (H*F*F'*H' + (Ns/rho)*I)^{-1} * H*F
    % Or equivalently (using matrix inversion lemma):
    %   W_MMSE = (1/sqrt(rho)) * ((1/rho)*F'*H'*H*F + (Ns/rho)*I)^{-1} * F'*H'
    %   then transpose.

    HFopt = H * Fopt;   % Nr x Ns

    Wopt = (1/sqrt(SNR_lin)) * ...
        ((1/SNR_lin) * (Fopt' * H' * H * Fopt) + ...
         (Ns / SNR_lin) * eye(Ns)) \ ...
        (Fopt' * H');
    Wopt = Wopt.';  % Convert to Nr x Ns

end
