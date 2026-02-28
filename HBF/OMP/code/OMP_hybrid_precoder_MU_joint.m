function [FRF, FBB, Heff] = OMP_hybrid_precoder_MU_joint(H_all, At_dict, Nt, NtRF, Ns_per_user, K)
% OMP_HYBRID_PRECODER_MU_JOINT  Multi-user hybrid precoding using JOINT OMP
%   for the analog RF stage and Zero-Forcing (ZF) for the digital baseband.
%
%   [FRF, FBB, Heff] = OMP_hybrid_precoder_MU_joint(H_all, At_dict, Nt, NtRF, Ns_per_user, K)
%
%   This implements APPROACH B (Joint OMP): Run OMP on the aggregate
%   stacked channel to find a shared FRF across all users.
%
%   For APPROACH A (Per-user OMP), see OMP_hybrid_precoder_MU.m
%
%   Inputs:
%       H_all       - Cell array of K channel matrices, H_all{k} is Nr x Nt
%                     (or equivalently: a (K*Nr) x Nt stacked channel matrix)
%       At_dict     - Nt x Ndict dictionary of Tx array response vectors
%       Nt          - Number of transmit antennas
%       NtRF        - Number of transmit RF chains
%       Ns_per_user - Number of data streams per user
%       K           - Number of users
%
%   Outputs:
%       FRF         - Nt x NtRF analog (RF) precoding matrix
%       FBB         - NtRF x (K*Ns_per_user) digital baseband precoding matrix
%       Heff        - (K*Ns_per_user) x NtRF effective channel after RF precoding
%                     (used for debugging / verification)
%
%   Two-Stage Multi-User Hybrid Precoding (Joint OMP):
%   ===================================================
%
%   STAGE 1: Analog RF Precoder (Joint OMP on aggregate channel)
%   ------------------------------------------------------------
%   1. Stack all user channels: H_stack = [H_1; H_2; ...; H_K]
%   2. Compute SVD of H_stack to get Fopt_agg = V(:, 1:NtRF)
%   3. Run OMP on Fopt_agg to select NtRF beams from the dictionary
%      that best approximate the aggregate optimal precoder.
%
%   STAGE 2: Digital Baseband Precoder (Zero-Forcing)
%   -------------------------------------------------
%   1. Compute effective channel: H_eff_k = H_k * FRF   for each user k
%   2. Stack: Heff = [H_eff_1; H_eff_2; ...; H_eff_K]
%   3. ZF precoder on effective channel:
%       FBB = Heff' * inv(Heff * Heff')
%   4. Normalize columns to satisfy power constraint.
%
%   References:
%       [1] O. El Ayach et al., "Spatially Sparse Precoding in Millimeter
%           Wave MIMO Systems," IEEE Trans. Wireless Commun., 2014.
%       [2] A. Alkhateeb et al., "Limited Feedback Hybrid Precoding for
%           Multi-Layer Millimeter Wave MIMO," IEEE JSAC, 2015.
%       [3] X. Yu et al., "Alternating Minimization Algorithms for Hybrid
%           Precoding in Millimeter Wave MIMO Systems," IEEE JSTSP, 2016.

    Ns_total = K * Ns_per_user;

    % ======================================================================
    % STAGE 1: Determine the analog RF precoder FRF via Joint OMP
    % ======================================================================
    % Strategy: Stack all user channels into one big channel, compute SVD,
    % and run OMP to find the best NtRF atoms from the dictionary.
    %
    % Alternative: Run per-user OMP with NtRF/K beams each (see
    % OMP_hybrid_precoder_MU.m for that approach).

    % Stack channels: H_stack is (K*Nr) x Nt
    if iscell(H_all)
        H_stack = cell2mat(H_all(:));
    else
        H_stack = H_all;  % Already stacked
    end

    % SVD of stacked channel to get optimal aggregate precoder
    [~, ~, V_stack] = svd(H_stack);
    Fopt_agg = V_stack(:, 1:NtRF);  % Nt x NtRF

    % Run OMP on the aggregate optimal precoder
    % (Selecting NtRF beams from the dictionary)
    FRF = [];
    Fres = Fopt_agg;

    for i = 1:NtRF
        % Correlation
        Psi = At_dict' * Fres;

        % Atom selection
        [~, k_best] = max(sum(abs(Psi).^2, 2));

        % Append atom
        FRF = [FRF, At_dict(:, k_best)]; %#ok<AGROW>

        % Least squares
        FBB_temp = (FRF' * FRF) \ (FRF' * Fopt_agg);

        % Residual update
        Fres_unnorm = Fopt_agg - FRF * FBB_temp;
        nrm = norm(Fres_unnorm, 'fro');
        if nrm > 1e-10
            Fres = Fres_unnorm / nrm;
        end
    end

    % ======================================================================
    % STAGE 2: Digital baseband ZF precoder on the effective channel
    % ======================================================================
    % Compute effective channel for each user
    if iscell(H_all)
        Heff = zeros(0, NtRF);
        for kk = 1:K
            Heff = [Heff; H_all{kk} * FRF]; %#ok<AGROW>
        end
    else
        Nr_per_user = size(H_stack, 1) / K;
        Heff = zeros(K * Nr_per_user, NtRF);
        for kk = 1:K
            row_idx = (kk-1)*Nr_per_user + (1:Nr_per_user);
            Heff(row_idx, :) = H_stack(row_idx, :) * FRF;
        end
    end

    % For ZF with single-antenna users (or Ns_per_user = 1):
    % If each user has Nr antennas but only Ns_per_user=1 stream,
    % we use the "strongest direction" for each user.
    % If Ns_per_user = Nr per user, use full Heff.

    % Assume single-antenna users or use only Ns_per_user rows per user
    if iscell(H_all)
        Nr_per_user = size(H_all{1}, 1);
    else
        Nr_per_user = size(H_stack, 1) / K;
    end

    if Ns_per_user == 1 && Nr_per_user > 1
        % Each user uses the strongest direction in the effective channel
        Heff_reduced = zeros(K, NtRF);
        for kk = 1:K
            row_idx = (kk-1)*Nr_per_user + (1:Nr_per_user);
            Heff_k = Heff(row_idx, :);
            [U_k, ~, ~] = svd(Heff_k);
            Heff_reduced(kk, :) = U_k(:, 1)' * Heff_k;
        end
        Heff_for_ZF = Heff_reduced;  % K x NtRF
    else
        Heff_for_ZF = Heff;  % (K*Ns_per_user) x NtRF (when Nr=Ns_per_user=1)
    end

    % Zero-Forcing baseband precoder:
    %   FBB = Heff' * inv(Heff * Heff')
    FBB = Heff_for_ZF' * inv(Heff_for_ZF * Heff_for_ZF');  % NtRF x K*Ns_per_user

    % Power normalization: ||FRF * FBB||_F^2 = Ns_total
    FBB = sqrt(Ns_total) * FBB / norm(FRF * FBB, 'fro');

    % Return the effective channel for analysis
    Heff = Heff_for_ZF;

end
