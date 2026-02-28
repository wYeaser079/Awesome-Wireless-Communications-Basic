function [F_bs, W_bs] = beam_steering(H, At, Ar, Ns)
% BEAM_STEERING  Select the best beam-steering precoder and combiner.
%
%   [F_bs, W_bs] = beam_steering(H, At, Ar, Ns)
%
%   Inputs:
%       H   - Nr x Nt channel matrix
%       At  - Nt x Npaths matrix of Tx array response vectors
%       Ar  - Nr x Npaths matrix of Rx array response vectors
%       Ns  - Number of data streams (typically 1 for beam steering)
%
%   Outputs:
%       F_bs - Nt x Ns beam steering precoder
%       W_bs - Nr x Ns beam steering combiner
%
%   For beam steering with Ns = 1, we select the Tx and Rx array response
%   vectors that maximize the effective channel gain |ar' * H * at|.
%
%   For Ns > 1, we greedily select the Ns best beam pairs.
%
%   Reference:
%       O. El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., 2014.

    Npaths = size(At, 2);

    if Ns == 1
        % Exhaustive search over all Tx-Rx beam pairs
        best_gain = 0;
        best_t = 1;
        best_r = 1;

        for t = 1:Npaths
            for r = 1:Npaths
                gain = abs(Ar(:, r)' * H * At(:, t));
                if gain > best_gain
                    best_gain = gain;
                    best_t = t;
                    best_r = r;
                end
            end
        end

        F_bs = At(:, best_t);
        W_bs = Ar(:, best_r);

    else
        % Greedy selection for Ns > 1
        F_bs = zeros(size(At, 1), Ns);
        W_bs = zeros(size(Ar, 1), Ns);
        used_t = [];
        used_r = [];

        for s = 1:Ns
            best_gain = 0;
            best_t = 1;
            best_r = 1;

            for t = 1:Npaths
                if ismember(t, used_t), continue; end
                for r = 1:Npaths
                    if ismember(r, used_r), continue; end
                    gain = abs(Ar(:, r)' * H * At(:, t));
                    if gain > best_gain
                        best_gain = gain;
                        best_t = t;
                        best_r = r;
                    end
                end
            end

            F_bs(:, s) = At(:, best_t);
            W_bs(:, s) = Ar(:, best_r);
            used_t = [used_t, best_t]; %#ok<AGROW>
            used_r = [used_r, best_r]; %#ok<AGROW>
        end
    end
end
