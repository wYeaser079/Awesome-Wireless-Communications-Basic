# Comprehensive Guide: Optimization Problem Formulation in Wireless Communications
## Multi-User Massive MIMO-OFDM Systems with Hybrid Beamforming

---

## Table of Contents
1. [Understanding the Paper's Formulation](#1-understanding-the-papers-formulation)
2. [Four Major Optimization Objectives](#2-four-major-optimization-objectives)
3. [Mathematical Formulations](#3-mathematical-formulations)
4. [Solution Approaches](#4-solution-approaches)
5. [MATLAB Implementation](#5-matlab-implementation)
6. [Python Implementation](#6-python-implementation)
7. [Key References](#7-key-references)

---

## 1. Understanding the Paper's Formulation

### 1.1 System Model Overview

The paper considers a **multi-user MIMO-OFDM THz system** with:
- **Base Station (BS)**: Equipped with `Nt` transmit antennas and `NRF` RF chains
- **Users**: `U` single-antenna users
- **OFDM**: `K` sub-carriers
- **Hybrid Beamforming**: Analog precoder `W` (Nt × NRF) + Digital precoder `F` (NRF × U × K)

### 1.2 Signal Model

The received signal at the k-th sub-carrier for the q-th user:

```
y[k] = S₀ Hq[k] W Fq[k] sq[k] + Δq,k + n[k]
```

Where:
- `S₀`: ICI coefficient for the desired signal
- `Hq[k]`: Channel vector for user q at subcarrier k (1 × Nt)
- `W`: Analog beamformer (Nt × NRF)
- `Fq[k]`: Digital beamformer for user q at subcarrier k
- `sq[k]`: Transmitted symbol
- `Δq,k`: Total interference (inter-user + ICI)
- `n[k]`: AWGN noise

### 1.3 SINR Expression

```
            |S₀|² |Hq[k] W Fq[k]|²
SINRq,k = ─────────────────────────
            |Δq,k|² + ψ
```

Where `ψ = Kσ²n/P` (noise power normalized)

### 1.4 Spectral Efficiency

Per-user per-subcarrier: `γq,k = log₂(1 + SINRq,k)`

Total system: `SE = Σq Σk γq,k`

---

## 2. Four Major Optimization Objectives

### Objective 1: Maximize Total Sum Rate
**Goal**: Maximize the aggregate data rate across all users

### Objective 2: Minimize Total Power Consumption
**Goal**: Achieve target QoS (SINR) with minimum transmit power

### Objective 3: Max-Min Fairness (Maximize Minimum Rate)
**Goal**: Ensure fairness by maximizing the worst user's rate

### Objective 4: Maximize Per-User Spectral Efficiency
**Goal**: Optimize individual user performance (can be weighted)

---

## 3. Mathematical Formulations

### 3.1 PROBLEM 1: Maximize Total Sum Rate

**Objective**: Maximize the sum of achievable rates across all users and subcarriers.

```
┌─────────────────────────────────────────────────────────────────┐
│  MAXIMIZE:                                                       │
│                U    K                                            │
│     R_sum =   Σ    Σ   log₂(1 + SINRq,k)                        │
│              q=1  k=1                                            │
│                                                                  │
│  SUBJECT TO:                                                     │
│                                                                  │
│  C1: |(W)m,n| = 1/√Nt,  ∀m ∈ [1,Nt], ∀n ∈ [1,NRF]              │
│      (Constant modulus constraint for analog beamformer)         │
│                                                                  │
│  C2: ||W Fq[k]||²F = 1,  ∀q ∈ [1,U], ∀k ∈ [1,K]                │
│      (Power normalization per user per subcarrier)               │
│                                                                  │
│  C3: Σq Σk ||W Fq[k]||²F ≤ P_total                              │
│      (Total power constraint)                                    │
│                                                                  │
│  VARIABLES: W ∈ C^(Nt×NRF), F[k] ∈ C^(NRF×U) ∀k                 │
└─────────────────────────────────────────────────────────────────┘
```

**Extended SINR with ICI:**
```
                |S₀|² |Hq[k] W Fq[k]|²
SINRq,k = ───────────────────────────────────────────────────────
          Σ(u≠q) Σi |Si-k|² |Hq[i] W Fu[i]|² + Σ(i≠k) |Si-k|² |Hq[i] W Fq[i]|² + σ²
```

---

### 3.2 PROBLEM 2: Minimize Total Power Consumption

**Objective**: Minimize transmit power while guaranteeing minimum SINR for each user.

```
┌─────────────────────────────────────────────────────────────────┐
│  MINIMIZE:                                                       │
│                U    K                                            │
│     P_total = Σ    Σ   ||W Fq[k]||²F                            │
│              q=1  k=1                                            │
│                                                                  │
│  SUBJECT TO:                                                     │
│                                                                  │
│  C1: SINRq,k ≥ γ_min,  ∀q ∈ [1,U], ∀k ∈ [1,K]                  │
│      (Minimum SINR requirement per user per subcarrier)          │
│                                                                  │
│  C2: |(W)m,n| = 1/√Nt,  ∀m, ∀n                                  │
│      (Constant modulus constraint)                               │
│                                                                  │
│  C3: ||W Fq[k]||²F ≤ P_max_per_user                             │
│      (Optional: per-user power constraint)                       │
│                                                                  │
│  VARIABLES: W, F[k]                                              │
└─────────────────────────────────────────────────────────────────┘
```

**Alternative Formulation (with per-antenna constraints):**
```
MINIMIZE:  Σm [W F F^H W^H]m,m    (sum of per-antenna powers)

SUBJECT TO:
  C1: SINRq,k ≥ γq,min,  ∀q,k
  C2: [W F F^H W^H]m,m ≤ P_ant_max,  ∀m  (per-antenna power limit)
  C3: |(W)m,n| = 1/√Nt
```

---

### 3.3 PROBLEM 3: Max-Min Fairness (Maximize Minimum Rate)

**Objective**: Ensure user fairness by maximizing the minimum achievable rate.

```
┌─────────────────────────────────────────────────────────────────┐
│  MAXIMIZE:                                                       │
│                                                                  │
│     min   { Σk log₂(1 + SINRq,k) }                              │
│    q∈[1,U]                                                       │
│                                                                  │
│  SUBJECT TO:                                                     │
│                                                                  │
│  C1: |(W)m,n| = 1/√Nt,  ∀m, ∀n                                  │
│                                                                  │
│  C2: Σq Σk ||W Fq[k]||²F ≤ P_total                              │
│                                                                  │
│  C3: ||W Fq[k]||²F ≤ P_max_per_user,  ∀q,k                      │
│                                                                  │
│  VARIABLES: W, F[k]                                              │
└─────────────────────────────────────────────────────────────────┘
```

**Epigraph Form (Standard Transformation):**
```
MAXIMIZE:  t

SUBJECT TO:
  C1: Σk log₂(1 + SINRq,k) ≥ t,  ∀q ∈ [1,U]
  C2: |(W)m,n| = 1/√Nt
  C3: Σq Σk ||W Fq[k]||²F ≤ P_total

VARIABLES: W, F[k], t
```

**Weighted Max-Min (Prioritized Fairness):**
```
MAXIMIZE:  min { wq · Σk log₂(1 + SINRq,k) }
          q∈[1,U]

Where wq is the priority weight for user q
```

---

### 3.4 PROBLEM 4: Maximize Per-User Weighted Sum Rate

**Objective**: Maximize weighted combination of individual user rates.

```
┌─────────────────────────────────────────────────────────────────┐
│  MAXIMIZE:                                                       │
│           U                                                      │
│          Σ   αq · Rq                                            │
│         q=1                                                      │
│                                                                  │
│  Where Rq = Σk log₂(1 + SINRq,k) is user q's total rate         │
│        αq is the priority weight for user q                      │
│                                                                  │
│  SUBJECT TO:                                                     │
│                                                                  │
│  C1: Rq ≥ Rq,min,  ∀q  (minimum rate guarantee per user)        │
│                                                                  │
│  C2: |(W)m,n| = 1/√Nt,  ∀m, ∀n                                  │
│                                                                  │
│  C3: Σq Σk ||W Fq[k]||²F ≤ P_total                              │
│                                                                  │
│  VARIABLES: W, F[k]                                              │
└─────────────────────────────────────────────────────────────────┘
```

---

## 4. Solution Approaches

### 4.1 Key Challenges

All four problems are **non-convex** due to:
1. **Constant modulus constraint**: `|(W)m,n| = 1/√Nt` (unit circle constraint)
2. **Coupling between W and F**: Objective and constraints involve products
3. **Log function with interference**: Creates non-convex SINR expressions

### 4.2 Solution Framework: Alternating Optimization

```
┌────────────────────────────────────────────────────────────────────┐
│                   GENERAL SOLUTION FRAMEWORK                        │
├────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  STEP 1: Initialize W randomly (satisfying constant modulus)       │
│                                                                     │
│  STEP 2: Fix W, solve for F (Digital Beamformer)                   │
│          → Use Zero-Forcing (ZF) or MMSE                           │
│                                                                     │
│  STEP 3: Fix F, solve for W (Analog Beamformer)                    │
│          → Use Manifold Optimization                                │
│                                                                     │
│  STEP 4: Repeat Steps 2-3 until convergence                        │
│                                                                     │
└────────────────────────────────────────────────────────────────────┘
```

### 4.3 Digital Beamformer Design (Zero-Forcing)

**Effective Channel:**
```
He[k]^H = H[k]^H W
```

**ZF Precoder:**
```
F̃u[k] = (He[k] He[k]^H)^(-1) He[k]

# Normalized:
Fu[k] = F̃u[k] / ||W F̃u[k]||F
```

**Why ZF Works:**
- Nullifies inter-user interference: `Hq[k] W Fu[k] = 0` for u ≠ q
- Simplifies SINR by removing cross-user interference terms

### 4.4 Analog Beamformer Design (Manifold Optimization)

**Problem Reformulation:**
Due to constant modulus constraint, the problem lies on a **complex circle manifold**:
```
M = {x ∈ C^(Nt·NRF) : |xi| = 1/√Nt, ∀i}
```

**Riemannian Conjugate Gradient Algorithm:**

```
Algorithm: Manifold Optimization for Analog Beamformer
────────────────────────────────────────────────────────
Input: F[k], initial x₀, step size α, threshold ε₁
Output: Optimized W

1: d₀ = Tx₀(∇f(x)|x=x₀), t = 0
2: while |f(xt) - f(xt-1)| > ε₁ do:
3:    xt+1 = xt + α·dt           # Update position
4:    xt+1 = R(xt+1)             # Retraction onto manifold
5:    gt+1 = Txt+1(∇f(x)|x=xt+1) # Riemannian gradient
6:    βt+1 = (gt+1^H(gt+1 - Txt+1(gt))) / ||Txt+1(gt)||²F  # Polak-Ribiere
7:    dt+1 = gt+1 + βt+1·Txt+1(dt)  # Conjugate direction
8:    t ← t + 1
9: return vec⁻¹(xt)  # Reshape to matrix W
────────────────────────────────────────────────────────

Retraction function:
R(x) = (1/√Nt) · vec[x₁/|x₁|, x₂/|x₂|, ..., xN/|xN|]
```

### 4.5 WMMSE Approach (Alternative for Sum Rate Maximization)

The **Weighted Minimum Mean Square Error (WMMSE)** approach converts the non-convex sum rate maximization into an equivalent problem:

**Key Relationship:**
```
max Σq log₂(1 + SINRq)  ⟺  min Σq (wq·eq - log(wq))
```

Where:
- `eq`: MSE for user q
- `wq`: Weight for user q (optimally wq = eq^(-1))

**WMMSE Algorithm:**
```
1: Initialize beamformers W, F
2: repeat:
3:    Update receive filters: uq = ...
4:    Update MSE weights: wq = eq^(-1)
5:    Update transmit beamformers: solve weighted MMSE problem
6: until convergence
```

### 4.6 Semidefinite Relaxation (SDR) Approach

For power minimization problems, SDR transforms the non-convex problem:

**Original (Non-convex):**
```
min ||wq||²  s.t. SINRq ≥ γq
```

**SDR Form:**
Define `Wq = wq wq^H` (rank-1 matrix)
```
min Tr(Wq)
s.t. Tr(Hq Wq Hq^H) / (Σu≠q Tr(Hq Wu Hq^H) + σ²) ≥ γq
     Wq ⪰ 0  (positive semidefinite)
     rank(Wq) = 1  (relaxed!)
```

**Solve convex SDP, then extract rank-1 solution via eigenvalue decomposition.**

---

## 5. MATLAB Implementation

### 5.1 Problem 1: Sum Rate Maximization

```matlab
%% Sum Rate Maximization with Hybrid Beamforming
% System Parameters
Nt = 64;          % Transmit antennas
NRF = 9;          % RF chains
U = 16;           % Number of users
K = 64;           % Subcarriers
P_total = 1;      % Total power (normalized)
sigma2 = 0.01;    % Noise variance

% Channel Generation (Saleh-Valenzuela Model)
function H = generate_channel(Nt, U, K, Nc, Nray)
    H = zeros(U, Nt, K);
    for q = 1:U
        for k = 1:K
            for c = 1:Nc
                for r = 1:Nray
                    phi = 2*pi*rand();
                    theta = pi*rand() - pi/2;
                    alpha = (randn() + 1j*randn())/sqrt(2);
                    a = exp(1j*2*pi*(0:Nt-1)'*0.5*sin(phi)*sin(theta))/sqrt(Nt);
                    H(q,:,k) = H(q,:,k) + sqrt(Nt/(Nc*Nray))*alpha*a.';
                end
            end
        end
    end
end

% Zero-Forcing Digital Precoder
function F = compute_ZF_precoder(H, W, K, U)
    F = zeros(size(W,2), U, K);
    for k = 1:K
        Hk = squeeze(H(:,:,k));  % U x Nt
        He = Hk * W;              % Effective channel: U x NRF
        F_zf = He' / (He * He' + 1e-6*eye(U));  % ZF precoder
        % Normalize
        for u = 1:U
            F(:,u,k) = F_zf(:,u) / norm(W * F_zf(:,u));
        end
    end
end

% Sum Rate Calculation
function R = compute_sum_rate(H, W, F, U, K, sigma2, S)
    R = 0;
    for q = 1:U
        for k = 1:K
            Hqk = squeeze(H(q,:,k));
            signal = abs(S(1) * Hqk * W * F(:,q,k))^2;

            % Interference calculation
            interference = 0;
            for u = 1:U
                for i = 1:K
                    if ~(u==q && i==k)
                        S_coef = S(abs(i-k)+1);
                        interference = interference + ...
                            abs(S_coef)^2 * abs(squeeze(H(q,:,i)) * W * F(:,u,i))^2;
                    end
                end
            end

            SINR = signal / (interference + sigma2);
            R = R + log2(1 + SINR);
        end
    end
end

% Manifold Optimization for Analog Beamformer
function W = optimize_analog_manifold(W_init, H, F, U, K, sigma2, S, max_iter, alpha)
    W = W_init;
    Nt = size(W, 1);
    NRF = size(W, 2);

    for iter = 1:max_iter
        % Compute gradient (numerical)
        grad = zeros(size(W));
        epsilon = 1e-6;

        for m = 1:Nt
            for n = 1:NRF
                W_plus = W;
                W_plus(m,n) = W_plus(m,n) + epsilon;
                W_minus = W;
                W_minus(m,n) = W_minus(m,n) - epsilon;

                % Note: We minimize interference, which is equivalent to maximizing rate
                I_plus = compute_total_interference(H, W_plus, F, U, K, S);
                I_minus = compute_total_interference(H, W_minus, F, U, K, S);
                grad(m,n) = (I_plus - I_minus) / (2*epsilon);
            end
        end

        % Riemannian gradient (project onto tangent space)
        riem_grad = grad - real(grad .* conj(W)) .* W;

        % Update
        W_new = W - alpha * riem_grad;

        % Retraction (project back onto manifold)
        W = (1/sqrt(Nt)) * W_new ./ abs(W_new);
    end
end

function I_total = compute_total_interference(H, W, F, U, K, S)
    I_total = 0;
    for q = 1:U
        for k = 1:K
            for u = 1:U
                for i = 1:K
                    if ~(u==q && i==k)
                        S_coef = S(abs(i-k)+1);
                        I_total = I_total + abs(S_coef)^2 * ...
                            abs(squeeze(H(q,:,i)) * W * F(:,u,i))^2;
                    end
                end
            end
        end
    end
end

% Main Hybrid Beamforming Optimization
function [W_opt, F_opt, R_history] = hybrid_beamforming_sumrate(H, Nt, NRF, U, K, sigma2, S, max_outer_iter)
    % Initialize analog beamformer (random phases)
    W = (1/sqrt(Nt)) * exp(1j * 2 * pi * rand(Nt, NRF));

    R_history = zeros(max_outer_iter, 1);

    for outer = 1:max_outer_iter
        % Step 1: Optimize digital beamformer (ZF)
        F = compute_ZF_precoder(H, W, K, U);

        % Step 2: Optimize analog beamformer (Manifold)
        W = optimize_analog_manifold(W, H, F, U, K, sigma2, S, 50, 0.01);

        % Calculate sum rate
        R_history(outer) = compute_sum_rate(H, W, F, U, K, sigma2, S);

        % Check convergence
        if outer > 1 && abs(R_history(outer) - R_history(outer-1)) < 1e-4
            break;
        end
    end

    W_opt = W;
    F_opt = F;
end
```

### 5.2 Problem 2: Power Minimization (CVX)

```matlab
%% Power Minimization with SINR Constraints using CVX
% Requires CVX toolbox: https://cvxr.com/cvx/

function [W_opt, P_min] = power_minimization_cvx(H, gamma_min, Nt, U)
    % H: Channel matrix (U x Nt)
    % gamma_min: Minimum SINR requirement

    cvx_begin sdp quiet
        variable W(Nt, U) complex

        % Objective: Minimize total power
        minimize(sum(sum_square_abs(W)))

        subject to
            for q = 1:U
                % SINR constraint for each user
                signal_power = square_abs(H(q,:) * W(:,q));
                interference_power = 0;
                for u = 1:U
                    if u ~= q
                        interference_power = interference_power + ...
                            square_abs(H(q,:) * W(:,u));
                    end
                end
                signal_power >= gamma_min * (interference_power + sigma2);
            end
    cvx_end

    W_opt = W;
    P_min = sum(sum(abs(W).^2));
end

%% SDR Approach for Power Minimization
function [w_opt, P_min] = power_min_SDR(H, gamma_min, sigma2, Nt, U)
    % Using Semidefinite Relaxation

    cvx_begin sdp quiet
        variable W_sdp(Nt, Nt, U) hermitian semidefinite

        % Objective
        obj = 0;
        for u = 1:U
            obj = obj + trace(W_sdp(:,:,u));
        end
        minimize(obj)

        subject to
            for q = 1:U
                hq = H(q,:).';
                Hq = hq * hq';

                signal = trace(Hq * W_sdp(:,:,q));
                interference = 0;
                for u = 1:U
                    if u ~= q
                        interference = interference + trace(Hq * W_sdp(:,:,u));
                    end
                end

                % SINR constraint
                signal >= gamma_min * (interference + sigma2);
            end
    cvx_end

    % Extract rank-1 solution
    w_opt = zeros(Nt, U);
    for u = 1:U
        [V, D] = eig(W_sdp(:,:,u));
        [~, idx] = max(diag(D));
        w_opt(:,u) = sqrt(D(idx,idx)) * V(:,idx);
    end

    P_min = sum(sum(abs(w_opt).^2));
end
```

### 5.3 Problem 3: Max-Min Fairness

```matlab
%% Max-Min Fairness Beamforming
function [W_opt, R_min] = maxmin_fairness(H, P_total, Nt, U, sigma2)
    % Bisection method combined with feasibility check

    R_low = 0;
    R_high = log2(1 + P_total * max(sum(abs(H).^2, 2)) / sigma2);
    tolerance = 0.01;

    while (R_high - R_low) > tolerance
        R_target = (R_low + R_high) / 2;
        gamma_target = 2^R_target - 1;

        % Check feasibility
        [is_feasible, W_feas] = check_feasibility(H, gamma_target, P_total, Nt, U, sigma2);

        if is_feasible
            R_low = R_target;
            W_opt = W_feas;
        else
            R_high = R_target;
        end
    end

    R_min = R_low;
end

function [is_feasible, W] = check_feasibility(H, gamma_target, P_total, Nt, U, sigma2)
    % Solve power minimization with target SINR
    % If P_min <= P_total, it's feasible

    cvx_begin sdp quiet
        variable W(Nt, U) complex

        minimize(sum(sum_square_abs(W)))

        subject to
            for q = 1:U
                signal_power = square_abs(H(q,:) * W(:,q));
                interference_power = 0;
                for u = 1:U
                    if u ~= q
                        interference_power = interference_power + ...
                            square_abs(H(q,:) * W(:,u));
                    end
                end
                signal_power >= gamma_target * (interference_power + sigma2);
            end
    cvx_end

    P_min = sum(sum(abs(W).^2));
    is_feasible = (cvx_status == "Solved") && (P_min <= P_total);
end

%% Alternative: Epigraph Formulation for Max-Min
function [W_opt, t_opt] = maxmin_epigraph(H, P_total, Nt, U, sigma2)
    cvx_begin
        variable W(Nt, U) complex
        variable t nonnegative

        maximize(t)

        subject to
            % Each user's SINR >= t
            for q = 1:U
                signal_power = square_abs(H(q,:) * W(:,q));
                interference_power = 0;
                for u = 1:U
                    if u ~= q
                        interference_power = interference_power + ...
                            square_abs(H(q,:) * W(:,u));
                    end
                end
                signal_power >= t * (interference_power + sigma2);
            end

            % Total power constraint
            sum(sum_square_abs(W)) <= P_total;
    cvx_end

    W_opt = W;
    t_opt = t;
end
```

### 5.4 Problem 4: Weighted Sum Rate (WMMSE)

```matlab
%% WMMSE Algorithm for Weighted Sum Rate Maximization
function [W_opt, R_weighted] = WMMSE_algorithm(H, alpha, P_total, Nt, U, sigma2, max_iter)
    % alpha: user weights (U x 1)

    % Initialize beamformers (MRT)
    W = zeros(Nt, U);
    for u = 1:U
        W(:,u) = H(u,:)' / norm(H(u,:));
    end
    W = sqrt(P_total/U) * W;

    for iter = 1:max_iter
        % Step 1: Update receive filters (MMSE receivers)
        u_rx = zeros(U, 1);
        for q = 1:U
            signal = H(q,:) * W(:,q);
            interference_plus_noise = sigma2;
            for u = 1:U
                interference_plus_noise = interference_plus_noise + abs(H(q,:) * W(:,u))^2;
            end
            u_rx(q) = conj(signal) / interference_plus_noise;
        end

        % Step 2: Update MSE weights
        e = zeros(U, 1);  % MSE
        w_mse = zeros(U, 1);  % MSE weights
        for q = 1:U
            e(q) = abs(1 - u_rx(q) * H(q,:) * W(:,q))^2;
            for u = 1:U
                if u ~= q
                    e(q) = e(q) + abs(u_rx(q) * H(q,:) * W(:,u))^2;
                end
            end
            e(q) = e(q) + abs(u_rx(q))^2 * sigma2;
            w_mse(q) = 1 / e(q);
        end

        % Step 3: Update transmit beamformers
        % Solve: min sum_q alpha_q * w_q * e_q  s.t. ||W||^2 <= P

        % Lagrangian method (simplified)
        A = zeros(Nt, Nt);
        b = zeros(Nt, U);
        for q = 1:U
            A = A + alpha(q) * w_mse(q) * abs(u_rx(q))^2 * (H(q,:)' * H(q,:));
            b(:,q) = alpha(q) * w_mse(q) * conj(u_rx(q)) * H(q,:)';
        end

        % Water-filling like solution with power constraint
        lambda = 0;  % Lagrange multiplier for power
        lambda_low = 0;
        lambda_high = 100;

        for bisect = 1:20
            lambda = (lambda_low + lambda_high) / 2;
            W_temp = (A + lambda * eye(Nt)) \ b;
            P_used = norm(W_temp, 'fro')^2;

            if P_used > P_total
                lambda_low = lambda;
            else
                lambda_high = lambda;
            end
        end
        W = W_temp;
    end

    W_opt = W;

    % Calculate weighted sum rate
    R_weighted = 0;
    for q = 1:U
        SINR = abs(H(q,:) * W(:,q))^2 / ...
               (sum(abs(H(q,:) * W).^2) - abs(H(q,:) * W(:,q))^2 + sigma2);
        R_weighted = R_weighted + alpha(q) * log2(1 + SINR);
    end
end
```

---

## 6. Python Implementation

### 6.1 Setup and Dependencies

```python
"""
Optimization Problem Formulations for MIMO-OFDM Systems
Dependencies: numpy, scipy, cvxpy, pymanopt
"""

import numpy as np
from scipy.linalg import pinv, eig
import cvxpy as cp

# For manifold optimization
try:
    import pymanopt
    from pymanopt.manifolds import ComplexCircle
    from pymanopt import Problem
    from pymanopt.optimizers import ConjugateGradient
    PYMANOPT_AVAILABLE = True
except ImportError:
    PYMANOPT_AVAILABLE = False
    print("pymanopt not installed. Install with: pip install pymanopt")
```

### 6.2 Channel Generation

```python
def generate_sv_channel(Nt, U, K, Nc=3, Nray=10, d=5, f=0.3e12):
    """
    Generate Saleh-Valenzuela THz channel

    Parameters:
    -----------
    Nt : int - Number of transmit antennas
    U : int - Number of users
    K : int - Number of subcarriers
    Nc : int - Number of clusters
    Nray : int - Number of rays per cluster
    d : float - Distance in meters
    f : float - Frequency in Hz

    Returns:
    --------
    H : ndarray - Channel matrix (U, Nt, K)
    """
    c = 3e8  # Speed of light
    lambda_wave = c / f
    kabs = 0.0033  # Absorption coefficient at 0.3 THz

    # Path loss
    path_loss = (c / (4 * np.pi * f * d))**2 * np.exp(-kabs * d)

    H = np.zeros((U, Nt, K), dtype=complex)

    for q in range(U):
        for k in range(K):
            h_qk = np.zeros(Nt, dtype=complex)
            for c_idx in range(Nc):
                for r in range(Nray):
                    # Random AoD
                    phi = 2 * np.pi * np.random.rand()
                    theta = np.pi * np.random.rand() - np.pi/2

                    # Channel gain
                    alpha = np.sqrt(path_loss) * (np.random.randn() + 1j * np.random.randn()) / np.sqrt(2)

                    # Array response (ULA)
                    n = np.arange(Nt)
                    a = np.exp(1j * np.pi * n * np.sin(phi)) / np.sqrt(Nt)

                    h_qk += np.sqrt(Nt / (Nc * Nray)) * alpha * a

            H[q, :, k] = h_qk

    return H
```

### 6.3 Digital Beamformer (Zero-Forcing)

```python
def compute_zf_precoder(H, W, K, U, NRF):
    """
    Compute Zero-Forcing digital precoder

    Parameters:
    -----------
    H : ndarray - Channel matrix (U, Nt, K)
    W : ndarray - Analog beamformer (Nt, NRF)
    K : int - Number of subcarriers
    U : int - Number of users
    NRF : int - Number of RF chains

    Returns:
    --------
    F : ndarray - Digital precoder (NRF, U, K)
    """
    F = np.zeros((NRF, U, K), dtype=complex)

    for k in range(K):
        Hk = H[:, :, k]  # U x Nt
        He = Hk @ W      # Effective channel: U x NRF

        # ZF precoder: F = H_e^H (H_e H_e^H)^{-1}
        try:
            F_zf = He.conj().T @ pinv(He @ He.conj().T)
        except:
            F_zf = pinv(He)

        # Normalize per user
        for u in range(U):
            norm_factor = np.linalg.norm(W @ F_zf[:, u])
            if norm_factor > 1e-10:
                F[:, u, k] = F_zf[:, u] / norm_factor
            else:
                F[:, u, k] = F_zf[:, u]

    return F
```

### 6.4 Problem 1: Sum Rate Maximization

```python
def compute_sinr(H, W, F, q, k, U, K, sigma2, S):
    """Compute SINR for user q at subcarrier k with ICI"""
    Hqk = H[q, :, k]

    # Desired signal
    signal = np.abs(S[0] * Hqk @ W @ F[:, q, k])**2

    # Interference (inter-user + ICI)
    interference = 0
    for u in range(U):
        for i in range(K):
            if not (u == q and i == k):
                s_idx = min(abs(i - k), len(S) - 1)
                s_coef = S[s_idx]
                interference += np.abs(s_coef)**2 * np.abs(H[q, :, i] @ W @ F[:, u, i])**2

    return signal / (interference + sigma2)


def compute_sum_rate(H, W, F, U, K, sigma2, S):
    """Compute total sum rate"""
    R_sum = 0
    for q in range(U):
        for k in range(K):
            sinr = compute_sinr(H, W, F, q, k, U, K, sigma2, S)
            R_sum += np.log2(1 + sinr)
    return R_sum


def compute_total_interference(H, W, F, U, K, S):
    """Compute total interference (for minimization)"""
    I_total = 0
    for q in range(U):
        for k in range(K):
            for u in range(U):
                for i in range(K):
                    if not (u == q and i == k):
                        s_idx = min(abs(i - k), len(S) - 1)
                        s_coef = S[s_idx]
                        I_total += np.abs(s_coef)**2 * np.abs(H[q, :, i] @ W @ F[:, u, i])**2
    return I_total


class ManifoldOptimizer:
    """Riemannian Manifold Optimization for Analog Beamformer"""

    def __init__(self, Nt, NRF):
        self.Nt = Nt
        self.NRF = NRF

    def retraction(self, W):
        """Project onto constant modulus manifold"""
        return (1 / np.sqrt(self.Nt)) * W / np.abs(W)

    def riemannian_gradient(self, W, euclidean_grad):
        """Project Euclidean gradient onto tangent space"""
        return euclidean_grad - np.real(euclidean_grad * np.conj(W)) * W

    def optimize(self, W_init, H, F, U, K, S, alpha=0.01, max_iter=100, tol=1e-6):
        """Conjugate Gradient on Manifold"""
        W = W_init.copy()

        for iteration in range(max_iter):
            # Compute numerical gradient
            grad = np.zeros_like(W, dtype=complex)
            eps = 1e-7

            for m in range(self.Nt):
                for n in range(self.NRF):
                    W_plus = W.copy()
                    W_plus[m, n] += eps
                    W_minus = W.copy()
                    W_minus[m, n] -= eps

                    I_plus = compute_total_interference(H, W_plus, F, U, K, S)
                    I_minus = compute_total_interference(H, W_minus, F, U, K, S)
                    grad[m, n] = (I_plus - I_minus) / (2 * eps)

            # Riemannian gradient
            riem_grad = self.riemannian_gradient(W, grad)

            # Update with step size
            W_new = W - alpha * riem_grad

            # Retraction
            W = self.retraction(W_new)

            # Check convergence
            if np.linalg.norm(riem_grad) < tol:
                break

        return W


def hybrid_beamforming_sumrate(H, Nt, NRF, U, K, sigma2, S, max_outer_iter=20):
    """
    Main function for sum rate maximization with hybrid beamforming

    Returns:
    --------
    W_opt : Analog beamformer
    F_opt : Digital beamformer
    R_history : Convergence history
    """
    # Initialize analog beamformer (random phases)
    W = (1 / np.sqrt(Nt)) * np.exp(1j * 2 * np.pi * np.random.rand(Nt, NRF))

    R_history = []
    optimizer = ManifoldOptimizer(Nt, NRF)

    for outer in range(max_outer_iter):
        # Step 1: Optimize digital beamformer (ZF)
        F = compute_zf_precoder(H, W, K, U, NRF)

        # Step 2: Optimize analog beamformer (Manifold)
        W = optimizer.optimize(W, H, F, U, K, S)

        # Calculate sum rate
        R = compute_sum_rate(H, W, F, U, K, sigma2, S)
        R_history.append(R)

        print(f"Iteration {outer + 1}: Sum Rate = {R:.4f} bps/Hz")

        # Check convergence
        if outer > 0 and abs(R_history[-1] - R_history[-2]) < 1e-4:
            break

    return W, F, R_history
```

### 6.5 Problem 2: Power Minimization (CVXPY)

```python
def power_minimization_cvxpy(H, gamma_min, sigma2, Nt, U):
    """
    Power minimization with SINR constraints using CVXPY

    Parameters:
    -----------
    H : ndarray - Channel matrix (U, Nt) for single subcarrier
    gamma_min : float - Minimum SINR requirement
    sigma2 : float - Noise variance

    Returns:
    --------
    W_opt : Optimal beamformer
    P_min : Minimum power
    """
    # Decision variable
    W = cp.Variable((Nt, U), complex=True)

    # Objective: minimize total power
    objective = cp.Minimize(cp.sum(cp.sum_squares(cp.abs(W))))

    # Constraints
    constraints = []
    for q in range(U):
        hq = H[q, :]

        # Signal power
        signal = cp.abs(hq @ W[:, q])**2

        # Interference power
        interference = 0
        for u in range(U):
            if u != q:
                interference = interference + cp.abs(hq @ W[:, u])**2

        # SINR constraint: signal >= gamma_min * (interference + noise)
        constraints.append(signal >= gamma_min * (interference + sigma2))

    # Solve
    prob = cp.Problem(objective, constraints)
    prob.solve(solver=cp.SCS)

    if prob.status == cp.OPTIMAL:
        W_opt = W.value
        P_min = np.sum(np.abs(W_opt)**2)
        return W_opt, P_min
    else:
        print(f"Problem status: {prob.status}")
        return None, None


def power_minimization_sdr(H, gamma_min, sigma2, Nt, U):
    """
    Power minimization using Semidefinite Relaxation (SDR)
    """
    # Decision variables: W_u = w_u @ w_u^H (positive semidefinite)
    W_sdp = [cp.Variable((Nt, Nt), hermitian=True) for _ in range(U)]

    # Objective: minimize sum of traces
    objective = cp.Minimize(sum(cp.trace(W_sdp[u]) for u in range(U)))

    # Constraints
    constraints = []

    # PSD constraints
    for u in range(U):
        constraints.append(W_sdp[u] >> 0)

    # SINR constraints
    for q in range(U):
        hq = H[q, :].reshape(-1, 1)
        Hq = hq @ hq.conj().T

        signal = cp.trace(Hq @ W_sdp[q])
        interference = sum(cp.trace(Hq @ W_sdp[u]) for u in range(U) if u != q)

        constraints.append(signal >= gamma_min * (interference + sigma2))

    # Solve
    prob = cp.Problem(objective, constraints)
    prob.solve(solver=cp.SCS)

    if prob.status == cp.OPTIMAL:
        # Extract rank-1 solutions via eigendecomposition
        w_opt = np.zeros((Nt, U), dtype=complex)
        for u in range(U):
            W_val = W_sdp[u].value
            eigenvalues, eigenvectors = eig(W_val)
            idx = np.argmax(np.real(eigenvalues))
            w_opt[:, u] = np.sqrt(np.real(eigenvalues[idx])) * eigenvectors[:, idx]

        P_min = np.sum(np.abs(w_opt)**2)
        return w_opt, P_min
    else:
        return None, None
```

### 6.6 Problem 3: Max-Min Fairness

```python
def maxmin_fairness_bisection(H, P_total, sigma2, Nt, U):
    """
    Max-min fairness using bisection method

    Returns:
    --------
    W_opt : Optimal beamformer
    R_min : Maximum minimum rate
    """
    # Bounds for bisection
    R_low = 0
    channel_norms = np.sum(np.abs(H)**2, axis=1)
    R_high = np.log2(1 + P_total * np.max(channel_norms) / sigma2)

    tolerance = 0.01
    W_opt = None

    while (R_high - R_low) > tolerance:
        R_target = (R_low + R_high) / 2
        gamma_target = 2**R_target - 1

        # Check feasibility
        is_feasible, W_feas = check_feasibility_cvxpy(H, gamma_target, P_total, sigma2, Nt, U)

        if is_feasible:
            R_low = R_target
            W_opt = W_feas
        else:
            R_high = R_target

    return W_opt, R_low


def check_feasibility_cvxpy(H, gamma_target, P_total, sigma2, Nt, U):
    """Check if target SINR is achievable within power budget"""
    W = cp.Variable((Nt, U), complex=True)

    objective = cp.Minimize(cp.sum(cp.sum_squares(cp.abs(W))))

    constraints = []
    for q in range(U):
        hq = H[q, :]
        signal = cp.abs(hq @ W[:, q])**2
        interference = sum(cp.abs(hq @ W[:, u])**2 for u in range(U) if u != q)
        constraints.append(signal >= gamma_target * (interference + sigma2))

    prob = cp.Problem(objective, constraints)

    try:
        prob.solve(solver=cp.SCS, verbose=False)
        if prob.status == cp.OPTIMAL:
            P_min = prob.value
            return P_min <= P_total, W.value
    except:
        pass

    return False, None


def maxmin_fairness_epigraph(H, P_total, sigma2, Nt, U):
    """
    Max-min fairness using epigraph formulation

    maximize t
    s.t. SINR_q >= t, for all q
         total_power <= P_total
    """
    W = cp.Variable((Nt, U), complex=True)
    t = cp.Variable(nonneg=True)

    objective = cp.Maximize(t)

    constraints = []

    # SINR constraints
    for q in range(U):
        hq = H[q, :]
        signal = cp.abs(hq @ W[:, q])**2
        interference = sum(cp.abs(hq @ W[:, u])**2 for u in range(U) if u != q)
        constraints.append(signal >= t * (interference + sigma2))

    # Power constraint
    constraints.append(cp.sum(cp.sum_squares(cp.abs(W))) <= P_total)

    prob = cp.Problem(objective, constraints)
    prob.solve(solver=cp.SCS)

    if prob.status == cp.OPTIMAL:
        return W.value, t.value, np.log2(1 + t.value)
    else:
        return None, None, None
```

### 6.7 Problem 4: Weighted Sum Rate (WMMSE)

```python
def wmmse_algorithm(H, alpha, P_total, sigma2, Nt, U, max_iter=100):
    """
    WMMSE Algorithm for Weighted Sum Rate Maximization

    Parameters:
    -----------
    H : ndarray - Channel matrix (U, Nt)
    alpha : ndarray - User weights (U,)
    P_total : float - Total power budget
    sigma2 : float - Noise variance

    Returns:
    --------
    W_opt : Optimal beamformer
    R_weighted : Weighted sum rate
    """
    # Initialize with MRT (Maximum Ratio Transmission)
    W = np.zeros((Nt, U), dtype=complex)
    for u in range(U):
        W[:, u] = H[u, :].conj() / np.linalg.norm(H[u, :])
    W = np.sqrt(P_total / U) * W

    R_history = []

    for iteration in range(max_iter):
        # Step 1: Update receive filters (MMSE)
        u_rx = np.zeros(U, dtype=complex)
        for q in range(U):
            signal = H[q, :] @ W[:, q]
            total_power = sigma2 + np.sum(np.abs(H[q, :] @ W)**2)
            u_rx[q] = np.conj(signal) / total_power

        # Step 2: Compute MSE and update weights
        e = np.zeros(U)
        w_mse = np.zeros(U)
        for q in range(U):
            # MSE = E[|s - u*y|^2]
            e[q] = np.abs(1 - u_rx[q] * H[q, :] @ W[:, q])**2
            for u in range(U):
                if u != q:
                    e[q] += np.abs(u_rx[q] * H[q, :] @ W[:, u])**2
            e[q] += np.abs(u_rx[q])**2 * sigma2

            # Optimal weight
            w_mse[q] = 1 / e[q]

        # Step 3: Update transmit beamformers
        # Solve: min sum_q alpha_q * w_q * e_q s.t. ||W||^2 <= P_total

        # Build matrices for the update
        A = np.zeros((Nt, Nt), dtype=complex)
        B = np.zeros((Nt, U), dtype=complex)

        for q in range(U):
            hq = H[q, :].reshape(-1, 1)
            A += alpha[q] * w_mse[q] * np.abs(u_rx[q])**2 * (hq @ hq.conj().T)
            B[:, q] = alpha[q] * w_mse[q] * np.conj(u_rx[q]) * H[q, :].conj()

        # Water-filling like solution with power constraint (bisection on lambda)
        lambda_low, lambda_high = 0, 100

        for _ in range(30):
            lam = (lambda_low + lambda_high) / 2
            try:
                W_temp = np.linalg.solve(A + lam * np.eye(Nt), B)
            except:
                W_temp = pinv(A + lam * np.eye(Nt)) @ B

            P_used = np.linalg.norm(W_temp, 'fro')**2

            if P_used > P_total:
                lambda_low = lam
            else:
                lambda_high = lam

        W = W_temp

        # Compute weighted sum rate
        R = 0
        for q in range(U):
            signal = np.abs(H[q, :] @ W[:, q])**2
            interference = np.sum(np.abs(H[q, :] @ W)**2) - signal
            sinr = signal / (interference + sigma2)
            R += alpha[q] * np.log2(1 + sinr)

        R_history.append(R)

        # Check convergence
        if iteration > 0 and abs(R_history[-1] - R_history[-2]) < 1e-6:
            break

    return W, R_history[-1], R_history
```

### 6.8 Complete Example Usage

```python
def main():
    """Complete example demonstrating all four optimization problems"""

    # System parameters
    Nt = 64          # Transmit antennas
    NRF = 9          # RF chains
    U = 4            # Number of users (reduced for demo)
    K = 16           # Subcarriers (reduced for demo)
    P_total = 1.0    # Total power
    sigma2 = 0.01    # Noise variance

    # ICI coefficients (simplified)
    S = np.array([1.0, 0.1, 0.05, 0.02])

    print("=" * 60)
    print("Optimization Problem Formulations Demo")
    print("=" * 60)

    # Generate channel
    print("\nGenerating channel...")
    H_full = generate_sv_channel(Nt, U, K)
    H_single = H_full[:, :, 0]  # Single subcarrier for some problems

    # ========================================
    # Problem 1: Sum Rate Maximization
    # ========================================
    print("\n" + "-" * 40)
    print("PROBLEM 1: Sum Rate Maximization")
    print("-" * 40)

    W_sr, F_sr, R_hist = hybrid_beamforming_sumrate(
        H_full, Nt, NRF, U, K, sigma2, S, max_outer_iter=10
    )
    print(f"Final Sum Rate: {R_hist[-1]:.4f} bps/Hz")

    # ========================================
    # Problem 2: Power Minimization
    # ========================================
    print("\n" + "-" * 40)
    print("PROBLEM 2: Power Minimization")
    print("-" * 40)

    gamma_min = 5.0  # Minimum SINR = 5 (linear)
    W_pm, P_min = power_minimization_cvxpy(H_single, gamma_min, sigma2, Nt, U)
    if W_pm is not None:
        print(f"Minimum Power: {P_min:.4f}")
        print(f"Target SINR: {gamma_min} ({10*np.log10(gamma_min):.2f} dB)")
    else:
        print("Problem infeasible with given SINR constraint")

    # ========================================
    # Problem 3: Max-Min Fairness
    # ========================================
    print("\n" + "-" * 40)
    print("PROBLEM 3: Max-Min Fairness")
    print("-" * 40)

    W_mm, R_min = maxmin_fairness_bisection(H_single, P_total, sigma2, Nt, U)
    if W_mm is not None:
        print(f"Maximum Minimum Rate: {R_min:.4f} bps/Hz")

        # Verify all users achieve similar rates
        rates = []
        for q in range(U):
            signal = np.abs(H_single[q, :] @ W_mm[:, q])**2
            interference = np.sum(np.abs(H_single[q, :] @ W_mm)**2) - signal
            sinr = signal / (interference + sigma2)
            rates.append(np.log2(1 + sinr))
        print(f"User rates: {[f'{r:.3f}' for r in rates]}")

    # ========================================
    # Problem 4: Weighted Sum Rate (WMMSE)
    # ========================================
    print("\n" + "-" * 40)
    print("PROBLEM 4: Weighted Sum Rate (WMMSE)")
    print("-" * 40)

    # Equal weights
    alpha = np.ones(U) / U
    W_wsr, R_weighted, R_hist_wmmse = wmmse_algorithm(
        H_single, alpha, P_total, sigma2, Nt, U, max_iter=50
    )
    print(f"Weighted Sum Rate: {R_weighted:.4f} bps/Hz")

    # Priority weights (user 0 has higher priority)
    alpha_priority = np.array([0.4, 0.3, 0.2, 0.1])
    W_wsr_p, R_weighted_p, _ = wmmse_algorithm(
        H_single, alpha_priority, P_total, sigma2, Nt, U, max_iter=50
    )
    print(f"Weighted Sum Rate (prioritized): {R_weighted_p:.4f} bps/Hz")

    print("\n" + "=" * 60)
    print("Demo Complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
```

---

## 7. Key References

### Foundational Papers

1. **Hybrid Beamforming**: Yuan et al., "Hybrid Beamforming for MIMO-OFDM Terahertz Wireless Systems," IEEE Trans. Commun., 2020.

2. **Manifold Optimization**: Yu et al., "Alternating Minimization Algorithms for Hybrid Precoding in Millimeter Wave MIMO Systems," IEEE JSTSP, 2016.

3. **WMMSE Algorithm**: Shi et al., "An Iteratively Weighted MMSE Approach to Distributed Sum-Utility Maximization," IEEE Trans. Signal Process., 2011.

4. **Optimal Beamforming**: Björnson et al., "Optimal Multiuser Transmit Beamforming: A Difficult Problem with a Simple Solution Structure," IEEE SPM, 2014.

5. **SDR Techniques**: Luo et al., "Semidefinite Relaxation of Quadratic Optimization Problems," IEEE SPM, 2010.

### Implementation Resources

- **CVX (MATLAB)**: https://cvxr.com/cvx/
- **CVXPY (Python)**: https://www.cvxpy.org/
- **Pymanopt**: https://pymanopt.org/
- **Geoopt (PyTorch)**: https://github.com/geoopt/geoopt
- **Emil Björnson's Code**: https://github.com/emilbjornson/optimal-beamforming

---

## Summary Table: Problem Formulations

| Problem | Objective | Key Constraints | Solution Method |
|---------|-----------|-----------------|-----------------|
| **Sum Rate Max** | max Σ log(1+SINR) | Power, Constant Modulus | Alternating Opt + Manifold |
| **Power Min** | min Σ\|\|w\|\|² | SINR ≥ γ_min | SDR, CVX |
| **Max-Min Fair** | max min_q R_q | Power | Bisection + Feasibility |
| **Weighted Sum** | max Σ α_q R_q | Power, Rate Guarantees | WMMSE |

---

*Document created for optimization problem formulation study in wireless communications.*
*Last updated: January 2026*
