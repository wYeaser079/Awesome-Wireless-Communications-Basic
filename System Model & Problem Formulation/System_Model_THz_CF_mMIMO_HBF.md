# System Model: THz Cell-Free Massive MIMO with Hybrid Beamforming

## II. SYSTEM MODEL AND PROBLEM FORMULATION

### A. System Architecture

We consider a downlink cell-free massive MIMO (CF-mMIMO) system operating in the
terahertz (THz) band, as illustrated in Fig. 1. The network consists of L access
points (APs) that coherently serve K single-antenna user equipments (UEs) over
the same time-frequency resources, where L >> K. Each AP is equipped with M
antennas and N_RF radio frequency (RF) chains, where K <= N_RF << M.
All APs are connected to a central processing unit (CPU) via high-capacity
fronthaul links. The CPU is responsible for channel state information (CSI)
collection, beamforming design, and data distribution. The system operates in
time-division duplex (TDD) mode, exploiting channel reciprocity for downlink
CSI acquisition from uplink pilot training.

To support wideband transmission in the THz band, orthogonal frequency-division
multiplexing (OFDM) is employed with N subcarriers over a total bandwidth B.
The center frequency is f_c, the subcarrier spacing is Delta_f = B/N, and
the frequency of the n-th subcarrier is:

    f_n = f_c + (n - (N-1)/2) * Delta_f,    n = 0, 1, ..., N-1         (1)

Each AP employs a hybrid beamforming (HBF) architecture comprising a
frequency-flat analog precoder and frequency-selective digital precoders,
which significantly reduces hardware cost and power consumption compared
to fully digital architectures at THz frequencies.


### B. THz Channel Model

The THz propagation environment differs fundamentally from conventional
microwave and mmWave bands due to severe free-space spreading loss,
frequency-dependent molecular absorption, and highly sparse multipath
structure. We adopt the extended Saleh-Valenzuela channel model adapted
for THz frequencies.

#### B.1 THz Path Gain

The path gain from AP l to UE k at frequency f_n incorporates both
free-space spreading loss and molecular absorption loss:

    G_{l,k}(f_n) = G_spread(f_n, d_{l,k}) * G_abs(f_n, d_{l,k})       (2)

where d_{l,k} is the distance between AP l and UE k.

The free-space spreading loss (inverse Friis formula) is:

    G_spread(f_n, d_{l,k}) = ( c / (4 * pi * f_n * d_{l,k}) )^2        (3)

where c = 3 x 10^8 m/s is the speed of light.

The molecular absorption loss follows the Beer-Lambert law:

    G_abs(f_n, d_{l,k}) = exp( -kappa_abs(f_n) * d_{l,k} )             (4)

where kappa_abs(f_n) [m^{-1}] is the frequency-dependent molecular
absorption coefficient, which can be computed from the HITRAN spectroscopic
database for given atmospheric conditions (temperature, pressure, humidity).

Therefore, the total THz path gain is:

    G_{l,k}(f_n) = ( c / (4*pi*f_n*d_{l,k}) )^2 * exp(-kappa_abs(f_n)*d_{l,k})   (5)

Remark 1: Unlike lower frequency bands where the path gain is approximately
constant across the bandwidth, the THz path gain G_{l,k}(f_n) varies across
subcarriers due to: (i) the f_n^{-2} dependence in the spreading loss, and
(ii) the frequency-dependent molecular absorption kappa_abs(f_n). This
frequency selectivity is a distinguishing feature of THz wideband systems.

#### B.2 THz Multipath Channel with Beam Squint

The THz channel is characterized by a dominant line-of-sight (LoS) path and
a small number of non-LoS (NLoS) scattered paths, modeled using Rician fading.
The channel vector from AP l to UE k on subcarrier n is:

    h_{l,k}[n] = sqrt(G_{l,k}(f_n)) * ( sqrt(kappa_R / (kappa_R + 1)) * h_{l,k}^{LoS}[n]
                                        + sqrt(1 / (kappa_R + 1)) * h_{l,k}^{NLoS}[n] )     (6)

where kappa_R is the Rician K-factor (linear scale), which is typically high
at THz frequencies (kappa_R = 10-20 dB) due to the sparse scattering environment.

The LoS component is deterministic and given by:

    h_{l,k}^{LoS}[n] = a_M(phi_{l,k}^{LoS}, f_n) * exp(-j*2*pi*f_n*tau_{l,k}^{LoS})     (7)

where phi_{l,k}^{LoS} is the angle of departure (AoD) of the LoS path from
AP l to UE k, and tau_{l,k}^{LoS} = d_{l,k}/c is the LoS propagation delay.

The NLoS component captures P scattered multipath clusters:

    h_{l,k}^{NLoS}[n] = (1/sqrt(P)) * sum_{p=1}^{P} alpha_{l,k,p} * a_M(phi_{l,k,p}, f_n)
                          * exp(-j*2*pi*f_n*tau_{l,k,p})                                     (8)

where alpha_{l,k,p} ~ CN(0, 1) is the complex gain of the p-th NLoS path,
phi_{l,k,p} is the AoD of the p-th path, and tau_{l,k,p} is the
corresponding propagation delay.

#### B.3 Frequency-Dependent Array Response (Beam Squint)

For a uniform linear array (ULA) with M elements and half-wavelength spacing
at the center frequency (d_s = c/(2*f_c) = lambda_c/2), the array response
vector at frequency f_n is:

    a_M(phi, f_n) = (1/sqrt(M)) * [1, e^{j*pi*(f_n/f_c)*sin(phi)}, ...,
                     e^{j*pi*(M-1)*(f_n/f_c)*sin(phi)}]^T                                   (9)

The critical feature is the frequency ratio xi_n = f_n / f_c, which causes
the effective spatial frequency to scale across subcarriers. At the center
frequency (f_n = f_c), this reduces to the standard array response:

    a_M(phi, f_c) = (1/sqrt(M)) * [1, e^{j*pi*sin(phi)}, ...,
                     e^{j*pi*(M-1)*sin(phi)}]^T                                             (10)

For edge subcarriers where f_n != f_c, the beam direction deviates from the
intended direction, a phenomenon known as beam squint. The fractional
bandwidth parameter beta = B/f_c quantifies this effect. For THz systems
(e.g., f_c = 300 GHz, B = 30 GHz), beta = 0.1, meaning the beam direction
can shift by up to 10% of the spatial frequency — a severe degradation
that fundamentally limits the performance of frequency-flat analog precoders.

Remark 2: The beam squint effect is naturally captured by the frequency-dependent
array response in (9). Since the analog precoder F_{RF,l} is designed for the
center frequency f_c but the actual channel response varies with f_n, the
beamforming gain degrades on edge subcarriers. The digital precoder F_{BB,l}[n]
can partially compensate for this mismatch, but its ability is limited by the
number of RF chains N_RF.


### C. Hybrid Beamforming Architecture

Each AP l employs a hybrid beamforming architecture consisting of an analog
precoder F_{RF,l} and per-subcarrier digital precoders F_{BB,l}[n].

#### C.1 Analog Precoder

The analog precoder F_{RF,l} in C^{M x N_RF} is implemented using a network
of phase shifters and is shared across all subcarriers (frequency-flat). Each
element satisfies the constant-modulus constraint:

    |[F_{RF,l}]_{m,r}| = 1/sqrt(M),   for all m = 1,...,M,  r = 1,...,N_RF,  l = 1,...,L    (11)

The factor 1/sqrt(M) ensures proper power normalization.

#### C.2 Digital Precoder

The digital precoder on subcarrier n for AP l is:

    F_{BB,l}[n] = [f_{BB,l,1}[n], f_{BB,l,2}[n], ..., f_{BB,l,K}[n]]  in  C^{N_RF x K}    (12)

where f_{BB,l,k}[n] in C^{N_RF x 1} is the digital beamforming vector from
AP l to UE k on subcarrier n.

#### C.3 Transmitted Signal

The transmitted signal from AP l on subcarrier n is:

    x_l[n] = F_{RF,l} * F_{BB,l}[n] * s[n]
           = F_{RF,l} * sum_{k=1}^{K} f_{BB,l,k}[n] * s_k[n]                               (13)

where s[n] = [s_1[n], ..., s_K[n]]^T is the vector of transmitted data symbols
with E[s[n]*s^H[n]] = I_K.

#### C.4 Per-AP Power Constraint

The average transmit power of AP l across all subcarriers must satisfy:

    (1/N) * sum_{n=0}^{N-1} ||F_{RF,l} * F_{BB,l}[n]||_F^2  <=  P_l,    for all l = 1,...,L   (14)

where P_l is the maximum transmit power budget of AP l.

Note: Since F_{RF,l}^H * F_{RF,l} != I_{N_RF} in general (the columns of the
analog precoder are not orthogonal), the power constraint couples the analog
and digital precoders, making the optimization more challenging.


### D. Phase Noise Model and Inter-Carrier Interference

At THz frequencies, oscillator imperfections cause significant phase noise that
degrades OFDM performance. In cell-free systems, each AP has an independent
local oscillator, making phase noise particularly challenging as it varies
across APs.

#### D.1 Wiener Process Phase Noise Model

The phase noise at AP l is modeled as a discrete-time Wiener process over the
N time-domain samples within one OFDM symbol:

    theta_l[t] = theta_l[t-1] + Delta_theta_l[t],    t = 0, 1, ..., N-1                     (15)

where Delta_theta_l[t] ~ N(0, sigma_{theta,l}^2) are independent and identically
distributed (i.i.d.) Gaussian increments with variance:

    sigma_{theta,l}^2 = 4 * pi * beta_{3dB,l} / (N * Delta_f)                               (16)

Here, beta_{3dB,l} [Hz] is the 3-dB linewidth of the local oscillator at AP l.
Similarly, UE k has phase noise psi_k[t] with variance sigma_{psi,k}^2 =
4*pi*beta_{3dB,k}^{UE} / (N*Delta_f).

For simplicity, we define the composite phase noise on the link from AP l to
UE k as:

    Phi_{l,k}[t] = theta_l[t] - psi_k[t]                                                    (17)

with composite innovation variance sigma_{Phi}^2 = sigma_{theta,l}^2 + sigma_{psi,k}^2.

#### D.2 DFT Decomposition: CPE and ICI Coefficients

The phase noise in the time domain creates multiplicative distortion. After
DFT processing at the receiver, this distortion decomposes into a Common
Phase Error (CPE) and Inter-Carrier Interference (ICI).

Define the DFT coefficients of the phase noise process on link (l, k):

    J_{l,k}[p] = (1/N) * sum_{t=0}^{N-1} exp(j*Phi_{l,k}[t]) * exp(-j*2*pi*p*t/N)          (18)

for p = 0, 1, ..., N-1.

- J_{l,k}[0] represents the Common Phase Error (CPE): a complex scalar
  rotation that is approximately constant across subcarriers.

- J_{l,k}[p] for p != 0 represents the ICI coefficients: they cause power
  leakage from subcarrier (n-p) to subcarrier n.

#### D.3 Statistical Properties

For the Wiener phase noise model, the following statistical properties hold:

    E[exp(j*(Phi_{l,k}[t1] - Phi_{l,k}[t2]))] = exp(-sigma_{Phi}^2 * |t1-t2| / 2)          (19)

The power of the CPE coefficient:

    E[|J_{l,k}[0]|^2] = (1/N^2) * sum_{t1=0}^{N-1} sum_{t2=0}^{N-1}
                          exp(-sigma_{Phi}^2*|t1-t2|/2) = 1 - epsilon_{l,k}                  (20)

where epsilon_{l,k} is the fraction of signal power leaked to ICI:

    epsilon_{l,k} = 1 - E[|J_{l,k}[0]|^2]                                                   (21)

For small phase noise variance (sigma_{Phi}^2*N << 1), the first-order
approximation gives:

    epsilon_{l,k} approx (pi * (beta_{3dB,l} + beta_{3dB,k}^{UE})) / Delta_f                (22)

The total ICI power equals: sum_{p != 0} E[|J_{l,k}[p]|^2] = epsilon_{l,k}.

Remark 3: In cell-free systems, the ICI coefficient epsilon_{l,k} is
INDEPENDENT across different APs l, since each AP has its own oscillator.
This is fundamentally different from co-located MIMO where all antennas
share a common oscillator. Consequently, the ICI from different APs adds
incoherently, which makes phase noise compensation more challenging but
also limits the worst-case ICI accumulation.

#### D.4 Assumption on Phase Noise Homogeneity

For tractability, we assume all APs and UEs employ oscillators with identical
3-dB linewidths: beta_{3dB,l} = beta_{3dB}^{AP} for all l, and
beta_{3dB,k}^{UE} = beta_{3dB}^{UE} for all k. The composite ICI parameter
then simplifies to:

    epsilon = pi * (beta_{3dB}^{AP} + beta_{3dB}^{UE}) / Delta_f                            (23)


### E. Downlink Received Signal Model

Combining the channel model, hybrid beamforming, and phase noise, the received
signal at UE k on subcarrier n is obtained by summing the contributions from
all L APs across all N subcarriers (due to ICI):

    y_k[n] = sum_{l=1}^{L} sum_{n'=0}^{N-1} J_{l,k}[<n-n'>_N] * h_{l,k}^H[n']
             * F_{RF,l} * F_{BB,l}[n'] * s[n']  +  n_k[n]                                   (24)

where <.>_N denotes the modulo-N operation and n_k[n] ~ CN(0, sigma^2) is
additive white Gaussian noise (AWGN).

Separating the same-subcarrier term (n' = n) from cross-subcarrier terms
(n' != n), and further separating the desired user from interfering users,
we obtain:

    y_k[n] = J_{l,k}^{sum}[0,n] * s_k[n]      ... (DS: Desired Signal with CPE)
           + IUI_k[n]                            ... (IUI: Inter-User Interference)
           + ICI_k[n]                            ... (ICI: Inter-Carrier Interference)
           + n_k[n]                              ... (Noise)                                  (25)

where each term is defined as follows.

#### E.1 Desired Signal (DS)

    DS_k[n] = sum_{l=1}^{L} J_{l,k}[0] * h_{l,k}^H[n] * F_{RF,l} * f_{BB,l,k}[n] * s_k[n]
                                                                                              (26)

This is the intended signal for UE k on subcarrier n, scaled by the
AP-dependent CPE J_{l,k}[0]. The CPE varies across APs, which impairs
the coherent combining gain.

#### E.2 Inter-User Interference (IUI)

    IUI_k[n] = sum_{l=1}^{L} J_{l,k}[0] * h_{l,k}^H[n]
               * F_{RF,l} * sum_{j=1, j!=k}^{K} f_{BB,l,j}[n] * s_j[n]                     (27)

This is the interference from the K-1 other users on the same subcarrier n,
arising from imperfect spatial separation by the beamformer.

#### E.3 Inter-Carrier Interference (ICI)

    ICI_k[n] = sum_{l=1}^{L} sum_{n'=0, n'!=n}^{N-1} J_{l,k}[<n-n'>_N]
               * h_{l,k}^H[n'] * F_{RF,l} * sum_{j=1}^{K} f_{BB,l,j}[n'] * s_j[n']         (28)

This is the interference leaking from all other subcarriers n' != n to the
current subcarrier n, caused by the phase noise of the THz oscillators. The
ICI couples signals across ALL subcarriers and ALL users.


### F. SINR and Achievable Rate Derivation

#### F.1 Effective Channel Definition

Define the effective baseband channel from AP l to UE k on subcarrier n
after analog beamforming:

    q_{l,k}[n] = h_{l,k}^H[n] * F_{RF,l}    in  C^{1 x N_RF}                               (29)

This effective channel captures both the physical THz channel (including
beam squint through the frequency-dependent array response) and the
frequency-flat analog precoder.

#### F.2 SINR Expression

We assume that the effective channel including CPE can be estimated through
downlink pilot training, and that the ICI is treated as unstructured noise.

The SINR of UE k on subcarrier n is:

                     |sum_{l=1}^{L} q_{l,k}[n] * f_{BB,l,k}[n]|^2
    SINR_k[n] = --------------------------------------------------------                     (30)
                 I_{IUI,k}[n]  +  I_{ICI,k}[n]  +  sigma^2

where the interference terms are:

  Inter-User Interference power:

    I_{IUI,k}[n] = sum_{j=1, j!=k}^{K} |sum_{l=1}^{L} q_{l,k}[n] * f_{BB,l,j}[n]|^2       (31)

  Inter-Carrier Interference power:

    I_{ICI,k}[n] = epsilon * sum_{l=1}^{L} sum_{j=1}^{K} |q_{l,k}[n] * f_{BB,l,j}[n]|^2   (32)

                 = epsilon * sum_{l=1}^{L} ||q_{l,k}[n] * F_{BB,l}[n]||^2                   (32')


Derivation of the ICI power (32):
-----------------------------------
The exact ICI power at UE k on subcarrier n is:
  E[|ICI_k[n]|^2] = sum_l sum_{n'!=n} E[|J_{l,k}[n-n']|^2]
                     * sum_j |q_{l,k}[n']*f_{BB,l,j}[n']|^2

Using the uniform ICI approximation (valid when the channel and beamformers
are relatively flat across subcarriers) and the identity
sum_{p!=0} E[|J_{l,k}[p]|^2] = epsilon, the ICI power simplifies to (32).
-----------------------------------

Remark 4: The SINR in (30) reveals the fundamental tension in THz CF-mMIMO
beamforming design:

  (i) IUI (31) is a COHERENT interference: it depends on the cross-AP
      beamforming coordination (sum over l inside the absolute value).
      It can be suppressed by designing the digital precoders to create
      spatial nulls toward unintended users.

  (ii) ICI (32) is an INCOHERENT interference: it depends on the per-AP
       signal power (sum over l outside the absolute value). Reducing
       ICI requires controlling the total signal power from each AP,
       which may conflict with maximizing the desired signal.

  This IUI-ICI tradeoff is a unique feature of THz CF-mMIMO systems and
  motivates the joint beamforming design proposed in this work.

#### F.3 Compact Matrix Form of SINR

For compact notation, define the stacked effective channel and digital
precoder across all APs:

    q_k[n] = [q_{1,k}[n], q_{2,k}[n], ..., q_{L,k}[n]]  in  C^{1 x L*N_RF}                (33)

    f_k[n] = [f_{BB,1,k}^T[n], f_{BB,2,k}^T[n], ..., f_{BB,L,k}^T[n]]^T  in  C^{L*N_RF x 1}
                                                                                              (34)

Then the SINR can be written as:

                     |q_k[n] * f_k[n]|^2
    SINR_k[n] = ------------------------------------                                         (35)
                 sum_{j!=k} |q_k[n]*f_j[n]|^2
                 + epsilon * sum_{l=1}^{L} ||q_{l,k}[n]*F_{BB,l}[n]||^2
                 + sigma^2

#### F.4 Achievable Rate

The achievable downlink rate of UE k (in bits/s/Hz), accounting for pilot
overhead and cyclic prefix, is:

    R_k = zeta * (1/N) * sum_{n=0}^{N-1} log_2(1 + SINR_k[n])                              (36)

where zeta = (1 - tau_p/tau_c) * N/(N + N_cp) is the pre-log factor that
accounts for: (i) the uplink pilot training overhead (tau_p pilot symbols
out of tau_c coherence block symbols), and (ii) the cyclic prefix overhead
(N_cp guard samples per OFDM symbol).

The system sum-rate is:

    R_sum = sum_{k=1}^{K} R_k                                                                (37)


### G. Optimization Problem Formulation

We aim to jointly design the analog precoders {F_{RF,l}} and digital precoders
{F_{BB,l}[n]} to maximize the system sum-rate while ensuring a minimum quality
of service (QoS) for every user.

#### G.1 Two-Stage Optimization Problem

The joint sum-rate maximization with QoS guarantee is formulated as:

    (P1)  max_{ {F_{RF,l}}, {F_{BB,l}[n]} }    sum_{k=1}^{K} R_k

          subject to:

          (C1)  R_k  >=  R_min,    for all k = 1, ..., K
                [Minimum rate QoS constraint for every user]

          (C2)  (1/N) * sum_{n=0}^{N-1} ||F_{RF,l} * F_{BB,l}[n]||_F^2  <=  P_l,
                for all l = 1, ..., L
                [Per-AP average transmit power constraint]

          (C3)  |[F_{RF,l}]_{m,r}|  =  1/sqrt(M),
                for all m = 1,...,M,  r = 1,...,N_RF,  l = 1,...,L
                [Constant-modulus constraint on analog precoder]                              (38)

where R_min > 0 is the minimum rate threshold.

#### G.2 Analysis of Problem Difficulty

Problem (P1) is non-convex and highly challenging due to:

  1. Coupling between F_{RF,l} and F_{BB,l}[n]: The SINR and power constraint
     both involve the product F_{RF,l} * F_{BB,l}[n], creating bilinear coupling.

  2. Constant-modulus constraint (C3): The set of unit-modulus matrices is
     non-convex, making direct optimization intractable.

  3. Fractional SINR structure: The beamforming vectors appear in both the
     numerator and denominator of the SINR, rendering the objective non-concave.

  4. ICI dependence on beamformers: The ICI term (32) couples the beamforming
     design across ALL users and ALL subcarriers, increasing the dimensionality.

  5. QoS constraint (C1): This constraint is non-convex because R_k involves
     the logarithm of a ratio of beamforming-dependent quantities.

  6. Large-scale dimensionality: The total number of optimization variables
     is L*(M*N_RF + N*N_RF*K), which is very large for practical THz systems.

#### G.3 Equivalent Reformulation

To facilitate the solution, we reformulate (P1). First, note that for fixed
analog precoders {F_{RF,l}}, the effective channels {q_{l,k}[n]} in (29) are
determined, and (P1) reduces to a digital precoder design problem. This
motivates an alternating optimization approach.

Substituting the SINR expression (30) and rate expression (36) into (P1),
the problem can be equivalently written as:

    (P2)  max_{ {F_{RF,l}}, {f_{BB,l,k}[n]} }

          sum_{k=1}^{K} (zeta/N) * sum_{n=0}^{N-1} log_2( 1 +

              |sum_l q_{l,k}[n]*f_{BB,l,k}[n]|^2
              / ( sum_{j!=k} |sum_l q_{l,k}[n]*f_{BB,l,j}[n]|^2
                + epsilon*sum_l ||q_{l,k}[n]*F_{BB,l}[n]||^2 + sigma^2 ) )

          subject to:  (C1), (C2), (C3)                                                      (39)

#### G.4 Decomposition via Alternating Optimization

Problem (P2) can be decomposed into two subproblems that are solved
alternately until convergence:

  Subproblem 1 — Analog Precoder Design:

    For fixed {F_{BB,l}[n]}, optimize {F_{RF,l}} to maximize the sum-rate.
    The key challenge is the unit-modulus constraint (C3).
    Applicable methods: Manifold optimization, alternating minimization,
    codebook-based search.

  Subproblem 2 — Digital Precoder Design:

    For fixed {F_{RF,l}}, optimize {f_{BB,l,k}[n]} for all l, k, n.
    This is a multi-user MIMO beamforming problem over effective channels
    q_{l,k}[n], with the additional ICI term in the SINR denominator.
    Applicable methods: WMMSE (weighted minimum mean square error)
    transformation, fractional programming, successive convex approximation.


### H. Summary of Key Notation

| Symbol              | Dimension                  | Description                                      |
|---------------------|----------------------------|--------------------------------------------------|
| L                   | scalar                     | Number of APs                                    |
| M                   | scalar                     | Number of antennas per AP                        |
| N_RF                | scalar                     | Number of RF chains per AP                       |
| K                   | scalar                     | Number of single-antenna UEs                     |
| N                   | scalar                     | Number of OFDM subcarriers                       |
| B                   | scalar [Hz]                | Total bandwidth                                  |
| f_c                 | scalar [Hz]                | Center frequency (THz)                           |
| Delta_f             | scalar [Hz]                | Subcarrier spacing = B/N                         |
| f_n                 | scalar [Hz]                | Frequency of n-th subcarrier (Eq. 1)             |
| d_{l,k}             | scalar [m]                 | Distance from AP l to UE k                       |
| kappa_abs(f)        | scalar [m^{-1}]            | Molecular absorption coefficient                 |
| G_{l,k}(f_n)        | scalar                     | THz path gain (Eq. 5)                            |
| kappa_R             | scalar                     | Rician K-factor (linear)                         |
| h_{l,k}[n]          | C^{M x 1}                 | Channel vector from AP l to UE k on subcarrier n |
| a_M(phi, f_n)       | C^{M x 1}                 | Frequency-dependent array response vector (Eq. 9)|
| xi_n = f_n/f_c      | scalar                     | Beam squint ratio                                |
| F_{RF,l}            | C^{M x N_RF}              | Analog precoder at AP l (frequency-flat)         |
| F_{BB,l}[n]         | C^{N_RF x K}              | Digital precoder at AP l on subcarrier n         |
| f_{BB,l,k}[n]       | C^{N_RF x 1}              | Digital precoder: AP l, UE k, subcarrier n       |
| P_l                 | scalar [W]                 | Maximum transmit power of AP l                   |
| theta_l[t]          | scalar [rad]               | Phase noise at AP l (Wiener process)             |
| beta_{3dB,l}        | scalar [Hz]                | 3-dB linewidth of AP l's oscillator              |
| sigma_{theta,l}^2   | scalar [rad^2]             | Phase noise innovation variance (Eq. 16)         |
| J_{l,k}[p]          | scalar (complex)           | Phase noise DFT coefficient (Eq. 18)             |
| epsilon              | scalar                     | ICI power fraction (Eq. 23)                      |
| q_{l,k}[n]          | C^{1 x N_RF}              | Effective channel after analog BF (Eq. 29)       |
| SINR_k[n]           | scalar                     | SINR of UE k on subcarrier n (Eq. 30)            |
| R_k                 | scalar [bits/s/Hz]         | Achievable rate of UE k (Eq. 36)                 |
| R_min               | scalar [bits/s/Hz]         | Minimum QoS rate threshold                       |
| sigma^2             | scalar [W]                 | AWGN noise power                                 |


### I. Typical Parameter Values for THz CF-mMIMO Simulation

| Parameter                       | Symbol        | Typical Value              |
|---------------------------------|---------------|----------------------------|
| Center frequency                | f_c           | 100 - 300 GHz              |
| Total bandwidth                 | B             | 1 - 30 GHz                 |
| Number of subcarriers           | N             | 64 - 2048                  |
| Subcarrier spacing              | Delta_f       | 120 kHz - 480 kHz (5G NR)  |
| Fractional bandwidth            | beta = B/f_c  | 0.03 - 0.1                 |
| Number of APs                   | L             | 16 - 128                   |
| Antennas per AP                 | M             | 16 - 256                   |
| RF chains per AP                | N_RF          | 2 - 16                     |
| Number of UEs                   | K             | 4 - 32                     |
| Rician K-factor                 | kappa_R       | 5 - 20 dB                  |
| Number of NLoS clusters         | P             | 2 - 5                      |
| Path loss exponent              | -             | 1.8 - 2.2 (near free-space)|
| AP oscillator 3-dB linewidth    | beta_{3dB}^AP | 10 kHz - 1 MHz             |
| UE oscillator 3-dB linewidth    | beta_{3dB}^UE | 100 kHz - 10 MHz           |
| Per-AP transmit power           | P_l           | 20 - 30 dBm                |
| Noise figure                    | NF            | 7 - 10 dB                  |
| Coherence block length          | tau_c         | 200 - 500 symbols          |
| Cyclic prefix length            | N_cp          | N/4 (typical)              |


### J. Key References for Verification

[R1] Jornet & Akyildiz, "Channel Modeling and Capacity Analysis for THz Band,"
     IEEE Trans. Wireless Commun., 2011. (THz path loss model)

[R2] Kokkoniemi & Lehtomaki, "Simplified Molecular Absorption Loss Model for
     275-400 GHz," 2020. (Absorption coefficient approximation)

[R3] Tarboush et al., "TeraMIMO: A Channel Simulator for Wideband
     Ultra-Massive MIMO THz Communications," IEEE Trans. Veh. Tech., 2021.

[R4] arXiv:2405.09905, "Cell-Free THz Massive MIMO: A Novel Paradigm Beyond
     Ultra-Massive MIMO," 2024.

[R5] arXiv:2501.16306, "GNN-Based Hybrid Beamforming for Wideband THz
     MIMO-OFDM," 2025. (Beam squint model)

[R6] arXiv:2406.16303, "Hybrid Precoding with Low-Resolution Phase Shifters
     for Wideband THz with Beam Squint," 2024.

[R7] arXiv:2410.18722, "Uplink Cell-Free mMIMO OFDM with Phase Noise-Aware
     Channel Estimation," 2024. (Phase noise in CF-mMIMO OFDM)

[R8] arXiv:2405.04099, "Effect of Realistic Oscillator Phase Noise on
     Cell-Free Massive MIMO Systems," 2024.

[R9] arXiv:2305.12484, "Impact of Phase Noise on Uplink CF-mMIMO OFDM," 2023.

[R10] PLOS One 2025, "RNN-Based Hybrid Precoding for CF-mMIMO under THz," 2025.
