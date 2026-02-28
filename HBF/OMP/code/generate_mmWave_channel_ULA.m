function [H, At, Ar, Alpha, Phi_AOD, Phi_AOA] = generate_mmWave_channel_ULA(Nt, Nr, Ncl, Nray, sigma_AS, d)
% GENERATE_MMWAVE_CHANNEL_ULA  Generate a narrowband Saleh-Valenzuela
%   clustered mmWave MIMO channel with ULA (Uniform Linear Array).
%
%   [H, At, Ar, Alpha, Phi_AOD, Phi_AOA] = generate_mmWave_channel_ULA(Nt, Nr, Ncl, Nray, sigma_AS)
%   [H, At, Ar, Alpha, Phi_AOD, Phi_AOA] = generate_mmWave_channel_ULA(Nt, Nr, Ncl, Nray, sigma_AS, d)
%
%   Inputs:
%       Nt       - Number of transmit antennas
%       Nr       - Number of receive antennas
%       Ncl      - Number of scattering clusters
%       Nray     - Number of rays per cluster
%       sigma_AS - Angular spread (radians) for Laplacian distribution
%       d        - Antenna spacing in wavelengths (default: 0.5)
%
%   Outputs:
%       H        - Nr x Nt channel matrix
%       At       - Nt x (Ncl*Nray) transmit steering vector matrix
%       Ar       - Nr x (Ncl*Nray) receive steering vector matrix
%       Alpha    - Ncl x Nray complex path gains
%       Phi_AOD  - Ncl x Nray angles of departure (radians)
%       Phi_AOA  - Ncl x Nray angles of arrival (radians)
%
%   Channel model:
%       H = sqrt(Nt*Nr/(Ncl*Nray)) * sum_{i,l} alpha_{il} * a_r(phi_AOA) * a_t(phi_AOD)^H
%
%   ULA steering vector: a(phi) = (1/sqrt(N)) * [1, e^{-j*pi*sin(phi)}, ..., e^{-j*pi*(N-1)*sin(phi)}]^T
%
%   Reference:
%       ge99210/Hybrid-Precoding-Combining- (mm_wave_channel_v3_2D.m)
%       Georgios K. Papageorgiou, 2019
%
% =========================================================================

if nargin < 6
    d = 0.5;
end

Npath = Ncl * Nray;

% --- ULA steering vector function ---
a = @(phi, N) exp(-1j * 2 * pi * d * sin(phi) * (0:N-1)).' / sqrt(N);

% --- Generate random cluster mean angles (uniform over [-pi, pi]) ---
Phi_AOD_mean = -pi + 2 * pi * rand(Ncl, 1);
Phi_AOA_mean = -pi + 2 * pi * rand(Ncl, 1);

% --- Complex path gains (CN(0,1)) ---
Alpha = (1/sqrt(2)) * (randn(Ncl, Nray) + 1j * randn(Ncl, Nray));

% --- Generate ray angles with truncated Laplacian distribution ---
Phi_AOD = zeros(Ncl, Nray);
Phi_AOA = zeros(Ncl, Nray);
H = zeros(Nr, Nt);

for i = 1:Ncl
    phi_AOD = truncated_laplacian_rnd(Nray, 1, Phi_AOD_mean(i), sigma_AS);
    phi_AOA = truncated_laplacian_rnd(Nray, 1, Phi_AOA_mean(i), sigma_AS);
    for l = 1:Nray
        Phi_AOD(i, l) = phi_AOD(l);
        Phi_AOA(i, l) = phi_AOA(l);
        alpha_il = Alpha(i, l);
        A_t = a(phi_AOD(l), Nt);
        A_r = a(phi_AOA(l), Nr);
        H = H + alpha_il * A_r * A_t';
    end
end

H = sqrt(Nt * Nr / Npath) * H;

% --- Also return full steering vector matrices for dictionary use ---
At = zeros(Nt, Npath);
Ar = zeros(Nr, Npath);
idx = 0;
for i = 1:Ncl
    for l = 1:Nray
        idx = idx + 1;
        At(:, idx) = a(Phi_AOD(i, l), Nt);
        Ar(:, idx) = a(Phi_AOA(i, l), Nr);
    end
end

end

% =========================================================================
function y = truncated_laplacian_rnd(m, n, mu, sigma)
% Truncated Laplacian random number generator.
%
% Reference: ge99210/Hybrid-Precoding-Combining- (mm_wave_channel_v3_2D.m)

if nargin < 4
    sigma = 1;
end
if nargin < 3
    mu = 0;
end

u = rand(m, n) - 0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u) .* log(1 - 2 * (1 - exp(-pi/b)) * abs(u));

end
