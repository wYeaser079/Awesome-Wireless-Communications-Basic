function [H, At, Ar, Alpha] = generate_mmWave_channel_UPA(Nt, Nr, Ncl, Nray, AS)
% GENERATE_MMWAVE_CHANNEL_UPA  Generate a narrowband Saleh-Valenzuela
%   clustered mmWave MIMO channel with UPA (Uniform Planar Array).
%
%   [H, At, Ar, Alpha] = generate_mmWave_channel_UPA(Nt, Nr, Ncl, Nray, AS)
%
%   Inputs:
%       Nt     - Number of transmit antennas (must be a perfect square for UPA)
%       Nr     - Number of receive antennas  (must be a perfect square for UPA)
%       Ncl    - Number of scattering clusters
%       Nray   - Number of rays per cluster
%       AS     - Angular spread in degrees (standard: 7.5)
%
%   Outputs:
%       H      - Nr x Nt channel matrix
%       At     - Nt x (Ncl*Nray) transmit array response matrix (dictionary)
%       Ar     - Nr x (Ncl*Nray) receive array response matrix  (dictionary)
%       Alpha  - 1 x (Ncl*Nray) complex path gains
%
%   Channel model:
%       H = sqrt(Nt*Nr/(Ncl*Nray)) * sum_{i,l} alpha_{il} * a_r(...) * a_t(...)^H
%
%   Ray angles are generated from a Laplacian distribution around cluster
%   mean angles.
%
%   Reference:
%       El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., vol. 13, no. 3,
%       pp. 1499-1513, Mar. 2014.
%
%   Source: Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%          metegenez/spatially-sparse-precoding
%
% =========================================================================

Npath = Ncl * Nray;

% --- Generate random cluster mean angles ---
% AoD/EoD: Angle/Elevation of Departure
% AoA/EoA: Angle/Elevation of Arrival
minAOD = -30;  maxAOD = 30;
minEOD = 80;   maxEOD = 100;

Cluster_AOD = rand(Ncl, 1) * (maxAOD - minAOD) + minAOD;
Cluster_EOD = rand(Ncl, 1) * (maxEOD - minEOD) + minEOD;
Cluster_AOA = (rand(Ncl, 1) - 0.5) * 360;
Cluster_EOA = rand(Ncl, 1) * 180;

% --- Generate ray angles with Laplacian distribution ---
b = AS / sqrt(2);  % Laplacian scale parameter
Randomness = rand(Npath, 1) - 0.5;

Ray_AOD = repelem(Cluster_AOD, Nray, 1) - b * sign(Randomness) .* log(1 - 2 * abs(Randomness));
Ray_EOD = repelem(Cluster_EOD, Nray, 1) - b * sign(Randomness) .* log(1 - 2 * abs(Randomness));
Ray_AOA = repelem(Cluster_AOA, Nray, 1) - b * sign(Randomness) .* log(1 - 2 * abs(Randomness));
Ray_EOA = repelem(Cluster_EOA, Nray, 1) - b * sign(Randomness) .* log(1 - 2 * abs(Randomness));

% --- Antenna element positions for UPA (y-z plane) ---
Nt_H = sqrt(Nt);  Nt_V = sqrt(Nt);
X_Tx = zeros(1, Nt);
[Y_Tx, Z_Tx] = meshgrid(0:Nt_H-1, 0:Nt_V-1);
TxPos = [X_Tx; Y_Tx(:).'; Z_Tx(:).'];

Nr_H = sqrt(Nr);  Nr_V = sqrt(Nr);
X_Rx = zeros(1, Nr);
[Y_Rx, Z_Rx] = meshgrid(0:Nr_H-1, 0:Nr_V-1);
RxPos = [X_Rx; Y_Rx(:).'; Z_Rx(:).'];

% --- Spherical unit vectors ---
SphericalUnitVecTx = get_spherical_unit_vector(Ray_EOD, Ray_AOD);
SphericalUnitVecRx = get_spherical_unit_vector(Ray_EOA, Ray_AOA);

% --- Array response matrices ---
At = (1/sqrt(Nt)) * exp(1j * pi * TxPos.' * SphericalUnitVecTx);  % Nt x Npath
Ar = (1/sqrt(Nr)) * exp(1j * pi * RxPos.' * SphericalUnitVecRx);  % Nr x Npath

% --- Complex path gains (CN(0,1)) ---
Alpha = sqrt(1/2) * (randn(1, Npath) + 1j * randn(1, Npath));

% --- Generate channel matrix ---
H = sqrt(Nt * Nr / Npath) * Ar * diag(Alpha) * At';

end

% =========================================================================
function SphericalUnitVector = get_spherical_unit_vector(theta, phi)
% theta and phi in DEGREES
SphericalUnitVector = [sind(theta) .* cosd(phi), ...
                       sind(theta) .* sind(phi), ...
                       cosd(theta)]';  % 3 x Npath
end
