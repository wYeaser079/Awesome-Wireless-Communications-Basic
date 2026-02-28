function a = steering_vec_UPA(theta, phi, Nh, Nv, d)
% STEERING_VEC_UPA  Compute UPA (Uniform Planar Array) steering vector.
%
%   a = steering_vec_UPA(theta, phi, Nh, Nv)
%   a = steering_vec_UPA(theta, phi, Nh, Nv, d)
%
%   Inputs:
%       theta - Elevation angle in RADIANS (from z-axis, 0 to pi)
%       phi   - Azimuth angle in RADIANS   (-pi to pi)
%       Nh    - Number of horizontal antenna elements
%       Nv    - Number of vertical antenna elements
%       d     - Antenna spacing in wavelengths (default: 0.5)
%
%   Output:
%       a     - (Nh*Nv) x 1 steering vector (normalized by 1/sqrt(Nh*Nv))
%
%   The UPA lies in the y-z plane. The total number of antennas is N = Nh * Nv.
%   For half-wavelength spacing:
%       a(theta, phi) = (1/sqrt(N)) * kron(a_v(theta), a_h(theta, phi))
%
%   where:
%       a_h(m) = exp(j * 2*pi*d * m * sin(theta)*sin(phi)),  m = 0,...,Nh-1
%       a_v(n) = exp(j * 2*pi*d * n * cos(theta)),           n = 0,...,Nv-1
%
%   Reference:
%       El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., vol. 13, no. 3,
%       pp. 1499-1513, Mar. 2014.
%
% =========================================================================

if nargin < 5
    d = 0.5;
end

N = Nh * Nv;  % total number of elements

% Horizontal (azimuth) array response
m = (0:Nh-1).';
a_h = exp(1j * 2 * pi * d * m * sin(theta) * sin(phi));

% Vertical (elevation) array response
n = (0:Nv-1).';
a_v = exp(1j * 2 * pi * d * n * cos(theta));

% UPA steering vector = Kronecker product (normalized)
a = (1/sqrt(N)) * kron(a_v, a_h);

end
