function a = steering_vec_ULA(phi, N, d)
% STEERING_VEC_ULA  Compute ULA steering vector for a given angle.
%
%   a = steering_vec_ULA(phi, N)
%   a = steering_vec_ULA(phi, N, d)
%
%   Inputs:
%       phi  - Angle of departure/arrival in RADIANS
%       N    - Number of antenna elements
%       d    - Antenna spacing in wavelengths (default: 0.5)
%
%   Output:
%       a    - N x 1 steering vector (normalized by 1/sqrt(N))
%
%   The ULA steering vector for half-wavelength spacing is:
%       a(phi) = (1/sqrt(N)) * [1, e^{j*pi*sin(phi)}, ..., e^{j*pi*(N-1)*sin(phi)}]^T
%
%   Reference:
%       El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., vol. 13, no. 3,
%       pp. 1499-1513, Mar. 2014.
%
% =========================================================================

if nargin < 3
    d = 0.5;  % half-wavelength spacing by default
end

n = (0:N-1).';  % column vector of element indices
a = (1/sqrt(N)) * exp(1j * 2 * pi * d * n * sin(phi));

end
