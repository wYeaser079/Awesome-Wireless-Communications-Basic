function [A_dict, angles] = build_dictionary_ULA(N, G, d)
% BUILD_DICTIONARY_ULA  Construct a dictionary matrix of ULA steering vectors
%   sampled on a uniform angular grid.
%
%   [A_dict, angles] = build_dictionary_ULA(N, G)
%   [A_dict, angles] = build_dictionary_ULA(N, G, d)
%
%   Inputs:
%       N      - Number of antenna elements
%       G      - Number of grid points (dictionary size / resolution)
%                Typically G >= N (oversampled dictionary when G > N)
%       d      - Antenna spacing in wavelengths (default: 0.5)
%
%   Outputs:
%       A_dict - N x G dictionary matrix (each column is a steering vector)
%       angles - 1 x G vector of quantized angles in radians
%
%   The dictionary samples the angular space uniformly. Two common
%   parameterizations are available:
%
%   Method 1 (Angle-domain): Uniformly sample sin(phi) in [-1, 1)
%       This gives a DFT-like dictionary when G = N.
%
%   Method 2 (Physical angle): Uniformly sample phi in [-pi/2, pi/2)
%
%   Here we use Method 1 (sin-domain), which is standard in the OMP
%   hybrid precoding literature. When G = N, this reduces to a DFT matrix.
%
%   Reference:
%       El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., 2014.
%
% =========================================================================

if nargin < 3
    d = 0.5;
end

% Uniform grid in the spatial frequency domain: sin(phi) in [-1, 1)
% This is the standard parameterization for the OMP dictionary
sin_grid = -1 + (0:G-1) * (2/G);  % G points uniformly in [-1, 1)
angles = asin(sin_grid);           % corresponding physical angles

% Build dictionary: each column is a(sin_val) = (1/sqrt(N)) * exp(j*2*pi*d*n*sin_val)
n = (0:N-1).';   % column vector of antenna indices
A_dict = (1/sqrt(N)) * exp(1j * 2 * pi * d * n * sin_grid);   % N x G

end
