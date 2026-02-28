function [A_dict, theta_grid, phi_grid] = build_dictionary_UPA(Nh, Nv, Gh, Gv)
% BUILD_DICTIONARY_UPA  Construct a dictionary matrix of UPA steering vectors
%   sampled on a uniform 2D angular grid.
%
%   [A_dict, theta_grid, phi_grid] = build_dictionary_UPA(Nh, Nv, Gh, Gv)
%
%   Inputs:
%       Nh         - Number of horizontal antenna elements
%       Nv         - Number of vertical antenna elements
%       Gh         - Number of azimuth grid points
%       Gv         - Number of elevation grid points
%
%   Outputs:
%       A_dict     - (Nh*Nv) x (Gh*Gv) dictionary matrix
%       theta_grid - 1 x Gv elevation angles in degrees
%       phi_grid   - 1 x Gh azimuth angles in degrees
%
%   The dictionary is formed by evaluating the UPA steering vector
%   on a uniform grid of (elevation, azimuth) angles. The total
%   dictionary size is G = Gh * Gv.
%
%   Reference:
%       El Ayach et al., "Spatially Sparse Precoding in Millimeter Wave
%       MIMO Systems," IEEE Trans. Wireless Commun., 2014.
%
% =========================================================================

N = Nh * Nv;
G = Gh * Gv;

% Define angular grids
theta_grid = linspace(0, 180, Gv);     % elevation: 0 to 180 degrees
phi_grid   = linspace(-180, 180, Gh);  % azimuth: -180 to 180 degrees

% Antenna positions on y-z plane
X = zeros(1, N);
[Y, Z] = meshgrid(0:Nh-1, 0:Nv-1);
pos = [X; Y(:).'; Z(:).'];  % 3 x N

% Build dictionary
A_dict = zeros(N, G);
idx = 0;
for iv = 1:Gv
    for ih = 1:Gh
        idx = idx + 1;
        % Spherical unit vector
        e = [sind(theta_grid(iv)) * cosd(phi_grid(ih));
             sind(theta_grid(iv)) * sind(phi_grid(ih));
             cosd(theta_grid(iv))];
        % Steering vector (half-wavelength spacing -> 2*pi*d = pi)
        A_dict(:, idx) = (1/sqrt(N)) * exp(1j * pi * pos.' * e);
    end
end

end
