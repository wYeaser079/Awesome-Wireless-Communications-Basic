function a = steering_vec_UPA_3D(theta, phi, Nh, Nv, d)
% STEERING_VEC_UPA_3D  Compute UPA steering vector using the spherical unit
%   vector approach (as in metegenez/spatially-sparse-precoding).
%
%   a = steering_vec_UPA_3D(theta, phi, Nh, Nv)
%   a = steering_vec_UPA_3D(theta, phi, Nh, Nv, d)
%
%   Inputs:
%       theta - Elevation angle in DEGREES (from z-axis)
%       phi   - Azimuth angle in DEGREES
%       Nh    - Number of horizontal elements
%       Nv    - Number of vertical elements
%       d     - Antenna spacing in wavelengths (default: 0.5)
%
%   Output:
%       a     - (Nh*Nv) x 1 steering vector
%
%   This uses the 3D spherical unit vector:
%       e = [sin(theta)*cos(phi); sin(theta)*sin(phi); cos(theta)]
%   and antenna positions on a planar grid in the y-z plane.
%
%   Reference:
%       Naren920421/Narrowband-mmWave-hybrid-precoding-algorithm
%       metegenez/spatially-sparse-precoding (ChannelGeneration.m)
%
% =========================================================================

if nargin < 5
    d = 0.5;
end

N = Nh * Nv;

% Spherical unit vector (using degrees)
e = [sind(theta) * cosd(phi);
     sind(theta) * sind(phi);
     cosd(theta)];

% Antenna element position vectors (y-z plane, x = 0)
X = zeros(1, N);
[Y, Z] = meshgrid(0:Nh-1, 0:Nv-1);
pos = [X; Y(:).'; Z(:).'];  % 3 x N position matrix

% Steering vector: a = (1/sqrt(N)) * exp(j * 2*pi*d * pos^T * e)
% For d = 0.5, 2*pi*d = pi
a = (1/sqrt(N)) * exp(1j * 2 * pi * d * pos.' * e);

end
