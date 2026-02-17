function p = params()
%PARAMS Default simulation parameters for Zhang et al. (JCN, 2019).
%   "Cell-Free Massive MIMO: Zero Forcing and Conjugate Beamforming Receivers"
%
%   Usage: p = params();
%
%   All parameters match the paper's Section IV (Numerical Results).

    % --- Network geometry ---
    p.area_side = 1000;         % Area side length [m] (1 km x 1 km square)

    % --- Network dimensions ---
    p.M   = 60;                 % Number of Access Points (APs)
    p.K   = 10;                 % Number of single-antenna users
    p.N   = 6;                  % Number of antennas per AP (paper tests N=2,4,6)

    % --- Coherence interval ---
    p.T   = 200;                % Coherence interval length [samples]
    p.tau = 10;                 % Pilot sequence length [samples]
    p.prelog = 1 - p.tau/p.T;   % Pre-log factor for spectral efficiency

    % --- Transmit powers (Section IV) ---
    % Paper: rho_p = rho_u = 200 mW, noise power = -96 dBm (BW = 20 MHz)
    p.pilot_power_mW = 200;     % Pilot transmit power [mW]
    p.data_power_mW  = 200;     % Uplink data transmit power [mW]
    p.noise_dBm      = -96;     % Noise power [dBm] (BW = 20 MHz)

    % Normalized SNR = transmit_power / noise_power (linear scale)
    noise_mW  = 10^(p.noise_dBm / 10);       % Noise in mW
    p.rho_p   = p.pilot_power_mW / noise_mW;  % Normalized pilot SNR
    p.rho_u   = p.data_power_mW  / noise_mW;  % Normalized data SNR

    % --- Path loss model (Ngo et al. 2017, Eq. 46-47) ---
    % Three-slope model parameters
    p.d0 = 10;                  % Reference distance [m]
    p.d1 = 50;                  % Breakpoint distance [m]
    p.hAP = 15;                 % AP antenna height [m]
    p.hUE = 1.65;               % User antenna height [m]
    p.sigma_sf = 8;             % Shadow fading std dev [dB]
    p.freq_GHz = 1.9;           % Carrier frequency [GHz]

    % Noise figure
    p.noise_figure_dB = 9;      % Noise figure [dB]
    p.bandwidth_MHz   = 20;     % System bandwidth [MHz]

    % --- Monte Carlo ---
    p.num_setups       = 100;   % Number of random network realizations
    p.num_channel_real = 300;   % Channel realizations per setup

end
