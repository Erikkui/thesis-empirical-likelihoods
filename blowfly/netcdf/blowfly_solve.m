function [N, N_burned] = blowfly_solve( theta, t, varargin )
    % Optional arguments
    p = inputParser;
    addParameter(p, 'N_init', 180);
    addParameter(p, 'burn_in', 20);
    addParameter(p, 'mu', 1.0);
    parse(p, varargin{:});
    N_init = p.Results.N_init;
    burn_in = p.Results.burn_in;
    mu = p.Results.mu;

    % Unpack parameters
    delta     = theta(1);
    P         = theta(2);
    N_0       = theta(3);
    sigma2_p  = theta(4);
    tau       = theta(5);
    sigma2_d  = theta(6);

    % Initialize variables
    lag = floor(tau);
    if lag < 1
        lag = 1;
    end
    total_time = t + lag + burn_in;
    N = zeros(1, total_time);
    N(1:lag) = N_init;

    % Gamma distributions for noise
    shape_p = mu^2 / sigma2_p;
    scale_p = sigma2_p / mu;
    shape_d = mu^2 / sigma2_d;
    scale_d = sigma2_d / mu;

    ee = gamrnd(shape_p, scale_p, 1, total_time);    % Reproduction noise
    epsilon = gamrnd(shape_d, scale_d, 1, total_time); % Death rate noise

    % Time loop
    for ii = (lag+1):total_time
        Nlag = N(ii - lag);          % Delayed population
        Nprev = N(ii - 1);           % Previous time step
        ee_t = ee(ii - 1);           % Reproduction noise
        epsilon_t = epsilon(ii - 1); % Death noise

        R_t = P * Nlag * exp(-Nlag / N_0) * ee_t;
        S_t = Nprev * exp(-delta * epsilon_t);
        N(ii) = R_t + S_t;
    end

    % Remove burn-in
    N_burned = N((1 + lag + burn_in):end);
end
