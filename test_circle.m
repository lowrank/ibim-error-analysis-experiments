opt.type = 'circle'; % type of the domain
opt.R  = 0.75;       % radius of the circle
opt.random_rot = false;
opt.acc = 0.992902336403640;

% regularity
opt.q = 1; % order of the regularity 

opt.f =  @(x, y)( cos(x.^2 - y)) .* sin(y.^2 - x.^3); % test integrand 

fprintf('accurate integral is %1.15f by Mathematica.\n', opt.acc);

% tube width
opt.alpha = 1;  % parameter for the tube width
opt.beta  = 0.5 + (opt.q + 1) * (1 - opt.alpha);  % theoretical value of decay rate.

K = 1; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)