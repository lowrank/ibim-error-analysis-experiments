opt.type = 'semicircle'; % type of the domain
opt.R  = 0.75;       % radius of the circle
opt.random_rot = true;
opt.acc = 0.539370143939280;

opt.f =  @(x, y)( cos(x.^2 - y)) .* sin(y.^2 - x.^3) .* (y > 0); % test integrand 
fprintf('accurate integral is %1.15f by Mathematica.\n', opt.acc);

% regularity
opt.q = 1; % order of the regularity 

% tube width
opt.alpha = 1;  % parameter for the tube width
opt.beta  = 0.5 + (opt.q + 1) * (1 - opt.alpha);  % theoretical value of decay rate.
opt.upsilon = 1 + (opt.q + 1) * (1 - opt.alpha) / (opt.q + 2);

K = 24; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)