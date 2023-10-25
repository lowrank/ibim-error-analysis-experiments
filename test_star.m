opt.type = 'star';
opt.random_rot = true;
opt.R = 0.75;
opt.r = 0.2;
opt.m = 3;
angles = linspace(0, 2*pi, 100);
opt.angles = angles(1:end-1);
opt.anchor = [ (opt.R + opt.r * cos(opt.m * opt.angles) ) .* cos(opt.angles) ;...
    (opt.R + opt.r * cos(opt.m * opt.angles) ) .* sin(opt.angles)]';
opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

opt.f =  @(x, y)( cos(x.^2 - y)) .* sin(y.^2 - x.^3); % test integrand 
opt.acc = 0.986770621149293;
% regularity
opt.q = 1;
% tube width
opt.alpha = 1;  % parameter for the tube width
opt.beta  = 0.5 + (opt.q + 1) * (1 - opt.alpha);  % theoretical value of decay rate.

K = 10; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)