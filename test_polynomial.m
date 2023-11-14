opt.type = 'polynomial'; % type of the domain
opt.random_rot = true;
opt.pts = linspace(-0.75, 0.75, 100);
opt.anchor = [opt.pts ; (opt.pts).^6]';
opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

opt.acc = 0.994644710747744;

% integrand function
opt.f =  @(x, y)( cos(x.^2 - y)) .* cos(y.^2 - x.^3) .* (x > -0.5) .* (x < 0.5); % test integrand 
% regularity
opt.q = 1;

% tube width
opt.alpha = 1;  % parameter for the tube width
opt.beta  = 0.5 + (1 - opt.alpha);  % theoretical value of decay rate.

K = 24; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)