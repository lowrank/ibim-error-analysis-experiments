opt.type = 'ellipse'; % type of the domain
opt.Rx = 0.75;
opt.Ry = 0.6;
angles = linspace(0, 2*pi, 100);
opt.angles = angles(1:end-1);
opt.anchor = [opt.Rx * cos(opt.angles); opt.Ry * sin(opt.angles)]';
opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);
opt.acc = 0.672602390335232;

% integrand function
opt.f =  @(x, y)( cos(x.^2 - y)) .* sin(y.^2 - x.^3); % test integrand 

% tube width
opt.alpha = 1;  % parameter for the tube width
opt.beta  = 0.5 + (opt.q + 1) * (1 - opt.alpha);  % theoretical value of decay rate.

K = 24; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)