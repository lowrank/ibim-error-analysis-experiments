opt.type = 'line'; % type of the domain
opt.random_rot = true;

opt.acc = 0.25;

% integrand function
opt.f =  @(x, y)(3*(15*x+x^2)).* (x > -0.5) .* (x < 0.5); % test integrand 
% regularity
opt.q = 2;

% tube width
opt.alpha= 1;
opt.beta = (3 - opt.alpha)/2; % variance of random slope
opt.upsilon = 0.5;

K = 24; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)