% uncomment the following code block for experiments on convex curve with '
% zero curvature points.
opt.type = 'degenerate'; % type of the domain
opt.random_rot = true;
opt.Rx  = 0.75^2;       
opt.Ry  = 0.75;
angles = linspace(0, 2*pi, 100);
opt.angles = angles(1:end-1);
opt.anchor = [sqrt( opt.Rx * abs(cos(opt.angles)) ) .* sign(cos(opt.angles)) ; opt.Ry * sin(opt.angles)]';
opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);
opt.acc = 1.15335526133127;

% integrand function
opt.f =  @(x, y)( cos(x.^2 - y)) .* sin(y.^2 - x.^3); % test integrand 
% regularity
opt.q = 2;
% tube width
opt.alpha = 0;  % parameter for the tube width
opt.beta  = 0.5 + (opt.q + 1) * (1 - opt.alpha);  % theoretical value of decay rate.

if isfield(opt, 'upsilon')
    opt = rmfield(opt, 'upsilon');
end


K = 24; % number of grid sizes.
S = 32;  % number of sampled translations.

ibim_2d_experiments(K, S, opt)