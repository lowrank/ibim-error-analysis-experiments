opt.type = 'sphere';
opt.R = 0.75;
opt.f = @(x,y,z)(( cos(x.^2 - y - z.^3)) .* sin(y.^2 - x.^3 - z)); 

% accurate integral
opt.acc = 1.71628896296959  * 0.75^2;

% regularity of weight function
opt.q = 1;
    
% tube width
opt.alpha = 1;
opt.beta = 1 + (opt.q + 1)* (1 - opt.alpha);

if isfield(opt, 'upsilon')
    opt = rmfield(opt, 'upsilon');
end

K = 16;
S = 32;

ibim_3d_experiments(opt, K, S)