
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tests for Jacobian evaluation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../'));

opt.type = 'degenerate'; % type of the domain
opt.random_rot = false;
opt.Rx  = 0.75^2;       
opt.Ry  = 0.75;
angles = linspace(0, 2*pi, 100);
opt.angles = angles(1:end-1);
opt.anchor = [sqrt( opt.Rx * abs(cos(opt.angles)) ) .* sign(cos(opt.angles)) ; opt.Ry * sin(opt.angles)]';
opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

% opt.type = 'ellipse'; % type of the domain
% opt.Rx = 0.75;
% opt.Ry = 0.6;
% angles = linspace(0, 2*pi, 100);
% opt.angles = angles(1:end-1);
% opt.anchor = [opt.Rx * cos(opt.angles); opt.Ry * sin(opt.angles)]';
% opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

% opt.type = 'star';
% opt.random_rot = true;
% opt.R = 0.75;
% opt.r = 0.2;
% opt.m = 3;
% angles = linspace(0, 2*pi, 100);
% opt.angles = angles(1:end-1);
% opt.anchor = [ (opt.R + opt.r * cos(opt.m * opt.angles) ) .* cos(opt.angles) ;...
%     (opt.R + opt.r * cos(opt.m * opt.angles) ) .* sin(opt.angles)]';

x = [0.95, 0.2];
[dist0, normal0, jac0] = dist2curve(x(1), x(2), opt);

h = 1e-4;

dist1 = dist2curve(x(1)+h, x(2), opt);
dist2 = dist2curve(x(1), x(2)+h, opt);
dist3 = dist2curve(x(1)-h, x(2), opt);
dist4 = dist2curve(x(1), x(2)-h, opt);

jac_h = 1 - dist0*(dist1 + dist2 + dist3 + dist4 - 4 * dist0)/(h^2);

fprintf( 'jacobian error using discrete laplacian %6.4e\n', jac_h - jac0);
