%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% Capsule Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt.type = 'capsule';
opt.q = 1;

% integrand function
opt.f = @(x, y) ( cos(x.^2 - y)) .* sin(y.^2 - x.^3);

% In the implementation, instead of rotating the geometry, we rotate the
% lattice points for keep the integrand function unmodified.

% inclination of capsule.
opt.slope = sqrt(2);  

% radius of the semi-circles in capsule.
%
% This automatically cuts the domain into 3 parts based on x-values: 
% [-\inf, -0.5], [-0.5, 0.5], [0.5, \inf].
%
opt.R = 0.2;

% compute accurate integral.
acc_semi_circles = 0.02272974686046358;
acc_rectangle    = 0.0778101035506957;
acc = acc_semi_circles + acc_rectangle;

%% tube width
alpha= 1;
beta = 2-alpha; % for lines
kappa = 0.5 + (opt.q + 1) * (1 - alpha); % for curves
delta = 3 - alpha;% for lines, random

K = 24; % number of grid sizes
S = 32; % number of sampled rigid transforms

ret = zeros(K, S);

base_grid = 100;
grow_rate = 1.2;

if S == 1
    opt.random = false;
    v  = [opt.slope,  1]/sqrt(1 + opt.slope^2);

    parfor k = 1:K
        N = floor(grow_rate^k * base_grid);
        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * (2/N)^alpha;
        end

        pts = linspace(-1, 1, N + 1);    
        h = 2 / N;   
        [X, Y] = meshgrid(pts);
        
        M = (N + 1)^2;
    
        for i = 1:M
            rx = X(i) * v(1) + Y(i) * v(2);
            ry = X(i) * v(2) - Y(i) * v(1);
            
            if rx > 0.5
                % in right semi-cicle centered at (0.5, 0),
                dist   =  sqrt((rx -0.5)^2 + ry^2) - opt.R;
                normal = [rx - 0.5, ry] / sqrt((rx-0.5)^2 + ry^2);
                Jac    = opt.R / (dist + opt.R);
                px = rx - dist * normal(1);
                py = ry - dist * normal(2);
                if abs(dist) <= EPS
                    ret(k) = ret(k) + opt.f(px, py) *...
                        weight_func(dist, EPS, opt) * Jac;
                end

            elseif rx < -0.5
                % in left semi-circle centered at (-0.5, 0)
                dist   =  sqrt((rx + 0.5)^2 + ry^2) - opt.R;
                normal = [rx + 0.5, ry] / sqrt((rx + 0.5)^2 + ry^2);
                Jac    = opt.R / (dist + opt.R);
                px = rx - dist * normal(1);
                py = ry - dist * normal(2);
                if abs(dist) <= EPS
                    ret(k) = ret(k) + opt.f(px, py) *...
                        weight_func(dist, EPS, opt) * Jac;
                end
            else
                % in rectangle
                if abs(ry - opt.R) <= EPS
                    ret(k) = ret(k) +  opt.f(rx, opt.R) *...
                        weight_func(ry - opt.R, EPS, opt);
                elseif abs(ry + opt.R) <= EPS
                    ret(k) = ret(k) +  opt.f(rx, -opt.R) *...
                        weight_func(ry + opt.R, EPS, opt);
                end
            end
        end
        ret(k) = ret(k) * h^2;
    end
    %% plot the convergence rate
    g = 2.0./(floor(grow_rate.^(1:K) * base_grid));
    
    err = abs(ret - acc);
    gamma = ( g.^(beta) * err) / norm(g.^(beta))^2;
    eta = ( g.^(delta) * err) / norm(g.^(delta))^2;
    loglog(g, err, '-bo',  g, gamma * g.^(beta), '-r',...
        g, eta * g.^(delta), '--k');
    
    legend_handler = legend('quadrature error',...
        sprintf('O(h^{%1.1f})', beta),...
        sprintf('O(h^{%1.1f})', delta),...
        'Location', 'southeast');
    fontsize(legend_handler,18,'points');
    
    grid on;
else
    for s = 1:S
        % random rotation
   
        theta = rand() * 2 * pi;
        v = [cos(theta), sin(theta)];

        parfor k = 1:K
            N = floor(grow_rate^k * base_grid);
            if alpha == 0
                EPS = 0.1;
                shift = EPS/10 * [rand(), rand()];
            else
                EPS = 2 * (2/N)^alpha;
                shift = EPS/2 * [rand(), rand()];
            end
    
            pts = linspace(-1, 1, N + 1);    
            h = 2 / N;   
            [X, Y] = meshgrid(pts);

            X = X + shift(1);
            Y = Y + shift(2);
            
            M = (N + 1)^2;
        
            for i = 1:M
                rx = X(i) * v(1) + Y(i) * v(2);
                ry = X(i) * v(2) - Y(i) * v(1);
                
                if rx > 0.5
                    % in right semi-cicle centered at (0.5, 0),
                    dist   =  sqrt((rx -0.5)^2 + ry^2) - opt.R;
                    normal = [rx - 0.5, ry] / sqrt((rx-0.5)^2 + ry^2);
                    Jac    = opt.R / (dist + opt.R);
                    px = rx - dist * normal(1);
                    py = ry - dist * normal(2);
                    if abs(dist) <= EPS
                        ret(k,s) = ret(k,s) + opt.f(px, py) *...
                            weight_func(dist, EPS, opt) * Jac;
                    end
    
                elseif rx < -0.5
                    % in left semi-circle centered at (-0.5, 0)
                    dist   =  sqrt((rx + 0.5)^2 + ry^2) - opt.R;
                    normal = [rx + 0.5, ry] / sqrt((rx + 0.5)^2 + ry^2);
                    Jac    = opt.R / (dist + opt.R);
                    px = rx - dist * normal(1);
                    py = ry - dist * normal(2);
                    if abs(dist) <= EPS
                        ret(k,s) = ret(k,s) + opt.f(px, py) *...
                            weight_func(dist, EPS, opt) * Jac;
                    end
                else
                    % in rectangle
                    if abs(ry - opt.R) <= EPS
                        ret(k,s) = ret(k,s) +  opt.f(rx, opt.R) *...
                            weight_func(ry - opt.R, EPS, opt);
                    elseif abs(ry + opt.R) <= EPS
                        ret(k,s) = ret(k,s) +  opt.f(rx, -opt.R) *...
                            weight_func(ry + opt.R, EPS, opt);
                    end
                end
            end
            ret(k,s) = ret(k,s) * h^2;
        end
    end

    %% plot the convergence rate
    g = 2./( floor(base_grid * grow_rate.^(1:K)));

    var_err = sum( (ret - acc).^2, 2)/S;
    gamma = ( g.^(delta) * var_err) / norm(g.^(delta))^2;% line

    eta = ( g.^(2*kappa) * var_err) / norm(g.^(2*kappa))^2; %curve

    loglog(g, var_err, '-bo',  g, gamma * g.^(delta), '-r', g, eta * g.^(2*kappa), '--k');
    legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', delta),...
        sprintf('O(h^{%1.1f})', 2*kappa), 'Location', 'southeast');
    fontsize(legend_handler,18,'points');

    grid on;
end
