%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% Line Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('./'));

opt.type = 'line';
opt.q = 2;

% integrand is defined along the "line".
opt.f = @(x) (3*x^2);

% quadratically irrational slope, used only if opt.random = false
opt.slope = sqrt(2);  
% opt.slope = (sqrt(5)+1)/2;

% accurate integral over [-0.5, 0.5]
acc = 0.25;

%% tube width
alpha= 0;
beta = 2 - alpha; % irrational slope
delta = 3 - alpha; % variance of random slope

K = 24; % number of grid sizes
S = 64; % number of sampled rigid transforms

ret = zeros(K*S, 1);

base_grid = 100;
grow_rate = 1.2;

if S == 1
    opt.random = false;
    v  = [opt.slope,  1]/sqrt(1 + opt.slope^2);
    progress = PoolWaitbar(K, 'Starting');
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

        Data = zeros(M , 1);
    
        for i = 1:M
            rx = X(i) * v(1) + Y(i) * v(2);
            ry = X(i) * v(2) - Y(i) * v(1);
            
            if rx >= -0.5 && rx <= 0.5
                if abs(ry) <= EPS
                    ret(k) = ret(k) +  opt.f(rx) * weight_func(ry, EPS, opt);
                end
            end
        end

        ret(k) = ret(k) * h^2;
        increment(progress);
    end
    %% plot the convergence rate
    g = 2.0./(floor(grow_rate.^(1:K) * base_grid));
    
    err = abs(ret - acc);
    gamma = ( g.^(beta) * err) / norm(g.^(beta))^2;
    
    loglog(g, err, '-bo',  g, gamma * g.^(beta), '-r');
    
    legend_handler = legend('quadrature error', sprintf('O(h^{%1.1f})', beta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');
    fontsize(gca, 15,'points');
    grid on;
else
    progress = PoolWaitbar(K*S, 'Starting');
    parfor l = 1:K*S
        % random rotation
        theta = rand() * 2 * pi;
        v = [cos(theta), sin(theta)];

        [k, s] = ind2sub([K,S], l);

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

        % make a random shift of grid.
        X = X + shift(1);
        Y = Y + shift(2);
        
        M = (N + 1)^2;

        Data = zeros(M , 1);
    
        for i = 1:M
            rx = X(i) * v(1) + Y(i) * v(2);
            ry = X(i) * v(2) - Y(i) * v(1);
            
            if rx >= -0.5 && rx <= 0.5
                if abs(ry) <= EPS
                    ret(l) = ret(l) +  opt.f(rx) * weight_func(ry, EPS, opt);
                end
            end
        end

        ret(l) = ret(l) * h^2;
    
        increment(progress);
    end

    ret = reshape(ret, K, S);

    %% plot the convergence rate
    g = 2./( floor(base_grid * grow_rate.^(1:K)));

    var_err = sum( (ret - acc).^2, 2)/S;
    gamma = ( g.^(delta) * var_err) / norm(g.^(delta))^2;

    loglog(g, var_err, '-bo',  g, gamma * g.^(delta), '-r');
    legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', delta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');
    fontsize(gca, 15,'points');
    grid on;
end


F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

filename = sprintf('2D-T%s-Q%d-K%d-S%d-A%6.4f.png',opt.type, opt.q, K, S, alpha);
exportgraphics(gca,filename,'Resolution',300);
close all;

