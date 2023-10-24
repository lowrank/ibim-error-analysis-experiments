%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% 3D Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('./'));

opt.type = 'sphere';
opt.R = 0.75;
opt.f = @(x,y,z)(( cos(x.^2 - y - z.^3)) .* sin(y.^2 - x.^3 - z)); 

% accurate integral
acc = 1.71628896296959  * 0.75^2;

% regularity of weight function
opt.q = 1;
    
% tube width
alpha = 1;
beta = 1 + (opt.q + 1)* (1 - alpha);

K = 16;
S = 32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ret = zeros(K*S, 1);

base_grid = 10;
grow_rate = 1.2;

if S == 1
    opt.random = false;
    progress = PoolWaitbar(K, 'Starting');
    for k = 1:K
        
        N = floor(grow_rate^k * base_grid);
        
        h = 2 / N;
    
        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * h^alpha;
        end
        ret(k) = ibim_quadrature_3d(N, EPS, opt);
        increment(progress);
    end
    
    %% plot the convergence rate
    % plot the convergence rate
    g = 2./(floor(base_grid * grow_rate.^(1:K)));

    err = abs(ret - acc);
    gamma = ( g.^(beta) * err) / norm(g.^(beta))^2;

    loglog(g, err, '-bo',  g, gamma * g.^(beta), '-r');

    legend_handler = legend('quadrature error', sprintf('O(h^{%1.1f})', beta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');
    fontsize(gca, 15,'points');
    grid on;
else 
    opt.random = true;
    progress = PoolWaitbar(K*S, 'Starting');

    for l = 1:K*S
        [k, s] = ind2sub([K, S], l);
        N = floor(grow_rate^k * base_grid);
        h = 2 / N;
    
        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * h^alpha;
        end

        ret(l) = ibim_quadrature_3d(N, EPS, opt);
        increment(progress);
    end
    

    ret = reshape(ret, K, S);
    %% plot the convergence rate
    g = 2./( floor(base_grid * grow_rate.^(1:K)));

    var_err = sum( (ret - acc).^2, 2)/S;
    gamma = ( g.^(2*beta) * var_err) / norm(g.^(2*beta))^2;

    loglog(g, var_err, '-bo',  g, gamma * g.^(2*beta), '-r');
    legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', 2*beta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');
    fontsize(gca, 15,'points');
    grid on;
    
end    


F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F);

filename = sprintf('3D-T%s-Q%d-K%d-S%d-A%6.4f.png',opt.type, opt.q, K, S, alpha);
exportgraphics(gca,filename,'Resolution',300);
close all;
