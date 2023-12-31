%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% 3D Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ibim_3d_experiments(K, S, opt)
    addpath(genpath('./'));
    
    ret = zeros(K*S, 1);
    
    if nargin == 3
        base_grid = 10;
    else
        base_grid = base;
    end
    grow_rate = 1.2;
    
    if S == 1
        opt.random = false;
        progress = PoolWaitbar(K, 'Starting');
        parfor k = 1:K
            
            N = floor(grow_rate^k * base_grid);
            
            h = 2 / N;
        
            if opt.alpha == 0
                EPS = 0.1;
            else
                EPS = 2 * h^opt.alpha;
            end
            ret(k) = ibim_quadrature_3d(N, EPS, opt);
            increment(progress);
        end
        
        %% plot the convergence rate
        % plot the convergence rate
        g = 2./(floor(base_grid * grow_rate.^(1:K)));
    
        err = abs(ret - opt.acc);
        gamma = ( g.^(opt.beta) * err) / norm(g.^(opt.beta))^2;
    
        loglog(g, err, '-bo',  g, gamma * g.^(opt.beta), '-r');
    
        legend_handler = legend('quadrature error', sprintf('O(h^{%1.1f})', opt.beta), 'Location', 'northwest');
        fontsize(legend_handler,18,'points');
        fontsize(gca, 15,'points');
        grid on;
    else 
        opt.random = true;
        progress = PoolWaitbar(K*S, 'Starting');
    
        parfor l = 1:K*S
            [k, ~] = ind2sub([K, S], l);
            N = floor(grow_rate^k * base_grid);
            h = 2 / N;
        
            if opt.alpha == 0
                EPS = 0.1;
            else
                EPS = 2 * h^opt.alpha;
            end
    
            ret(l) = ibim_quadrature_3d(N, EPS, opt);
            increment(progress);
        end
        
    
        ret = reshape(ret, K, S);
        %% plot the convergence rate
        g = 2./( floor(base_grid * grow_rate.^(1:K)));
    
        var_err = sum( (ret - opt.acc).^2, 2)/S;
        gamma = ( g.^(2*opt.beta) * var_err) / norm(g.^(2*opt.beta))^2;
    
        loglog(g, var_err, '-bo',  g, gamma * g.^(2*opt.beta), '-r');
        legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', 2*opt.beta), 'Location', 'northwest');
        fontsize(legend_handler,18,'points');
        fontsize(gca, 15,'points');
        grid on;
        
    end    
    
    
    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F);
    
    filename = sprintf('3D-T%s-Q%d-K%d-S%d-A%6.4f.png',opt.type, opt.q, K, S, opt.alpha);
    exportgraphics(gca,filename,'Resolution',300);
    close all;

end
