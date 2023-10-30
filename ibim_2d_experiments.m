
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% 2D Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ret=ibim_2d_experiments(K, S, opt, base)

    addpath(genpath('./'));
    
    ret = zeros(K, S); % return value of the integral.
    
    if nargin == 3
        base_grid = 10;
    else
        base_grid = base;
    end
    grow_rate = 1.2;
    
    if S == 1 % no random translation is applied
        opt.random = false;
        progress = PoolWaitbar(K, 'Starting');
        parfor k = 1:K
            N = floor(grow_rate^k * base_grid);
            if opt.alpha == 0
                EPS = 0.1;
            else
                EPS = 2 * (2/N)^opt.alpha;
            end
            ret(k) = ibim_quadrature_2d(N, EPS, opt);
            increment(progress);
        end
       
    
        %% plot the convergence rate
        g = 2./(floor(base_grid * grow_rate.^(1:K)));
    
        err = abs(ret - opt.acc);
        gamma = ( g.^(opt.beta) * err) / norm(g.^(opt.beta))^2;
    
        loglog(g, err, '-bo',  g, gamma * g.^(opt.beta), '-r');
    
        legend_handler = legend('quadrature error', sprintf('O(h^{%1.1f})', opt.beta), 'Location', 'northwest');
        fontsize(legend_handler,18,'points');
        fontsize(gca, 15,'points');
        grid on;
    
    
    else % random translations are used
        opt.random = true;
        progress = PoolWaitbar(K*S, 'Starting');
        parfor l = 1:K*S
            
            [k, ~] = ind2sub([K, S], l);
    
            N = floor(grow_rate^k * base_grid);
    
            if opt.alpha == 0
                EPS = 0.1;
            else
                EPS = 2 * (2/N)^opt.alpha;
            end
    
    
            ret(l) = ibim_quadrature_2d(N, EPS, opt);
            
            increment(progress);
        end
    
        ret = reshape(ret, K, S);
    
        %% plot the convergence rate
        g = 2./( floor(base_grid * grow_rate.^(1:K)));
    
        var_err = sum( (ret - opt.acc).^2, 2)/S;
        gamma = ( g.^(2*opt.beta) * var_err) / norm(g.^(2*opt.beta))^2;
    
        if isfield(opt, 'upsilon')
            tau = ( g.^(2*opt.upsilon) * var_err) / norm(g.^(2*opt.upsilon))^2;
            loglog(g, var_err, '-bo',  g, gamma * g.^(2*opt.beta), '-r', g, tau * g.^(2*opt.upsilon), '--k');
            legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', 2*opt.beta), sprintf('O(h^{%1.1f})', 2*opt.upsilon), 'Location', 'northwest');
            fontsize(legend_handler,18,'points');
            fontsize(gca, 15,'points');
            grid on;
        else
            loglog(g, var_err, '-bo',  g, gamma * g.^(2*opt.beta), '-r');
            legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', 2*opt.beta), 'Location', 'northwest');
            fontsize(legend_handler,18,'points');
            fontsize(gca, 15,'points');
            grid on;
        end
    
    end
    
    F = findall(0,'type','figure','tag','TMWWaitbar');
    delete(F);
    
    filename = sprintf('2D-T%s-Q%d-K%d-S%d-A%6.4f.png',opt.type, opt.q, K, S, opt.alpha);
    exportgraphics(gca,filename,'Resolution',300);
    close all;

end