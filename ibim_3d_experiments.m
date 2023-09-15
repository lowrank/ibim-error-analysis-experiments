%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% 3D Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opt.type = 'sphere';
opt.R = 0.75;
opt.f = @(x,y,z)(( cos(x.^2 - y - z.^3)) .* sin(y.^2 - x.^3 - z)); 

%% accurate integral
acc = 1.71628896296959  * 0.75^2;

%% regularity of weight function
opt.q = 2;
    
%% tube width
alpha = 1;
beta = 1 + (opt.q + 1)* (1 - alpha);

K = 20;
S = 32;
ret = zeros(K, S);

base_grid = 10;
grow_rate = 1.2;

if S == 1
    opt.random = false;
    for k = 1:K
        
        N = floor(grow_rate^k * base_grid);
        
        h = 2 / N;
    
        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * h^alpha;
        end
        ret(k) = ibim_quadrature_3d(N, EPS, opt);

    end
    
    %% plot the convergence rate
    g = 2./( floor(base_grid * grow_rate.^(1:K)));
    
    loglog(g,...
           abs(ret -acc), ...
           '-bo', ...
           g, ...
           0.03*g.^beta);
    
    legend_handler = legend('quadrature error', sprintf('O(h^{%1.1f})', beta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');

    grid on;
else 
    opt.random = true;
    for k = 1:K

        N = floor(grow_rate^k * base_grid);
        
        h = 2 / N;
    
        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * h^alpha;
        end

        parfor s = 1:S
            % add random translation
            
            ret(k, s) = ibim_quadrature_3d(N, EPS, opt);
        end

    end
    
    %% plot the convergence rate
    g = 2./( floor(base_grid * grow_rate.^(1:K)));
    
    loglog(g,...
           sum( (ret - acc).^2, 2)/S, ...
           '-bo', ...
           g, ...
           0.0003*g.^(2*beta));
    
    legend_handler = legend('error variance', sprintf('O(h^{%1.1f})', 2*beta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');

    grid on;
    
end    