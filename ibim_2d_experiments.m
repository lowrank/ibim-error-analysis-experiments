
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% 2D Experiments
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('./'));

opt.type = 'circle'; % type of the domain
opt.R  = 0.75;       % radius of the circle
opt.f =  @(x, y)( cos(x.^2 - y)) .* sin(y.^2 - x.^3); % test integrand 

% call Mathematica to produce accurate value.
if check_installation()
    F = sprintf(['''NIntegrate[%f * Cos[%f^2 * Cos[t]^2 - %f * Sin[t]]', ...
        '* Sin[%f^2 * Sin[t]^2 - %f^3 * Cos[t]^3],', ...
        '{t, 0, 2*Pi}, WorkingPrecision->15]'''],...
        opt.R, opt.R, opt.R, opt.R, opt.R);  % Mathematica Eval
    
    [status, cmdout] = system(sprintf(['wolframscript -code ', F]));
    acc = str2double( cmdout(1:16) );
else 
    acc = 0.992902336403640;
end

fprintf('accurate integral is %1.15f by Mathematica.\n', acc);

%% regularity
opt.q = 2; % order of the regularity 

%% tube width
alpha = 0.5;  % parameter for the tube width
beta  = 0.5 + (opt.q + 1) * (1 - alpha);  % theoretical value of decay rate.

K = 24; % number of grid sizes.
S = 32;  % number of sampled translations.

ret = zeros(K, S); % return value of the integral.

base_grid = 100;
grow_rate = 1.2;

if S == 1 % no random translation is applied
    opt.random = false;
    progress = PoolWaitbar(K, 'Starting');
    parfor k = 1:K
        N = floor(grow_rate^k * base_grid);
        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * (2/N)^alpha;
        end
        ret(k) = ibim_quadrature_2d(N, EPS, opt);
        increment(progress);
    end
   

    % plot the convergence rate
    g = 2./(floor(base_grid * grow_rate.^(1:K)));

    err = abs(ret - acc);
    gamma = ( g.^(beta) * err) / norm(g.^(beta))^2;

    loglog(g, err, '-bo',  g, gamma * g.^(beta), '-r');

    legend_handler = legend('quadrature error', sprintf('O(h^{%1.1f})', beta), 'Location', 'northwest');
    fontsize(legend_handler,18,'points');

    grid on;


else % random translations are used
    opt.random = true;
    progress = PoolWaitbar(K*S, 'Starting');
    parfor l = 1:K*S
        
        [k, s] = ind2sub([K, S], l);

        N = floor(grow_rate^k * base_grid);

        if alpha == 0
            EPS = 0.1;
        else
            EPS = 2 * (2/N)^alpha;
        end


        ret(l) = ibim_quadrature_2d(N, EPS, opt);
        
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

    grid on;
    

end




