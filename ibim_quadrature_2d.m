function ret = ibim_quadrature_2d(N, EPS, opt)
     
    if opt.random
        shift = EPS/10 * [rand(), rand()];
    else
        shift = [0, 0];
    end

    pts = linspace (-1, 1, N + 1);
    h   = 2 / N;

    [X, Y] = meshgrid (pts);
    M      = (N+1)^2;

    % make a shift of grid.
    X = X + shift(1);
    Y = Y + shift(2);

    ret = 0;  % initialize the integral.
    for i = 1:M

        [d, n, J] = dist2curve(X(i), Y(i), opt);
        % projection.
        px = X(i) - d * n(1);
        py = Y(i) - d * n(2);

        % compute the integral.
        if abs(d) < EPS 
            ret = ret + opt.f(px, py) * weight_func(d, EPS, opt) * J; 
        end

    end

    ret = ret * h^2;

end