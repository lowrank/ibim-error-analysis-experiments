function ret = ibim_quadrature_3d(N, EPS, opt)

    if opt.random
        shift = EPS/10 * [rand(), rand(), rand()];
    else
        shift = [0, 0, 0];
    end

    coordinates = linspace(-1, 1, N+1);
    h = 2/ N;
    [X, Y, Z] = meshgrid(coordinates);

    % make a shift of grid.
    X = X + shift(1);
    Y = Y + shift(2);
    Z = Z + shift(3);

    M = (N+1)^3;
    vol = h^3;
    ret = 0;
    for i = 1:M
        [d, n, J] = dist2surf(X(i), Y(i), Z(i), opt);
        px = X(i) - n(1) * d;
        py = Y(i) - n(2) * d;
        pz = Z(i) - n(3) * d;

        
        if abs(d) <= EPS
            ret = ret + opt.f(px, py, pz) * weight_func(d, EPS, opt) * J;
        end
    end

    ret = ret * vol;
end