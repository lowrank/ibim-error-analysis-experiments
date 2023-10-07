function ret = ibim_quadrature_2d(N, EPS, opt)

    h   = 2 / N;
    qt = qtree(2, 2/N, EPS, opt.random); % no rotation.
    qt.distFunc = @(x) (dist2curve(x(1), x(2), opt));
    qt.populate();
    M = length(qt.validList);

    ret = 0;  % initialize the integral.
    for i = 1:M
        id = qt.validList(i);

        X = qt.dict{id}.center(1);
        Y = qt.dict{id}.center(2);
        
        [d, n, J] = dist2curve(X, Y, opt);
        % projection.
        px = X - d * n(1);
        py = Y - d * n(2);

        % compute the integral.
        if abs(d) < EPS 
            ret = ret + opt.f(px, py) * weight_func(d, EPS, opt) * J; 
        end
        
    end

    ret = ret * h^2;

end