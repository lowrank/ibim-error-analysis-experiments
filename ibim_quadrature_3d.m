function ret = ibim_quadrature_3d(N, EPS, opt)

    h   = 2 / N;
    vol = h^3;
    qt = qtree(3, 2/N, EPS, opt.random); % no rotation.
    qt.distFunc = @(x) (dist2surf(x(1), x(2), x(3), opt));
    qt.populate();
    M = length(qt.validList);

    ret = 0;
    for i = 1:M
        id = qt.validList(i);
        
        X = qt.dict{id}.center(1);
        Y = qt.dict{id}.center(2);
        Z = qt.dict{id}.center(3);
        [d, n, J] = dist2surf(X, Y, Z, opt);
        px = X - n(1) * d;
        py = Y - n(2) * d;
        pz = Z - n(3) * d;

        
        if abs(d) <= EPS
            ret = ret + opt.f(px, py, pz) * weight_func(d, EPS, opt) * J;
        end
        
    end

    ret = ret * vol;
end