function y = weight_func(s, EPS, opt)
    if opt.q == 1
        y = (1 - abs(s)/EPS) / EPS;
    elseif opt.q == 2
        y = (1 + cos(pi*s/EPS)) / (2*EPS);
    end
end