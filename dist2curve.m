function [dist, normal, Jac] = dist2curve(x, y, opt)
% DIST2CURVE  Compute the distance between a point and a curve

    if strcmp(opt.type,  'circle')
        dist = sqrt(x^2 + y^2) - opt.R;
        normal = [x, y] / sqrt(x^2 + y^2);
        Jac  = opt.R / (dist + opt.R);
    end
end