function [dist, normal, Jac] = dist2surf(x, y, z, opt)
    if strcmp(opt.type, 'sphere')
        dist = sqrt(x.^2 + y.^2 + z.^2) - opt.R;
        normal = [x, y, z];
        normal = normal / norm(normal);
        Jac = ( opt.R / (opt.R + dist) )^2;
    end
end