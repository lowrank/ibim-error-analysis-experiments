function [dist, normal, Jac] = dist2curve(x, y, opt)
% DIST2CURVE  Compute the distance between a point and a curve

    if strcmp(opt.type,  'circle')
        dist = sqrt(x^2 + y^2) - opt.R;
        normal = [x, y] / sqrt(x^2 + y^2);
        Jac  = opt.R / (dist + opt.R);
    elseif strcmp(opt.type, 'ellipse')
        
        distSquared = (x - opt.anchor(:,1)).^2 + (y - opt.anchor(:,2)).^2;
        [~, idx] = min(distSquared);
        minAngle = opt.angles(idx);

        objFunc = @(angle) ( (opt.Rx * cos(angle) - x )^2 + ( opt.Ry * sin(angle) - y )^2 ) ;
        
        [projAngle, ~, ~, ~] =fminunc(objFunc, minAngle, opt.min_options);

        if x^2/opt.Rx^2 + y^2/opt.Ry^2 > 1
            dist = sqrt( (x - opt.Rx * cos(projAngle))^2 + (y - opt.Ry * sin(projAngle))^2 );
        else
            dist = -sqrt( (x - opt.Rx * cos(projAngle))^2 + (y - opt.Ry * sin(projAngle))^2 );
        end 

        normal = [opt.Rx * cos(projAngle) * opt.Ry^2, opt.Ry * sin(projAngle) * opt.Rx^2];
        normal = normal / sqrt(normal(1)^2 + normal(2)^2);

        % Jac =  1 + dist * Laplace (dist), a 4th order scheme should
        % suffice for low regularity weight functions and thin width.

        Jac = 1;
    end
end