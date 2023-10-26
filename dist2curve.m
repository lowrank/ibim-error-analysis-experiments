function [dist, normal, Jac] = dist2curve(x, y, opt)
% DIST2CURVE  Compute the distance between a point and a curve

    if strcmp(opt.type,  'circle') || strcmp(opt.type,  'semicircle')
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

        kappa = (opt.Rx * opt.Ry)/(opt.Rx^2 * sin(projAngle)^2 + opt.Ry^2 * cos(projAngle)^2)^(3/2);
        
        Jac = 1/(1 + kappa * dist);

    elseif strcmp(opt.type, 'degenerate')
        distSquared = (x - opt.anchor(:,1)).^2 + (y - opt.anchor(:,2)).^2;
        [~, idx] = min(distSquared);
        minAngle = opt.angles(idx);

        objFunc = @(angle) ( (sqrt( opt.Rx * abs(cos(angle)) ) .* sign(cos(angle))  - x )^2 + ( opt.Ry * sin(angle) - y )^2 ) ;
        
        [projAngle, ~, ~, ~] =fminunc(objFunc, minAngle, opt.min_options);

        if x^4/opt.Rx^2 + y^2/opt.Ry^2 >= 1
            dist = sqrt( (x - sqrt(opt.Rx * abs(cos(projAngle))) * sign(cos(projAngle)))^2 + (y - opt.Ry * sin(projAngle))^2 );
        else
            dist = -sqrt( (x - sqrt(opt.Rx * abs(cos(projAngle))) * sign(cos(projAngle)))^2 + (y - opt.Ry * sin(projAngle))^2 );
        end 

        normal = [2*(opt.Rx * cos(projAngle))*sqrt(opt.Rx * abs(cos(projAngle))) / opt.Rx^2, opt.Ry * sin(projAngle) / opt.Ry^2];
        normal = normal / sqrt(normal(1)^2 + normal(2)^2);

        % Jac =  1 + dist * Laplace (dist), a 4th order scheme should
        % suffice for low regularity weight functions and thin width.
        
        kappa = sqrt(opt.Rx) * opt.Ry * (1/2 + sin(projAngle)^2/4) * abs(cos(projAngle))/(opt.Rx*sin(projAngle)^2/4 + opt.Ry^2 * abs(cos(projAngle))^3)^(3/2);

        Jac = 1/(1 + kappa * dist);

    elseif strcmp(opt.type, 'star')
        distSquared = (x - opt.anchor(:,1)).^2 + (y - opt.anchor(:,2)).^2;
        [~, idx] = min(distSquared);
        minAngle = opt.angles(idx);
        objFunc = @(angle) ( ( (opt.R + opt.r * cos(opt.m * angle) ) * cos(angle) - x )^2 + ...
            ( (opt.R + opt.r * cos(opt.m * angle) ) * sin(angle) - y )^2 ) ;

        [projAngle, ~, ~, ~] =fminunc(objFunc, minAngle, opt.min_options);

        L = sqrt(x^2 + y^2);
        rho = opt.R + opt.r * cos(opt.m * projAngle);
        rhop = -opt.r * opt.m * sin(opt.m * projAngle);
        if L > rho
            dist = sqrt((x - (opt.R + opt.r * cos(opt.m * projAngle))*cos(projAngle))^2 +...
                (y - (opt.R + opt.r * cos(opt.m * projAngle))*sin(projAngle))^2   );
        else 
            dist = -sqrt((x - (opt.R + opt.r * cos(opt.m * projAngle))*cos(projAngle))^2 +...
                (y - (opt.R + opt.r * cos(opt.m * projAngle))*sin(projAngle))^2   );
        end
        normal = [rho * cos(projAngle) + rhop * sin(projAngle), rho * sin(projAngle) - rhop * cos(projAngle) ];
        normal = normal / sqrt(normal(1)^2 + normal(2)^2);


        rhopp = -opt.r * opt.m^2 * cos(opt.m * projAngle);
        kappa = (2 * rhop^2 - rho * rhopp + rho^2)/(rhop^2 + rho^2)^(3/2);

        Jac = 1/(1 + kappa * dist);

    elseif strcmp(opt.type, "polynomial")
        distSquared = (x - opt.anchor(:,1)).^2 + (y - opt.anchor(:,2)).^2;
        [~, idx] = min(distSquared);
        min_x = opt.pts(idx);

        objFunc = @(pts) ( ( pts - x )^2 + ( (pts)^6 - y )^2 ) ;
        [pts_m, ~, ~, ~] = fminunc(objFunc, min_x, opt.min_options);

        if y > (pts_m)^6
            dist = - sqrt((x - pts_m)^2 + (y - (pts_m)^6)^2 );
        else
            dist =   sqrt((x - pts_m)^2 + (y - (pts_m)^6)^2 );
        end
        normal = [1, 6 * (pts_m)^5];
        normal = normal / sqrt(normal(1)^2 + normal(2)^2);

        kappa = 30*(pts_m)^4/(1 + 36*(pts_m)^(10))^(3/2);
        Jac = 1/(1 + kappa * dist);
    end


end