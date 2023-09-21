function [dist, normal, Jac] = dist2curve(x, y, opt)
% DIST2CURVE  Compute the distance between a point and a curve

    if strcmp(opt.type,  'circle')
        dist = sqrt(x^2 + y^2) - opt.R;
        normal = [x, y] / sqrt(x^2 + y^2);
        Jac  = opt.R / (dist + opt.R);
    elseif strcmp(opt.type, 'ellipse')
        % it uses a polynomial rootfinding algorithm (should be analytic).
        [dist, normal] = ellipse_subroutine(x, y, opt.Rx, opt.Ry);
        % Jac under construction. For approximation at thin tube, Jac = 1
        % is acceptable, otherwise an error of (tube-width^2) is
        % introduced.
        Jac = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ELLIPSE SUBROUTINE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dist, normal] = ellipse_subroutine(x, y, a, b)
% The distance is only well-defined near the boundary, otherwise the result
% could be completely wrong.

    if x^2/a^2 + y^2/b^2 > 1
        sgn = 1;
    else
        sgn = -1;
    end

    % reference:  Convex functions: constructions, characterizations and
    % counterexamples, J.M. Borwein & J.D. Vanderwerff (2010). 
    
    % The projection must be in a form of 
    % (\frac{a^2 x}{a^2 - t}, \frac{b^2 y}{b^2- t}).
    
    % Then solve quartic polynomial.
    B0 = -1;
    B1 = 2*(a^2 + b^2);
    B2 = (-a^4 - 4*a^2*b^2 - b^4 + a^2 * x^2 + b^2 * y^2);
    B3 = 2* a^4*b^2 + 2*a^2*b^4 - 2*a^2*b^2*(x^2+y^2);
    B4 = -a^4*b^4 + a^2*b^2 * (b^2 * x^2 + a^2 * y^2);

    p1 = 2*B2^3 - 9 * B1 * B2 * B3 + 27 * B0 * B3^2 + 27*B1^2 *B4 - 72*B0*B2*B4;
    D0 = B2^2 - 3 * B1*B3 + 12*B0 * B4;
    p2 = p1 + (p1^2 - 4*D0^3)^(1/2);
    D1 = (p2/2)^(1/3);
    p3 = D0/3/B0/D1 + D1/3/B0;

    p4 = (B1^2/4/B0^2 - 2*B2/3/B0 + p3)^(1/2);
    p5 = B1^2/2/B0^2 - 4*B2/3/B0 - p3;
    p6 = (-B1^3/B0^3 + 4*B1*B2/B0^2 - 8 *B3/B0)/4/ p4;

    x1 = -B1/4/B0 - p4/2 + ((p5 - p6))^(0.5)/2;
    x2 = -B1/4/B0 - p4/2 - ((p5 - p6))^(0.5)/2;
    x3 = -B1/4/B0 + p4/2 + ((p5 + p6))^(0.5)/2;
    x4 = -B1/4/B0 + p4/2 - ((p5 + p6))^(0.5)/2;

    t_candidates = [x1, x2, x3, x4];
    t_candidates = real(t_candidates(abs(t_candidates - real(t_candidates)) < 1e-9));
    
    dist = inf;
    
    for candiate_id = 1:length(t_candidates)
        Px = a^2 * x/(a^2 - t_candidates(candiate_id));
        Py = b^2 * y/(b^2 - t_candidates(candiate_id));

        cur_dist =  sqrt((Px - x)^2 + (Py - y)^2);

        if dist > cur_dist
            dist = cur_dist;
            normal = [Px/a^2, Py/b^2];
            normal = normal/norm(normal);
        end
    end

    dist = dist * sgn;
end