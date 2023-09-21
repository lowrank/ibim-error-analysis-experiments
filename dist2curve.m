function [dist, normal, Jac] = dist2curve(x, y, opt)
% DIST2CURVE  Compute the distance between a point and a curve

    if strcmp(opt.type,  'circle')
        dist = sqrt(x^2 + y^2) - opt.R;
        normal = [x, y] / sqrt(x^2 + y^2);
        Jac  = opt.R / (dist + opt.R);
    elseif strcmp(opt.type, 'ellipse')
        % it uses a polynomial rootfinding algorithm (should be analytic).
        dist = ellipse_subroutine(x, y, opt.Rx, opt.Ry);
        % normal and Jac are under construction.
        normal = [];
        Jac = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ELLIPSE SUBROUTINE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = ellipse_subroutine(x, y, a, b)
    L = a^2 - b^2;
    G0 = x^2/a^2 + y^2/b^2 - 1;
    T0 = x^2 + y^2 - a^2 - b^2;
    S0 = x^2/a^4 + y^2/b^4;

    B0 = L^2;
    B1 = -2*L*(L*(a^2 + b^2 + x^2 + y^2) + a^2 * y^2 - b^2 * x^2);
    B2 = 6*L*(a^4*y^2+a^2*y^4 - b^4*x^2 - b^2*x^4 + L*(a^2*b^2 + x^2*y^2)) +...
        (L^2 - (a^2*x^2 + b^2*y^2))^2;
    B3 = -2*a^2*b^2 *(a^2*b^2 * T0 * G0^2 -...
        ((a^2+b^2)*T0^2 + 3*a^2*b^2*T0 - 6*a^4*b^4*S0)*G0 + 2*a^2*b^2*T0^2*S0);
    B4 = a^4*b^4*G0^2*(T0^2 + 4*a^2*b^2*G0);
   
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

    dist_candidates = sqrt( [x1, x2, x3, x4] );

    % find the smallest real number.
    dist = min(dist_candidates(dist_candidates == real(dist_candidates)));

end