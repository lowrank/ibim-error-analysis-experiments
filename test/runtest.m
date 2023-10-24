
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tests for QTREE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


qt = qtree(2, 0.05, 0.1, true, true);

opt.type = 'star';
opt.R = 0.75;
opt.r = 0.2;
opt.m = 3;
angles = linspace(0, 2*pi, 100);
opt.angles = angles(1:end-1);
opt.anchor = [ (opt.R + opt.r * cos(opt.m * opt.angles) ) .* cos(opt.angles) ;...
    (opt.R + opt.r * cos(opt.m * opt.angles) ) .* sin(opt.angles)]';
opt.min_options = optimoptions('fminunc','Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);


qt.distFunc = @(x) (dist2curve(x(1), x(2), opt));

qt.populate();
disp(qt)

pts = zeros(2, length(qt.validList));
for i = 1:length(qt.validList)
    id = qt.validList(i);
    pts(:,i) = qt.dict{id}.center;
end

pts = pts(:,all(pts,1));
scatter(pts(1, :), pts(2, :));
daspect([1 1 1]);
grid on;
grid minor