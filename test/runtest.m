
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tests for QTREE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


qt = qtree(2, 0.1, 0.2, false, true);

qt.distFunc = @(x) (norm(x) - 0.75);
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