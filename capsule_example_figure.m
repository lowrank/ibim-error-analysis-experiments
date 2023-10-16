

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% Generate Demonstration Figure/Capsule
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ylim([-1 1]);
xlim([-1 1]);
% gridlines ---------------------------
hold on
g_y=-1:0.4:1; % user defined grid Y [start:spaces:end]
g_x=-1:0.4:1; % user defined grid X [start:spaces:end]

for i=1:length(g_x)
   plot([g_x(i) g_x(i)],[g_y(1) g_y(end)],'k:') %y grid lines
   hold on    
end

for i=1:length(g_y)
   plot([g_x(1) g_x(end)],[g_y(i) g_y(i)],'k:') %x grid lines
   hold on    
end

[X,Y] = meshgrid(g_x, g_y);

scatter(X, Y, 40,'MarkerEdgeColor',[0 .5 .5],...
         'MarkerFaceColor',[0 .7 .7],...
         'LineWidth',1.5);
hold on;

theta = atan(sqrt(2));

pts = [0.5, 0.2; 0.5, -0.2; -0.5, -0.2; -0.5, 0.2]; % 4 x 2
rot = [cos(theta), sin(theta); -sin(theta), cos(theta)];

pts = pts * rot;

plot(pts([2,3], 1), pts([2,3], 2),'r','Linewidth',2);
hold on;

plot(pts([1,4], 1), pts([1,4], 2),'r','Linewidth',2);
hold on;

M = 100;
angles = linspace(theta - pi/2, theta + pi/2, M);
center = [0.5, 0] * rot;
pts = 0.2*[cos(angles);sin(angles)]';
for i = 1:M
   plot(center(1) + pts(:, 1), center(2) + pts(:, 2), 'r', 'Linewidth', 2)
end

angles = linspace(theta + pi/2, theta + 3*pi/2, M);
center = [-0.5, 0] * rot;
pts = 0.2*[cos(angles);sin(angles)]';
for i = 1:M
   plot(center(1) + pts(:, 1), center(2) + pts(:, 2), 'r', 'Linewidth', 2)
end

daspect([1 1 1])