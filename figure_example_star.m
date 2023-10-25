

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

M = 100;
angles = linspace(theta, theta + 2*pi, M);
center = [0, 0];
rho = 0.75 + 0.2*cos(3 * angles);
pts = [rho.* cos(angles); rho.*sin(angles)]';
for i = 1:M
   plot(center(1) + pts(:, 1), center(2) + pts(:, 2), 'r', 'Linewidth', 2)
end

daspect([1 1 1]);
fontsize(gca, 15,'points');
exportgraphics(gca,'star.png','Resolution',300);
close all;
