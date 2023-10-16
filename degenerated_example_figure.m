

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical Experiments for Validation of Numerical Analysis of IBIM
%
% Generate Demonstration Figure
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

plot(g_x,g_x,'r','Linewidth',2);
hold on;

plot(g_x,g_x + 0.6,'r--','Linewidth',2);
hold on;

plot(g_x,g_x - 0.6,'r--','Linewidth',2);
hold on;

for i = 1:length(X)^2
    if abs(X(i) - Y(i)) < 0.6
    plot([X(i) (X(i) + Y(i))/2], [Y(i) (X(i) + Y(i))/2], '--b');
    hold on;
    end
end

daspect([1 1 1])