function plot_traj(traj,x0,xf,R_planet)
    R_scale = 6378e3;
    R_planet = R_planet * R_scale
  %% Plot optimized trajectory with velocity arrows
    figure;
    plot(traj(1,:)*R_scale, traj(2,:)*R_scale, 'k-o','LineWidth',2); hold on;
    plot(x0(1)*R_scale, x0(2)*R_scale,'go','MarkerSize',10,'MarkerFaceColor','g');
    plot(xf(1)*R_scale, xf(2)*R_scale,'ro','MarkerSize',10,'MarkerFaceColor','r');
    
    
    xlabel('X [m]');
    ylabel('Y [m]');
    title('SCP Convexified 2D Trajectory with Gravity');
    legend('Trajectory','Start','Target','Location','best');

    % Velocity arrow at final node
    vx_f = traj(3,end);
    vy_f = traj(4,end);
    quiver(traj(1,end)*R_scale, traj(2,end)*R_scale, vx_f, vy_f, ...
       0.5, 'r', 'LineWidth',2, 'MaxHeadSize',3); 

    planet_center = [0,0];  % assume at origin
    
    % Draw planet
    theta = linspace(0,2*pi,200);
    x_planet = planet_center(1) + R_planet*cos(theta);
    y_planet = planet_center(2) + R_planet*sin(theta);
    fill(x_planet, y_planet, [0.6 0.8 1], 'FaceAlpha',0.5, 'EdgeColor','b'); % semi-transparent


    grid on; axis equal;


end