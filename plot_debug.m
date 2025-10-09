function plot_debug(traj,x0,xf,R_planet)
    R_scale = 1;
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
    

    planet_center = [0,0];  % assume at origin
    
    % Draw planet
    theta = linspace(0,2*pi,200);
    x_planet = planet_center(1) + R_planet*cos(theta);
    y_planet = planet_center(2) + R_planet*sin(theta);
    fill(x_planet, y_planet, [0.6 0.8 1], 'FaceAlpha',0.5, 'EdgeColor','b'); % semi-transparent


    grid on; axis equal;


end