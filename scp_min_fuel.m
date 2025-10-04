
%% Sequential Convex Programming (SCP) for 2D Minimum-Fuel STO

function scp_min_fuel(x0,xf,planet_radius,mu_real,periods) 
    
    N = 100;    % discretization points
    mu=1; %variables are scaled so that mu and R == 1
    
    T0 = sqrt(planet_radius^3/mu_real);
    Tdim = periods * 2 * pi * sqrt(sqrt(x0(1)^2 + x0(2)^2)^3 / mu_real);
    x0(1:2) = x0(1:2) / planet_radius; 
    xf(1:2) = xf(1:2) / planet_radius;

    T = Tdim / T0;
    a0 = planet_radius / T0^2;
    v0 = planet_radius / T0;
    x0(3:4) = x0(3:4) / v0; 
    xf(3:4) = xf(3:4) / v0;
    % Maximum thrust
    u_max = (2/a0);

    %% SCP Parameters
    max_iter = 15;
    tolerance = 0.01;
    
    
    
    dt = T/N;
   
    r = norm(x0(1:2));   % distance from planet center
    rf = norm(xf(1:2));   % distance from planet center
   % v_orb = abs(sqrt(mu / rf));  % scalar
   % v_current = norm(x0(3:4));
  %  v_dir = x0(3:4) / v_current;  % unit vector along current velocity
  %  delta_v = v_orb - v_current;
    t_hat = [x0(2); x0(1)] / r;
    v_orb = sqrt(mu/r) * t_hat;  % desired orbital velocity
    dv = v_orb - x0(3:4);              % velocity delta to align with orbit
    %U_mag = delta_v;              % approximate constant accel over transfer
    %U_mag = max(min(U_mag, u_max), -u_max);  % clamp to bounds
    %U_mat = repmat(U_mag * v_dir, 1, N);
    %tau = 1;                        % decay factor (tune: 0.1-0.3 works well)
    %weights = exp(-tau * (0:N-1)); 
    U_mat = zeros(2,N) * 1e-6;% U_mat .* weights;
    %U_mat(1:2) = dv/dt;
    % Straight line from x0 to xf
    x_ref = zeros(4,N);
    u = x0(1:2)';
    v = xf(1:2)';
    %theta = acos(dot(,) / (norm(x0(1:2)) * norm(xf(1:2)))); 
    theta = atan2(det([v; u]), dot(v,u));
    if theta < (periods-1) * 360 * pi/180
        theta = theta + 2*pi;
    end
    theta0 = atan(x0(1)/x0(2));
    if(x0(2)<0)
        theta0 = theta0 + pi;
    end

    for k = 1:N
        alpha = (k-1)/(N-1);
        alpha1=0;
        y = cos((theta ) * alpha+ theta0) * ((r*(1-alpha1)) + (rf*alpha1));
        x = sin((theta) * alpha + theta0) * ((r*(1-alpha1)) + (rf*alpha1));
        x_ref(1:2,k) = [x;y];
    end

    for k = 1:N-1
        x_ref(3:4,k) = (x_ref(1:2,k+1) - x_ref(1:2,k)) / dt;
    end
    % for k = 1:N
    %     r = norm(x_ref(1:2,k));
    %     v_orb = sqrt(mu / r);
    %     dir = x_ref(3:4,k) / norm(x_ref(3:4,k));
    %     x_ref(3:4,k) = dir * v_orb;  % set magnitude to local circular speed
    % end


    % For the last node, just copy the previous velocity
    x_ref(3:4,N) = x_ref(3:4,N-1);

 R_planet = 1; 
%    figure; hold on; axis equal; grid on;
% 
% % Plot the planet
%     % scaled radius
%  theta = linspace(0,2*pi,100);
%  fill(R_planet*cos(theta), R_planet*sin(theta), [0.5 0.5 1], 'FaceAlpha',0.3, 'EdgeColor','b');
% % 
% % % Plot the initial trajectory guess
%  plot(x_ref(1,:), x_ref(2,:), 'k-o','LineWidth',2, 'MarkerSize',4);
% % 
% % % Plot start and target
%  plot(x0(1), x0(2), 'go', 'MarkerSize',10,'MarkerFaceColor','g');
%  plot(xf(1), xf(2), 'ro', 'MarkerSize',10,'MarkerFaceColor','r');
% % 
% % % Optionally, draw velocity arrows at every few nodes
%  skip = 5;
%  quiver(x_ref(1,1:skip:end), x_ref(2,1:skip:end), ...
%         x_ref(3,1:skip:end), x_ref(4,1:skip:end), 0.5, 'r');
% 
%  xlabel('X [scaled]');
%  ylabel('Y [scaled]');
%  title('Initial Curved Trajectory Guess');
%  legend('Planet','Trajectory','Start','Target','Velocity');

    
    % clip magnitude per node to u_max (scale columns that exceed)
   % mags = sqrt(sum(U_mat.^2,1));
   % scale = min(1, u_max ./ (mags + 1e-12));
   % U_mat = U_mat .* scale;
    U_ref = reshape(U_mat, 2*N, 1);

    
    
    %% SCP Iteration
    for iter = 1:max_iter
        
        % 1. Linearize dynamics around x_ref
        % inside the SCP iteration, after x_ref exists
        [A,B,d] = linearize_dynamics(x_ref,N,dt,mu);
    
           %% --- Objective: fuel + squared-smoothness + total-variation (TV) on vector diffs
        epsilon = 1e-6;
        lambda_TV = 0.05;        % TV weight (promotes piecewise-constant thrusts => longer burns)
        
        fuel = 1;
        objective = @(U) ...        
            fuel* sum(sqrt(U(1:2:end).^2 + U(2:2:end).^2 + epsilon^2)) ...                   % fuel (approx l2)            
            + lambda_TV * sum(sqrt((U(1:2:end-2)-U(3:2:end)).^2 + (U(2:2:end-2)-U(4:2:end)).^2 + epsilon^2)); %...
    
        lb = -u_max*ones(2*N,1);
        ub = u_max*ones(2*N,1);
        options = optimoptions('fmincon', ...
            'Display','iter', ...
            'Algorithm','sqp', ...
            'MaxFunctionEvaluations',1e5, ...
            'StepTolerance', 1e-7, ...
            'ConstraintTolerance', 1e-4, ...
            'OptimalityTolerance',1.5); 
            % <- loosen this);
    
        [U_opt, fval, exitflag, output, lambda] = fmincon(objective, U_ref, [], [], [], [], lb, ub, ...
            @(U) dynamics_constraints_SCP(U, x0, xf, N, A, B,d,dt), options);
    
      %  [~, ceq_test] = dynamics_constraints_SCP(U_opt, x0, xf, N, A, B, d,R_planet);
      %  max_ceq = max(abs(ceq_test));

        % 3. Reconstruct trajectory
           
        U_opt(abs(U_opt) < 0.001) = 0;
        x = x0;
        traj = zeros(4,N);
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            
            traj(:,k) = x;
        end
    
        % 4. Check convergence
        final_error = max(abs(traj(:,end) - xf));
        fprintf('Iteration %d, final state error = %.6f\n', iter, final_error);

        
        %must recalculate dynamics because our path has changed
        [A,B,d] = linearize_dynamics(traj,N,dt,mu); 
        
        x = x0;
        traj = zeros(4,N);
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            traj(:,k) = x;
        end
        
        final_error_recalc = max(abs(traj(:,end) - xf));

        fprintf('Iteration %d fter recalculating gravity, final state error = %.6f\n', iter, final_error_recalc);
        if abs(final_error_recalc - final_error) < 1e-4
            fprintf('Converged!\n');
            break;
        end
    
        % 5. Update reference trajectory & control
        x_ref = traj;
        U_ref = U_opt;
        
    end
    
    plot_control(T,N,U_opt);
    plot_traj(traj,x0,xf,R_planet);
    
    %% Dynamics constraint function
    function [c, ceq] = dynamics_constraints_SCP(U, x0, xf, N, A, B,d,dt)
        R_planet = 1;
    x = x0;
    traj = zeros(2,N);
    for k = 1:N
        uk = U(2*k-1:2*k);
        x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
        traj(:,k) = x(1:2);
    end

%        figure; hold on; axis equal; grid on;
% R_planet = 1;      % scaled radius
%  theta = linspace(0,2*pi,100);
%  fill(R_planet*cos(theta), R_planet*sin(theta), [0.5 0.5 1], 'FaceAlpha',0.3, 'EdgeColor','b');
% plot(traj(1,:), traj(2,:), 'k-o','LineWidth',2, 'MarkerSize',4);


    % Planet constraint (unchanged)
    pos = traj;
    dist2 = sum(pos.^2,1);
    c = (R_planet^2 - dist2)';

    % --- Line intersection constraint ---
    p1 = traj(:,N-1);
    p2 = traj(:,N);
    pt = xf(1:2);

    v = p2 - p1;            % line direction
    
    v_hat = v / norm(v);
    r = pt - p2;            % vector to target
    
    
    eps_perp = 0.001;
    % 1) Perpendicular distance (cross product test)
    cross_val = abs(v(1)*r(2) - v(2)*r(1));   % scalar cross product magnitude
    c_perp = cross_val - eps_perp * norm(v)^2;

    % 2) Along-line tolerance (optional)
    c_along = abs(dot(v_hat, r)) - sqrt(v(1)^2 + v(2)^2) ;

    % Combine
    c = [c; c_perp; c_along];
    ceq = [];
    end
end