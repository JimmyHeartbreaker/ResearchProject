
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
    u_max = (5/a0);

    %% SCP Parameters
    max_iter = 15;
    tolerance = 0.03;
    
    x_ref = zeros(4,N); 
    for k = 1:N 
        alpha = k/N; 
        x_ref(:,k) = x0*(1-alpha) + xf*alpha; 
    end
    
    
    R_planet = 1;
    dt = T/N;
       
    
    %% Initial control guess: small random values
    U_ref = zeros(2*N,1);
    U_ref(1:2) = max(-u_max,min(u_max, x0(3:4)));
    U_ref(3:4) = max(-u_max,min(u_max, x0(3:4)));
    
    
    %% SCP Iteration
    for iter = 1:max_iter
        
        % 1. Linearize dynamics around x_ref
        % inside the SCP iteration, after x_ref exists
        [A,B,d] = linearize_dynamics(x_ref,N,dt,mu);
    
           %% --- Objective: fuel + squared-smoothness + total-variation (TV) on vector diffs
        epsilon = 1e-6;
        lambda_TV = 0.1;        % TV weight (promotes piecewise-constant thrusts => longer burns)
      
        fuel = 1;
        objective = @(U) ...        
            fuel* sum(sqrt(U(1:2:end).^2 + U(2:2:end).^2 + epsilon^2)) ...                   % fuel (approx l2)
            + lambda_TV * sum(sqrt((U(1:2:end-2)-U(3:2:end)).^2 + (U(2:2:end-2)-U(4:2:end)).^2 + epsilon^2)); %...
    
        lb = -u_max*ones(2*N,1);
        ub = u_max*ones(2*N,1);
        ternary = @(varargin) varargin{end - varargin{1}};
        options = optimoptions('fmincon', ...
            'Display','iter', ...
            'Algorithm','sqp', ...
            'MaxFunctionEvaluations',ternary(iter<5,1e4,1e5), ...
            'StepTolerance', 1e-6, ...
            'ConstraintTolerance', 1e-6);
    
        [U_opt, ~] = fmincon(objective, U_ref, [], [], [], [], lb, ub, ...
            @(U) dynamics_constraints_SCP(U, x0, xf, N, A, B,d,R_planet), options);
    
        U_opt(abs(U_opt) < 0.001) = 0;

        % 3. Reconstruct trajectory
           
        x = x0;
        traj = zeros(4,N+1);
        traj(:,1) = x0;
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            
            traj(:,k+1) = x;
        end
    
        % 4. Check convergence
        final_error = max(abs(traj(:,end) - xf));
        fprintf('Iteration %d, final state error = %.6f\n', iter, final_error);

        
        %must recalculate dynamics because our path has changed
        [A,B,d] = linearize_dynamics(traj(:,2:end),N,dt,mu); 
        
        x = x0;
        traj = zeros(4,N);
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            traj(:,k) = x;
        end
        
        final_error = max(abs(traj(:,end) - xf));
        fprintf('Iteration %d fter recalculating gravity, final state error = %.6f\n', iter, final_error);
        if final_error < tolerance
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
    function [c, ceq] = dynamics_constraints_SCP(U, x0, xf, N, A, B,d,R_planet)
        x = x0;
    
        traj = zeros(2,N);
        for k = 1:N
            uk = U(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            traj(:,k) = x(1:2);
          
        end
        ceq = x - xf;
        
        pos = traj; 
        dist2 = sum(pos.^2,1) ;      % squared distance from origin
        c = (R_planet^2 - dist2)';   % require dist^2 >= R^2 → c <= 0
    end
end