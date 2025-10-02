%% Sequential Convex Programming (SCP) for 2D Minimum-Fuel STO
clc; clear; close all;

%% Problem setup
N = 50;    % discretization points


% Initial and final states: [x; y; vx; vy]
x0 = [0; 1000; 1.8; -0.3];
xf = [-300; -1000; -0.01;  0];
v0_hat = x0(3:4) / norm(x0(3:4));


% Gravity (central body at origin)
mu = 1000000;
% Maximum thrust
u_max = 3;
T = 100;

%% SCP Parameters
max_iter = 15;
tolerance = 5;

%% Initial reference trajectory: linear interpolatin
%trag = trag_estimate ( atan(x0(1)/x0(2)),sqrt(x0(1)^2 + x0(2)^2),sqrt(xf(1)^2 + xf(2)^2),N,mu); %zeros(4,N);
%x_ref = trag(:,2:5).';
weights = exp(-3*(0:N-1)/(N-1));  % exponential decay


x_ref = zeros(4,N); 
for k = 1:N 
    alpha = k/N; 
    x_ref(:,k) = x0*(1-alpha) + xf*alpha; 
end

epsilon = 1e-6; 

R_planet = 500;
dt = T/N;

pairs = reshape(x_ref, 2, []);  % each column is a pair

% take every 2nd pair
selected_pairs = pairs(:, 1:2:end);
diff = ( selected_pairs(:, 2:end) - selected_pairs(:, 1:end-1)) / dt;

%% Initial control guess: small random values
U_ref =[diff(:);[0;0]]; % x_ref(1:2:end)'% selected_pairs(:);% rand(2*N,1);%  (rand(2*N,1) - 0.5) * 1e-3;
%U_ref(1:2) = x0(3:4);


%% SCP Iteration
for iter = 1:max_iter
    
    % 1. Linearize dynamics around x_ref
    % inside the SCP iteration, after x_ref exists
    [A,B,d] = linearize_dynamics(x_ref,N,dt,mu);

       %% --- Objective: fuel + squared-smoothness + total-variation (TV) on vector diffs
    epsilon = 1e-6;
    lambda_smooth = 0;    % squared difference weight (smoothness)
    lambda_TV = 1;        % TV weight (promotes piecewise-constant thrusts => longer burns)
    lambda_vel = 2;
    fuel = 1;
    objective = @(U) ...        
        fuel* sum(sqrt(U(1:2:end).^2 + U(2:2:end).^2 + epsilon^2)) ...                   % fuel (approx l2)
        + lambda_TV * sum(sqrt((U(1:2:end-2)-U(3:2:end)).^2 + (U(2:2:end-2)-U(4:2:end)).^2 + epsilon^2)) ...
        + sum(lambda_vel * sum( ...
        weights .* (1 - ((v0_hat(1)*U(1:2:end) + v0_hat(2)*U(2:2:end)) ...
                         ./ sqrt(U(1:2:end).^2 + U(2:2:end).^2 + epsilon))).^2 ));


    lb = -u_max*ones(2*N,1);
    ub = u_max*ones(2*N,1);

    options = optimoptions('fmincon', ...
        'Display','iter', ...
        'Algorithm','sqp', ...
        'MaxFunctionEvaluations',1e6, ...
        'StepTolerance', 1e-5, ...
        'ConstraintTolerance', 1e-4);

    [U_opt, ~] = fmincon(objective, U_ref, [], [], [], [], lb, ub, ...
        @(U) dynamics_constraints_SCP(U, x0, xf, N, A, B,d,R_planet), options);

    U_opt(abs(U_opt) < 0.001) = 0;
    % 3. Reconstruct trajectory
   %    
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
    if final_error < tolerance
        

        [A,B,d] = linearize_dynamics(traj(:,2:end),N,dt,mu); 
        
        x = x0;
        traj = zeros(4,N);
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            traj(:,k) = x;
        end
        
        final_error = max(abs(traj(:,end) - xf));
        fprintf('Iteration %d, final state error = %.6f\n', iter, final_error);
        if final_error < tolerance
            fprintf('Converged!\n');
            break;
        end
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
    % Enforce initial state (already fixed by how we integrate) and final target
    % Instead of storing x in history, we just need final equality and optionally per-step defects.
    % To enforce dynamics exactly, better to stack defects between steps if x were an optimization variable.
    % But because we only optimize U, we enforce final state equality:
    ceq = x - xf;
    
    pos = traj; 
    dist2 = sum(pos.^2,1) ;      % squared distance from origin
    c = (R_planet^2 - dist2)';   % require dist^2 >= R^2 → c <= 0
end
