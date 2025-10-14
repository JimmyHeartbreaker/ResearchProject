
%% Sequential Convex Programming (SCP) for 2D Minimum-Fuel STO

function final_vector = scp_min_fuel(x0,xf,planet_radius,mu_real,extra_periods,N) 
    u = x0(1:2)';
    v = xf(1:2)';
    theta = atan2(det([v; u]), dot(v,u));
    if theta < ((extra_periods) * 360 * pi/180) + 1e-6
        theta = theta + 2*pi;
    end
     u2 = u + x0(3:4)';
    theta2 = atan2(det([u2; u]), dot(u2,u));
    if(theta2<0)
        theta = theta - pi*2;
    end

    periods = abs(theta / ( pi*2));
   
    mu=1; %variables are scaled so that mu and R == 1
    
    T0 = sqrt(planet_radius^3/mu_real);
    Tdim = periods * 2 * pi * sqrt(sqrt(x0(1)^2 + x0(2)^2)^3 / mu_real);
    x0(1:2) = x0(1:2) / planet_radius; 
    xf(1:2) = xf(1:2) / planet_radius;

    T = Tdim / T0;
    dt_scale = 1.2;
    ds = 1/N;
    a0 = planet_radius / T0^2;
    v0 = planet_radius / T0;
    x0(3:4) = x0(3:4) / v0; 
    xf(3:4) = xf(3:4) / v0;
    % Maximum thrust
    u_max = (5/a0);
    t_f = T;
    R_planet = 1;
    dt =  ds * t_f;
    dt_nom = dt;
    r = norm(x0(1:2));   % distance from planet center  
   
    t_hat = [x0(2); x0(1)] / r;
    v_orb = sqrt(mu/r) * t_hat * sign(theta2) ;  % desired orbital velocity
    dv = v_orb - x0(3:4);              % velocity delta to align with orbit
    
    U_mat = zeros(2,N) * 1e-6;
    U_mat(1:2) = dv/dt;
 
    x_ref = zeros(4,N+1);
    
    theta0 = atan(x0(1)/x0(2));

    if(x0(2)<0)
        theta0 = theta0 + pi;
    end
    for k = 1:(N+1)
        alpha = (k-1)/(N);
        y = cos((theta ) * alpha+ theta0) * r;
        x = sin((theta) * alpha + theta0) * r;
        x_ref(1:2,k) = [x;y];
        vec = [y,-x];

        x_ref(3:4,k) = (vec/norm(vec))  * norm(v_orb);
    end

    for k = 1:N

        x_ref(3:4,k) = (x_ref(1:2,k+1) - x_ref(1:2,k)) / dt;
    end
    
    x_ref(3:4,N+1) = x_ref(3:4,N);
    
    U_ref = reshape(U_mat, 2*N, 1);
    x_ref_flat = x_ref(:);
    z0 = [x_ref_flat(5:end-4); U_ref(:); dt_scale];
    
    [Gx,Gu,d] = linearize_dynamics(x_ref,U_ref,N,dt_nom,mu);
    max_iter = 50;
    x_max = 10;
    lb = [];
    ub = [];

    for k = 1:(N-1)*4
        lb = [lb; -x_max];
        ub = [ub; x_max];
    end
    
    for k = 1:N*2
        lb = [lb; -u_max];
        ub = [ub; u_max];
    end
    dt_min = 0.9;
    dt_max = 1.5;
    % Time scale bounds
    lb = [lb; dt_min];
    ub = [ub; dt_max];
    rk4_mode = 0;
    linearize_dynamics_fn = @linearize_dynamics;
    forward_sim_fn = @forward_sim;
    %% SCP Iteration
    for iter = 1:max_iter         
    
        dt =  ds * t_f;
    
      
                jiter_control = abs(theta) *6;
             obj = @(U) objective(U,z0,jiter_control,N,rk4_mode);
           options = optimoptions('fmincon', ...            
               'Display','iter', ...
            'Algorithm','sqp', ...
            'MaxFunctionEvaluations',2e5, ...  %            'FiniteDifferenceType','central', ...
             'StepTolerance',  1/(10^(min(8,iter+2))), ...    
    'OptimalityTolerance', 2, ...
    'ConstraintTolerance', 1/(10^(min(7,iter+1))));
            [U_opt, fval, exitflag, output, lambda] = fmincon(obj, z0, [], [], [], [], lb, ub, ...
                @(U) dynamics_constraints(U, x0, xf, N, Gx,Gu,d), options);
     
       
        dt_range_scaler = sqrt(U_opt(end) / dt_scale);
        dt_scale = U_opt(end);
        if(dt_range_scaler > 1)
            dt_min = dt_scale / dt_range_scaler;
            dt_max = dt_scale *dt_range_scaler;
        else
            dt_max = dt_scale / dt_range_scaler;
            dt_min = dt_scale *dt_range_scaler;
        end
        ub(end) = dt_max;
        lb(end) = dt_min;

        u = reshape(U_opt(4*(N-1)+1 : 4*(N-1)+2*N), 2, N);
        x= reshape( [x0;U_opt(1: 4*(N-1));xf], 4, N+1);
      
      
        [Gx,Gu,d] = linearize_dynamics_fn(x,u,N,dt_nom,mu);
        x_sim = forward_sim_fn(x0,u,dt_nom * dt_scale,N);
        x_ref_flat = x(:);
        z0 = [x_ref_flat(5:end-4); u(:); dt_scale];
        plot(x(1,:), x(2,:), 'o-',x_sim(1,:), x_sim(2,:), 'x-'); axis equal;
        legend('opt states','forward-sim'); title('Trajectory: optimized vs forward-sim');

    
            
        x_ref=x_sim;

        v = x_ref(3:4,N+1) * (dt/1.5);
        r = xf(1:2) - x_ref(1:2,N+1) ;  
        a = r;
        b = v;
        a_proj_onto_b = (dot(a,b) / norm(b)^2) * b; 
        perp = a - a_proj_onto_b;
        c_perp =norm(perp);
        c_along = norm(a_proj_onto_b) -  norm(v);     
        vector_error =norm(v / norm(v)- xf(3:4)/ norm(xf(3:4)));
        if rk4_mode==0 && norm(x_ref(:,N+1) - xf) < 0.1
            linearize_dynamics_fn = @linearize_dynamics_RK4;
            forward_sim_fn = @forward_sim_RK4;
            dt_nom = dt_nom * dt_scale;
            dt_scale = 1;
            ub(end) = 1;
            lb(end) = 1;
            rk4_mode = 1;
            [Gx,Gu,d] = linearize_dynamics_fn(x,u,N,dt_nom,mu);
        end
        
        fprintf('Iteration %d, c_along = %.6f, c_perp = %.6f, vector_error = %.6f %.6f\n', iter, c_along,c_perp,vector_error);
        if c_along < 0 && c_perp <  1e-3  &&  vector_error < 1e-3  
            fprintf('Converged!  iter = %d\n',iter);
            U_ref = u;
            U_opt = u;
            break;             
        end
        
     
        U_ref = U_opt;
            
    end
    
    plot_control(T,N,U_opt);
    plot_traj(x_ref,U_ref(:),x0,xf,R_planet);
    R_scale = 6378e3;

    final_vector = x_ref(:,end);
    final_vector(1:2) = final_vector(1:2) * R_scale; 
    final_vector(3:4) = final_vector(3:4) * v0;
end