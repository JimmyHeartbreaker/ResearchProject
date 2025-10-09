
%% Sequential Convex Programming (SCP) for 2D Minimum-Fuel STO

function scp_min_fuel(x0,xf,planet_radius,mu_real,extra_periods,N) 
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
    a0 = planet_radius / T0^2;
    v0 = planet_radius / T0;
    x0(3:4) = x0(3:4) / v0; 
    xf(3:4) = xf(3:4) / v0;
    % Maximum thrust
    u_max = (5/a0);

    R_planet = 1;
    dt =  T/N;

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
        alpha = (k-1)/(N+1);
        y = cos((theta ) * alpha+ theta0) * r;
        x = sin((theta) * alpha + theta0) * r;
        x_ref(1:2,k) = [x;y];
    end

    for k = 1:N
        x_ref(3:4,k) = (x_ref(1:2,k+1) - x_ref(1:2,k)) / dt;
    end
    
    x_ref(3:4,N+1) = x_ref(3:4,N);
    U_ref = reshape(U_mat, 2*N, 1);
  

    [A,B,d] = linearize_dynamics(x_ref,N,dt,mu);
    max_iter = 50;
    line_search=0;
    %% SCP Iteration
    for iter = 1:max_iter         
    
        lb = -u_max*ones(2*N,1);
        ub = u_max*ones(2*N,1);
    
        if(line_search == 0)
                jiter_control = abs(theta) * 3;
             obj = @(U) objective(U,U_ref,jiter_control);
            options = optimoptions('fmincon', ...            'Display','iter', ...%
            'Algorithm','sqp', ...
            'MaxFunctionEvaluations',2e5, ...  % 'FiniteDifferenceType','central', ...
            'StepTolerance', 1e-8, ...
            'ConstraintTolerance',1e-4, ...
            'OptimalityTolerance',1.5); 
            [U_opt, fval, exitflag, output, lambda] = fmincon(obj, U_ref, [], [], [], [], lb, ub, ...
                @(U) dynamics_constraints_point(U, x0, xf, N, A, B,d), options);
        else
           
                jiter_control = abs(theta) * 6;

            obj = @(U) objective(U,U_ref,jiter_control);
            options = optimoptions('fmincon', ...%            'Display','iter', ...
            'Algorithm','sqp', ...
            'MaxFunctionEvaluations',2e5, ...  % 'FiniteDifferenceType','central', ...
            'StepTolerance', 1e-8, ...
            'ConstraintTolerance',1e-5, ...
            'OptimalityTolerance',1.5); 
            [U_opt, fval, exitflag, output, lambda] = fmincon(obj, U_ref, [], [], [], [], lb, ub, ...
                @(U) dynamics_constraints_line(U, x0, xf, N, A, B,d,dt), options);
        end

        % 3. Reconstruct trajectory
           
        x = x0;
        traj = zeros(4,N+1);
        traj(:,1) = x0;
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A(:,:,k)*x + B(:,:,k)*uk + d(:,k);
            
            traj(:,k+1) = x;
        end
    
         %must recalculate dynamics because our path has changed
        
       [A2,B2,d2] = linearize_dynamics(traj,N,dt,mu); 
        diffA =abs( sum(sum(sum(A2-A))));
        if(diffA < 0.005)
             matAlpha = 0.9;
        else
            matAlpha = 0.8;
        end
        A_opt = (A2*matAlpha + A * (1-matAlpha)) ; %bisect matrices to stop oscillation
        B_opt= (B2*matAlpha + B * (1-matAlpha)) ;
        d_opt= (d2*matAlpha + d * (1-matAlpha)) ;
        x = x0;
        traj = zeros(4,N+1);
        traj(:,1) = x0;
        for k = 1:N
            uk = U_opt(2*k-1:2*k);
            x = A2(:,:,k)*x + B2(:,:,k)*uk + d2(:,k);
            
            traj(:,k+1) = x;
        end
        
        total_error =norm(traj(:,end) - xf(:));
        if(total_error< dt * norm(xf(3:4)))
            line_search = 1;
        elseif total_error > 0.5
            line_search = 0;
        end

       
        A = A_opt;
        B = B_opt;
        d = d_opt;
            
        x_ref=traj;

        v = x_ref(3:4,N+1) * (dt/1.5);
        r = xf(1:2) - x_ref(1:2,N+1) ;  
        a = r;
        b = v;
        a_proj_onto_b = (dot(a,b) / norm(b)^2) * b; 
        perp = a - a_proj_onto_b;
        c_perp =norm(perp);
        c_along = norm(a_proj_onto_b) -  norm(v);     
        vector_error =norm(v / norm(v)- xf(3:4)/ norm(xf(3:4)));
        fprintf('Iteration %d, c_along = %.6f, c_perp = %.6f, vector_error = %.6f %.6f\n', iter, c_along,c_perp,vector_error);
        if c_along < 0 && c_perp <  1e-3  &&  vector_error < 1e-3  
            fprintf('Converged!  iter = %d\n',iter);
            break;             
        end
        
     
        U_ref = U_opt;
            
    end
    
    plot_control(T,N,U_opt);
    plot_traj(x_ref,x0,xf,R_planet);
    
 
end