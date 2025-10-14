function x_sim = forward_sim_RK4(x0,u,dt_i,N)
    x_sim = zeros(4,N+1);
    x_sim(:,1) = x0;
    mu=1;
    for k = 1:N
    
        %% --- Dynamics function ---
        f = @(x,u) [ ...
            x(3); ...
            x(4); ...
            u(1,k) - mu*x(1)/norm(x(1:2))^3; ...
            u(2,k) - mu*x(2)/norm(x(1:2))^3 ...
        ];
    
        %% --- RK4 integration ---
        k1 = f(x_sim(:,k), u);
        k2 = f(x_sim(:,k) + 0.5*dt_i*k1, u);
        k3 = f(x_sim(:,k) + 0.5*dt_i*k2, u);
        k4 = f(x_sim(:,k) + dt_i*k3, u);
    
        x_sim(:,k+1) = x_sim(:,k) + (dt_i/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end