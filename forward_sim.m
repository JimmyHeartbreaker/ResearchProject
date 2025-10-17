function x_sim = forward_sim(x0,u,dt_i,N)
    x_sim = zeros(4,N+1);
    x_sim(:,1) = x0;
    for k = 1:N
        pos = x_sim(1:2,k);
        r = max(norm(pos),1e-12);
        acc = -pos/r^3;          % gravity
        f = [x_sim(3,k); x_sim(4,k); u(1,k)+acc(1); u(2,k)+acc(2)];
        x_sim(:,k+1) = x_sim(:,k) + dt_i*f;  % forward Euler step
    end
end