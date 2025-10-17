function [Gx,Gu,d] = linearize_dynamics(x_ref,u_ref,N,dt,mu)
    nx = 4; 
    nu = 2;
     %% Linearize dynamics around current trajectory
    Gx  = zeros(nx,nx,N);
    Gu  = zeros(nx,nu,N);
    d   = zeros(nx,1,N);
    f_i = zeros(nx,N);
    u_ref= reshape(u_ref(1:N*2), 2, N);

    for k = 1:N
        xk = x_ref(:,k);
        uk = u_ref(:,k);

        % --- Gravity ---
        pos = xk(1:2);
        r = max(norm(pos), 1e-12);
        g = -mu * pos / r^3;

        % --- Continuous dynamics ---
        vx = xk(3); vy = xk(4);
        f_i(:,k) = [vx; vy; uk(1)+g(1); uk(2)+g(2)];
       
        x_next  = xk + dt * f_i(:,k);

        % --- Linearized discrete-time matrices ---
        G = [-mu*(1/r^3 - 3*pos(1)^2/r^5), 3*mu*pos(1)*pos(2)/r^5;
             3*mu*pos(1)*pos(2)/r^5, -mu*(1/r^3 - 3*pos(2)^2/r^5)];
        A_c = [0 0 1 0;
               0 0 0 1;
               G(1,1) G(1,2) 0 0;
               G(2,1) G(2,2) 0 0];
        B_c = [0 0;
               0 0;
               1 0;
               0 1];

        Gx(:,:,k) = eye(nx) + dt*A_c;
        Gu(:,:,k) = dt*B_c;

        % --- Affine term ensures exact match ---
        d(:,:,k) = x_next - Gx(:,:,k)*xk - Gu(:,:,k)*uk;
    end
end