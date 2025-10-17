function [Gx,Gu,d] = linearize_dynamics_RK4(x_ref,u_ref,N,dt,mu)
   nx = 4; nu = 2;
Gx = zeros(nx,nx,N);
Gu = zeros(nx,nu,N);
d  = zeros(nx,1,N);
f_i = zeros(nx,N);
u_ref = reshape(u_ref(1:N*2),2,N);

    for k = 1:N
    xk = x_ref(:,k);
        uk = u_ref(:,k);
    
        % --- Dynamics function ---
        f = @(x,u) [x(3); x(4); u(1)-mu*x(1)/max(norm(x(1:2)),1e-12)^3; u(2)-mu*x(2)/max(norm(x(1:2)),1e-12)^3];
    
        % --- Continuous-time Jacobians ---
        A_fun = @(x) [0 0 1 0; 0 0 0 1; 
            -mu*(1/norm(x(1:2))^3 - 3*x(1)^2/norm(x(1:2))^5), 3*mu*x(1)*x(2)/norm(x(1:2))^5, 0, 0;
             3*mu*x(1)*x(2)/norm(x(1:2))^5, -mu*(1/norm(x(1:2))^3 - 3*x(2)^2/norm(x(1:2))^5), 0, 0];
        B_fun = [0 0; 0 0; 1 0; 0 1];
    
        % --- RK4 step ---
        k1 = f(xk,uk);
        k2 = f(xk + 0.5*dt*k1, uk);
        k3 = f(xk + 0.5*dt*k2, uk);
        k4 = f(xk + dt*k3, uk);
        x_next = xk + dt/6*(k1 + 2*k2 + 2*k3 + k4);
        f_i(:,k) = k1;
    
        % --- RK4 sensitivities ---
        A1 = A_fun(xk);          B1 = B_fun;
        Kx1 = A1;                Ku1 = B1;
        A2 = A_fun(xk + 0.5*dt*k1); Kx2 = A2*(eye(nx)+0.5*dt*Kx1); Ku2 = B1 + 0.5*dt*A2*Ku1;
        A3 = A_fun(xk + 0.5*dt*k2); Kx3 = A3*(eye(nx)+0.5*dt*Kx2); Ku3 = B1 + 0.5*dt*A3*Ku2;
        A4 = A_fun(xk + dt*k3);      Kx4 = A4*(eye(nx)+dt*Kx3);   Ku4 = B1 + dt*A4*Ku3;
    
        Gx(:,:,k) = eye(nx) + dt/6*(Kx1 + 2*Kx2 + 2*Kx3 + Kx4);
        Gu(:,:,k) = dt/6*(Ku1 + 2*Ku2 + 2*Ku3 + Ku4);
    
        % --- Affine term ensures exact match ---
        d(:,:,k) = x_next - Gx(:,:,k)*xk - Gu(:,:,k)*uk;
    end
end