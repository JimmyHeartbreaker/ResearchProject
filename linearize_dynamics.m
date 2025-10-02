function [A,B,d] = linearize_dynamics(x_ref,N,dt,mu)

    A = zeros(4,4,N);
    B = zeros(4,2,N);
    d = zeros(4,N);   % affine term

    for k = 1:N
        pos = x_ref(1:2,k);
        r = max(norm(pos), 1e-3);            % avoid division by zero
        g_ref = -mu * pos / r^3;             % gravity at reference
        G = -mu*(eye(2)/r^3 - 3*(pos*pos')/r^5); % Jacobian (g_lin)
    
        A(:,:,k) = [eye(2), dt*eye(2);
                    dt*G,  eye(2)];
        B(:,:,k) = [zeros(2,2); dt*eye(2)];
    
        % affine term: only on velocity rows
        d(1:2,k) = [0;0];   % no constant term for pos-equation
        d(3:4,k) = dt*(g_ref - G*pos);   % important affine correction
    end
end