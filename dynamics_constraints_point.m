function [c, ceq] = dynamics_constraints_point(U, x0, xf, N, A, B,d)
        traj = zeros(2,N+1);  
        
        x = x0;   
        traj(:,1) = x0(1:2);
        U = reshape(U, 2, N);

        for k = 1:N
            x = A(:,:,k)*x + B(:,:,k)*U(:,k) + d(:,k);
            traj(:,k+1) = x(1:2);
        end
        
        c= 1 - vecnorm(traj);
        
        ceq = xf-x;
      %  ceq
    end