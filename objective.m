 function J = objective(U,U_prev,jiter)
        epsilon = 1e-6;
        lambda_TV = 1;        % TV weight (promotes piecewise-constant thrusts => longer burns)
        
        fuel = 1;      

        diffU = U - U_prev;
        J_iter = sum(diffU.^2)*jiter;

       J=   fuel* sum( sqrt(U(1:2:end).^2 + U(2:2:end).^2 + epsilon^2)) ...                           
            + lambda_TV * sum(sqrt((U(1:2:end-2)-U(3:2:end)).^2 + (U(2:2:end-2)-U(4:2:end)).^2 + epsilon^2)) ...      
            + J_iter;
    
    end