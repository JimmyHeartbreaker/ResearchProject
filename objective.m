 function J = objective(U,U_prev,jiter,N,reduce_jitter)
        epsilon = 1e-6;
        lambda_TV = 0.05;        % TV weight (promotes piecewise-constant thrusts => longer burns)
        
        fuel = 1;      
        dt_scale = U(end);
     %   dt_scale_prev = U_prev(end);
        U = U((N-1)*4+1:end-1);

       % if(reduce_jitter)
            U_prev = U_prev((N-1)*4+1:end-1);
            diffU = U - U_prev;
            J_iter = sum(diffU.^2)*jiter;
       % else
       %     J_iter = 0;
       % end
    delta = 0.01;
    alpha = 0.5;
    %small_thrusts = (abs(U) > 0) & (abs(U) < 0.01);
   % penalty = alpha * sum(exp(-abs(U)/delta) );
    
       J=   fuel* sum( sqrt(U(1:2:end).^2 + U(2:2:end).^2 + epsilon^2)) ...% + (dt_scale-dt_scale_prev)^2*1000;% ...%+  abs(dt_scale/dt_scale_prev - 1)*100;%   ...                            ; %
              + lambda_TV * sum(sqrt((U(1:2:end-2)-U(3:2:end)).^2 + (U(2:2:end-2)-U(4:2:end)).^2 + epsilon^2))...
              + J_iter...
              + dt_scale;% + penalty;   
    
    end