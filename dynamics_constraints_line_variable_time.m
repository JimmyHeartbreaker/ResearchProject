function [c, ceq] = dynamics_constraints_line(U, x0, xf, N, A, B,d,dt)
        % persistent i
        % if( isempty(i))
        %     i= 0;
        % end
        % i = i + 1;
        
        traj = zeros(2,N+1);  
        
        x = x0;   
        traj(:,1) = x0(1:2);
        U = reshape(U, 2, N);

        for k = 1:N
            x = A(:,:,k)*x + B(:,:,k)*U(:,k) + d(:,k);
            traj(:,k+1) = x(1:2);
         
        end

         v =x(3:4)*(dt/2); % line direction
            r = xf(1:2) - traj(1:2,N+1) ;  
            a = r;
            b = v;
            a_proj_onto_b = (dot(a,b) / norm(b)^2) * b; 
            perp = a - a_proj_onto_b;
            c_perp =norm(perp);
            c_along = norm(a_proj_onto_b) -  norm(v);     

            vector_error =norm(v / norm(v)- xf(3:4)/ norm(xf(3:4)));
       

        % if mod(i,100000)==0
        %     plot_debug(traj,x0,xf,1);
        % end
        in_planet = -1;
        if(any(1 - vecnorm(traj) >= 0))
            in_planet = 1;
        end
        c= [in_planet;c_along];
        
        ceq = [c_perp;vector_error];
      %  ceq
    end