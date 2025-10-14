function [c, ceq] = dynamics_constraints_point_variable_time(z, x0, xf, N, Gx,Gu,Gdt,d,dt_nom)
        nx = 4; 
        nu = 2;
        x = reshape(z(1:nx*(N-1)), nx, N-1);
        u = reshape(z(nx*(N-1)+1:nx*(N-1)+nu*N), nu, N);
        dt_scaler = z(end);
    
        ceq =  zeros(4,N);  
        I = eye(nx);
        
        Gx_eff = I + (Gx - I) * dt_scaler;
       Gu_u_plus_d = dt_scaler * (squeeze(pagemtimes(Gu, reshape(u, [size(u,1), 1, size(u,2)]))) + ...
                           reshape(d, size(d,1), size(d,3)));

       x_prev = [x0, x(:,1:end-1)];                  % nx × (N-1)
        x_prev3 = reshape(x_prev, [size(x_prev,1),1,size(x_prev,2)]); % nx × 1 × (N-1)
        
        % --- Step 4: Multiply all Gx_eff(:,:,k) * x_prev(:,k) in one batch ---
        Gx_xprev = squeeze(pagemtimes(Gx_eff(:,:,1:N-1), x_prev3)); % nx × (N-1)
        
        % --- Step 5: Compute ceq for middle steps k = 2:N-1 ---
        
        ceq(:,2:N-1) = x(:,2:N-1) - (Gx_xprev(:,2:N-1) + Gu_u_plus_d(:,2:N-1));

ceq(:,1) = x(:,1) - (Gx_eff(:,:,1)*x0 + Gu_u_plus_d(:,1));
ceq(:,N) = xf - (Gx_eff(:,:,N)*x(:,N-1) + Gu_u_plus_d(:,N));

     
        ceq = ceq(:);
        c= [];%1 - vecnorm(traj);
        
end


% function [c, ceq] = dynamics_constraints_point_variable_time(z, x0, xf, N, Gx,Gu,Gdt,d,dt_nom)
%         nx = 4; 
%         nu = 2;
%         x = reshape(z(1:nx*(N-1)), nx, N-1);
%         u = reshape(z(nx*(N-1)+1:nx*(N-1)+nu*N), nu, N);
%         dt_scaler = z(end);
% 
%         ceq =  zeros(4,N);  
%         I = eye(nx);
%        Gx_eff = I + (Gx(:,:,1) - I) * dt_scaler;
%        ceq(:,1) =x(:,1) - (Gx_eff*x0 + (dt_scaler)*(Gu(:,:,1)*u(:,1) +  d(:,:,1)));
% 
%         Gx_eff = I + (Gx - I) * dt_scaler;
%         for k = 2:N-1
%             ceq(:,k) =x(:,k) - (Gx_eff(:,:,k)*x(:,k-1) + (dt_scaler)*(Gu(:,:,k)*u(:,k) +  d(:,:,k)));
%         end
% 
%        Gx_eff = I + (Gx(:,:,N) - I) * dt_scaler;
%        ceq(:,N) =xf - (Gx_eff*x(:,N-1) + (dt_scaler)*(Gu(:,:,N)*u(:,N) +  d(:,:,N)));
% 
%         ceq = ceq(:);
%         c= [];%1 - vecnorm(traj);
% 
%     end