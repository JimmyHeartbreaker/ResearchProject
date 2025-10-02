function nodes = trag_estimate(theta0,r1, r2,N,mu)
    theta0 = theta0 - 90 * pi/180;
    
    %% Transfer orbit parameters
    a_t = (r1 + r2)/2;         % semi-major axis
    e = (r2 - r1)/(r2 + r1);   % eccentricity
    T_transfer = pi * sqrt(a_t^3 / mu); % half-period
    
    %% Time nodes
    t_nodes = linspace(0, T_transfer, N);
    
    %% Function to solve Kepler's equation: M = E - e*sin(E)
    function E = kepler_E(M, e)
        tol = 1e-10;
        max_iter = 100;
        if e < 0.8
            E = M;
        else
            E = pi;
        end
        for k = 1:max_iter
            f = E - e*sin(E) - M;
            f_prime = 1 - e*cos(E);
            dE = -f/f_prime;
            E = E + dE;
            if abs(dE) < tol
                break;
            end
        end
    end
    
    %% Compute nodes
    nodes = zeros(N,5); % columns: t, x, y, vx, vy
    
    for i = 1:N
        t = t_nodes(i);
        M = sqrt(mu / a_t^3) * t;        % mean anomaly
        E = kepler_E(M, e);              % eccentric anomaly
        % True anomaly
        theta = 2*atan(sqrt((1+e)/(1-e))*tan(E/2)); %+ theta0;
        % Radius
        r = a_t*(1 - e^2)/(1 + e*cos(theta));
        % Position
        x = r*cos(theta+theta0);
        y = r*sin(theta+theta0);
        % Orbital speed (vis-viva)
        v = sqrt(mu*(2/r - 1/a_t));
        % Velocity components in orbital plane
        vx = -v*sin(theta+theta0);
        vy = v*(e + cos(theta+theta0)) / (1 + e*cos(theta+theta0));
        
        nodes(i,:) = [t, x, y, vx, vy];
    end
    
    %% Display results
    fprintf('t[s]      x[km]      y[km]      vx[km/s]    vy[km/s]\n');
    for i = 1:N
        fprintf('%8.1f %10.1f %10.1f %10.3f %10.3f\n', nodes(i,:));
    end
end