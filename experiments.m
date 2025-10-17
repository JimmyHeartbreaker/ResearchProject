clc; clear; close all;

earth_r = 6378137 ;
earth_mu = 3.986004418e14;


for N = 4:5
    % % %test 1/2 orbit
    x0 = [0; earth_r*2; 5.59e3; 0];
    xf = [0;-earth_r*N;-5.59e3; 0];
    [final_vector, dt_guess, final_dt_scale] = scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);
N
dt_guess
final_dt_scale
end