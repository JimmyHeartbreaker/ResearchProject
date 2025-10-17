clc; clear; close all;

earth_r = 6378e3;
earth_mu = 3.986004418e14;

% first part of orbital manuver
x0 = [0; earth_r*2; 5.592e3; 0];
xf = [0; -earth_r*11; 2.38e3;  0];


final_vector = scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);


%second part
% x0 = final_vector;
% xf = [0; earth_r*3; 4.56e3;  0];
% 
% 
% 
% scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);