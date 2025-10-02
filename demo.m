clc; clear; close all;

earth_r = 500; %6378e3
earth_mu = 1000000; %3.986004418e14
y_pos = 10;
max_speed= 3;
x0 = [0; 1000; 1.8; -0.3];
xf = [-300; -1000; -0.01;  0];

%  x0 = [0; earth_r+y_pos; 0.001; 0];
 % xf = [0; -earth_r-y_pos;-0.001; 0];

  scp_min_fuel(x0,xf,earth_r,earth_mu,max_speed);