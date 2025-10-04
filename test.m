clc; clear; close all;

earth_r = 6378137 ;
earth_mu = 3.986004418e14;

% first part of orbital manuver
x0 = [0; earth_r*2; 5.591e3; 0];
xf = [earth_r*2; 10; 0;  4.61e3];


scp_min_fuel(x0,xf,earth_r,earth_mu,1.25);




% first part of orbital manuver
x0 = [0; earth_r*2; 5.61e3; 0];
xf = [0;-earth_r*2;  -5.61e3; 0];


%scp_min_fuel(x0,xf,earth_r,earth_mu,0.5);
