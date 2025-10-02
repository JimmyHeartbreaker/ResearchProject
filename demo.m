clc; clear; close all;

earth_r = 6378e3;
earth_mu = 3.986004418e14;

% first part of orbital manuver
x0 = [0; earth_r*2; 4.61e3; 0];
xf = [-earth_r*2; -earth_r*3; -5.61e3;  0];


  scp_min_fuel(x0,xf,earth_r,earth_mu,0.7);


%second part
x0 = [-earth_r*2; -earth_r*3; -5.61e3;  0];
xf = [earth_r*2; earth_r*3; 5.61e3;  0];



  scp_min_fuel(x0,xf,earth_r,earth_mu,0.5);