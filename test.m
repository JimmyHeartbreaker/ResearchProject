clc; clear; close all;

earth_r = 6378137 ;
earth_mu = 3.986004418e14;

%test 1.25 orbit
x0 = [0; earth_r*2; 5.58e3; 0];
xf = [earth_r*3; 0; 0; -4.66e3];
scp_min_fuel(x0,xf,earth_r,earth_mu,1,100);

%test 1/4 orbit
  x0 = [0; earth_r*2; -5e3; 0];
  xf = [-earth_r*3; 0; 0;-5e3];
  scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);

  %test 1/4 orbit, far
  x0 = [0; earth_r*10; -5e3; 0];
  xf = [-earth_r*3; -earth_r*3; 0;-5e3];
  scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);

%test 1/2 orbit
x0 = [0; earth_r*2; 5e3; 0];
xf = [0; -earth_r*2; -5e3; 0];
scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);

%test 1/2 orbit far
x0 = [0; earth_r*2; 5e3; 0];
xf = [0; -earth_r*6; -5e3; 0];
scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);


%test 3/4 orbit
x0 = [0; earth_r*2; 5e3; 0];
xf = [-earth_r*3; 0; 0;5e3];
scp_min_fuel(x0,xf,earth_r,earth_mu,0,100);


%test 1 orbit
x0 = [0; earth_r*2; 5.58e3; 0];
xf = [0; earth_r*3; 4.66e3; 0];
scp_min_fuel(x0,xf,earth_r,earth_mu,1,100);





