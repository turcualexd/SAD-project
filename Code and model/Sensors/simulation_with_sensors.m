clear; close all; clc; 

tol = 1e-2;

% sat data
I = diag([100.09 25.1 91.6]);
I_inv = inv(I);
w0 = [0.2; 0; 0.2];
s0 = [0; 0; 0];

% orbit data
OM = 0;
om = 0;
i = 0.9601;
e = 2.3103*1e-4;
a = 2.6584*1e4;
n = sqrt(astroConstants(13) / a^3);

R_OM =  [ cos(OM),    sin(OM),    0;
          -sin(OM),   cos(OM),    0;
          0,          0,          1 ];

R_i =   [ 1,          0,          0;
          0,          cos(i),     sin(i);
          0,          -sin(i),    cos(i) ];

R_om =  [ cos(om),    sin(om),    0;
          -sin(om),   cos(om),    0;
          0,          0,          1 ];

A_313 = ( R_om * R_i * R_OM )'; %from perifocal to inertial


% magnetic field

g01 = -29404.8; g11 = -1450.9; h11 = 4652.5;
H0 = sqrt(g01^2+g11^2+h11^2)*1e-9;
incl = deg2rad(11.5);
omegaE = deg2rad(15)/3600;




