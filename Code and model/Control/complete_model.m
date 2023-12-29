clear; close all; clc;

%------------------Space-Craft Geometry and Properties---------------------

I1 = 22972402681.94*1e-9;
I2 = 7776999842.85*1e-9;
I3 = 23189503421.86*1e-9;

I = diag([I1 I2 I3]);
Iinv = I\eye(3);

%
% TO ADD UNDEPLOYED CONFIGURATION
%

%-----------------------------Orbit properties-----------------------------

mu   = astroConstants(13);
a    = 6.851221970532768e+03;
%e    = 0.001830122216180;
e = 0;
i    = 1.699980862034725;
om   = 1.772848103192913;
OM   = 0.554268509489784;
th_0 = 3.109851878856139;
n    = sqrt(mu/a^3);
T    = 2*pi/n;

%--------------------------Position of Earth wrt Sun-----------------------

mu_sun = astroConstants(4);
e_E = 1.743742542622805E-02;
i_E = deg2rad(2.343620856296279E+01);
OM_E = deg2rad(3.599997935134700E+02);
om_E = deg2rad(1.030408771388982E+02);
th_E = deg2rad(3.403771999857636E+02);
a_E = 1.497081833259778E+08;

[rr, vv] = kep2car(a_E, e_E, i_E, OM_E, om_E, th_E, mu_sun);
rs = -rr;

%------------Initial Conditions % constant values for Simulation-----------

%--------------------Dynamics and Kinematics subsystems--------------------

omega0  = 0.06 * rand(3,1);

s0     = 2*pi*rand(3,1);

tol     = 0.2;  % tolerance for kinematic switch 312 - 313 and viceversa
%---------------------Keplerian Dynamic Sub-system-------------------------

R_OM =  [ cos(OM),    sin(OM),    0;
          -sin(OM),   cos(OM),    0;
          0,          0,          1 ];

R_i =   [ 1,          0,          0;
          0,          cos(i),     sin(i);
          0,          -sin(i),    cos(i) ];

R_om =  [ cos(om),    sin(om),    0;
          -sin(om),   cos(om),    0;
          0,          0,          1 ];


A_pn  = R_om * R_i * R_OM; %from inertial to perifocal

%-------------Magnetic Field Modeling and distrubance dipole j-------------

g01    = -29404.8; g11 = -1450.9; h11 = 4652.5;
H0     = sqrt(g01^2+g11^2+h11^2)*1e-9;
jb     = [0.01; 0.05; 0.01];
incl   = deg2rad(11.5);
omegaE = deg2rad(15)/3600;
Rt     = astroConstants(23);

%-------------------------------Sensors------------------------------------

FOV_earth  = deg2rad(33); % MEISEI ELECTRIC 
FOV_sun    = deg2rad(100); % ??? to choose
freq_earth = 30; 
freq_sun   = 30; 

% Magnetometer AAC ClydeSpace MM200
freq_mag     = 30;         % Can be modulated up to 500Hz
noise_sd_mag = 1.18*1e-9;  % T/sqrt(Hz) 
bias_mag     = rand(3,1);
bias_mag     = bias_mag./norm(bias_mag);

%--------------------------------------------------------------------------

freq_acds = 30;

%---------------------------Att. Determination-----------------------------

alpha1 = 0.4;
alpha2 = 0.1;
alpha3 = 0.5;

%-------------------------------Actuators----------------------------------

max_dip   =  150; %[A m^2] NMTR-X CUSTOM
freq_act  = 5;

%-----------------------------De-Tumbling----------------------------------

chsi = (i - deg2rad(11));
k_gain  = (4*pi/T) * (1 + sin(chsi)) * I2;

%-----------------------------Control--------------------------------------
Kx = (I3 - I2) /  I1; Ky = (I3 - I1) / I2;

A = [
    0 (1-Kx)*n 0 -n^2*Kx 0 0;
    (Ky-1)*n 0 0 0 -n^2*Ky 0;
    0 0 0 0 0 0;
    1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0];

B = [
    1/I1 0 0;
    0 1/I2 0;
    0 0 1/I3;
    0 0 0;
    0 0 0;
    0 0 0;];

C = [
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1];

D = zeros(3);

sys = ss(A, B, C, D);

Co = ctrb(A,B);

Ob = obsv(A,C);

% Kp = 1e-2*diag([0.4;0.1;0.56]);
% Kd = 1e-4*diag([0.1;0.05;0.105]);
% 
% K = [Kd Kp];
% 
p = [-0.02+0.1i;
    -0.02-0.1i;
    -0.01+0.06i;
    -0.01-0.06i;
    -0.05+0.1i;
    -0.05-0.1i];

K = place(A,B,p);

eig(A)

eig(A - B * K)