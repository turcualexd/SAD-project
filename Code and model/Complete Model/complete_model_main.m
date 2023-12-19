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
e    = 0.001830122216180;
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

omega0 = [0.01; 0.01; 0.01]*6;
s0     = zeros(3,1);
tol    = 0.2;  % tolerance for kinematic switch 312 - 313 and viceversa
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
FOV_sun    = deg2rad(85); % ??? to choose
freq_earth = 5; %Hertz to decrease
freq_sun   = 5; %Hertz to decrease
freq_mag   = 5;

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

%% SIMULATIONS 

%------------------------------Simulation----------------------------------

t0   = 0;
tf   = 500;
step = 0.05;
out  = sim("complete_model.slx", "Solver", "ode5", "StartTime", "t0", "StopTime", "tf", "SolverType", "Fixed-Step", "FixedStep", "step");

