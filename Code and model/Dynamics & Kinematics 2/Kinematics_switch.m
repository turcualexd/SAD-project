clear, clc, close all;

%% Initial conditions

% Constants
mu = 3.98600433e+5;

% Satellite data
Ix = 7776999842.85e-9;              % kg*m^2
Iy = 22972402681.94e-9;             % kg*m^2
Iz = 23189503421.86e-9;             % kg*m^2
I = diag([Ix Iy Iz]);
I_inv = inv(I);
w0 = [0.5; -0.2; 1.0];              % rad/s
s0 = [0; asin(0.4); pi/3];          % rad

% Orbit data
a = 6.851221970532768e3;            % km
e = 0.001830122216180;              % rad
i = 1.699980862034725;              % rad
OM = 0.554268509489784;             % rad
om = 1.772848103192913;             % rad
theta0 = 3.109851878856139;         % rad
n = sqrt(mu / a^3);                 % rad/s

%% Costants calculations

R_OM =  [ cos(OM),    sin(OM),    0;
          -sin(OM),   cos(OM),    0;
          0,          0,          1 ];

R_i =   [ 1,          0,          0;
          0,          cos(i),     sin(i);
          0,          -sin(i),    cos(i) ];

R_om =  [ cos(om),    sin(om),    0;
          -sin(om),   cos(om),    0;
          0,          0,          1 ];

A_313 = ( R_om * R_i * R_OM ).';    % from perifocal to inertial

%% Simulation

t0 = 0;
tf = 200;
step = 0.1;
tol = 0.2;

out = sim("Kinematics_switch_sim.slx", "Solver", "ode5", "StartTime", "t0", "StopTime", "tf", "SolverType", "Fixed-Step", "FixedStep", "step");

%% Output

t = out.tout;
w = out.w.data;
A = out.A.data;