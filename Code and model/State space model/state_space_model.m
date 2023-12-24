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


%--------------------Dynamics and Kinematics subsystems--------------------

omega0 = [0; 0; n];
s0     = [2.099, -0.2003, -1.703];
tol    = 0.2;  % tolerance for kinematic switch 312 - 313 and viceversa


%---------------------Linearized Dynamics Matrices-------------------------

ky = (I3 - I2)/I1;
kr = (I3 - I1)/I2;
kp = (I2 - I1)/I3;

A11 = [0 (1 - ky)*n 0; (kr-1)*n 0 0; 0 0 0];
A12 = diag([-ky*n^2 -4*kr*n^2 -3*kp*n^2 ]);
A21 = eye(3);

A   = zeros(6,6);
B   = zeros(6,3);

A(1:3,1:3)  = A11;
A(4:6, 1:3) = A21;
A(1:3, 4:6) = A12;

B(1:3,1:3)   = diag([1/I1 1/I2 1/I3]);

C   = [[0 0 0 1 0 0];[0 0 0 0 1 0]; [0 0 0 0 0 1]];

D   = zeros(3,3);


sys = ss(A, B, C, D);


Q = diag([1/deg2rad(3)^2 1/deg2rad(3)^2 1/deg2rad(3)^2 1/n^2 1/n^2 1/n^2]);
R = diag(100*ones(3,1));
N =  zeros(6,3);
[K,S,P] = lqr(sys,Q,R,N);


poles_cl = eig(A - B*K);
obs_poles = (min(real(poles_cl))*3)-(1:6);

Ltr = place(A', C', obs_poles);
L = Ltr';

obs = ss(A-L*C,[B L],C,0);
%bodeplot(obs)

%extended state space
A_ext11 = A - B*K;
A_ext12 = B*K;
A_ext21 = zeros(6,6);
A_ext22 = A - L*C;
A_ext(1:6, 1:6)  = A_ext11;
A_ext(7:12, 1:6) = A_ext21;
A_ext(7:12, 7:12) = A_ext22;
A_ext(1:6, 7:12)  = A_ext12;

