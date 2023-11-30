clear, clc, close all;

%% Initial conditions

I_vect = [0.07; 0.055; 0.025];
w0 = [0.5; -0.2; 1.0];
s0 = [0; asin(0.4); pi/3];

I = diag(I_vect);
I_inv = inv(I);

%% Simulation

t0 = 0;
tf = 200;
step = 0.1;
tol = 0.2;

out = sim("Kinematics_switch_sim.slx", "Solver", "ode5", "StartTime", "t0", "StopTime", "tf", "SolverType", "Fixed-Step", "FixedStep", "step");

%% Output

w = out.w.data;
t = out.tout;
b = out.b.data;
A = out.A.data;
h = out.h.data;

figure
hold on
plot(t, h(:,1))
plot(t, h(:,2))
plot(t, h(:,3))
plot(t, b*0.01, 'Marker', 'o')