clear; close all; clc; 
I1 = 1.009; I2 = 0.251; I3 = 0.916;
I = diag([I1 I2 I3]);
Iinv = I\eye(3);
Aln0 = eye(3);
Abn0 = Aln0;
r0 = [26578.137; 0; 0];
v0 = [0; 2.221; 3.173];
mu = astroConstants(13);%Km3/s2
[a, e, i, OM, om, theta0] = car2par(r0, v0, mu);
i_matrix = [1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];
n = sqrt(mu/a^3);
T = 2*pi/n;
t0 = 0;
tf = T;
omega0 = [0; 0; n];
e3 = [1; 0; 0];
rs = astroConstants(2); %Astronomical Unit (AU) Km
epsilon = astroConstants(8); %obliquity rad
Fe = astroConstants(31); %intensity of radiation at 1 AU W/m2
c = astroConstants(5); %speed of light in vacuum (Km/s)
c = c * 1e3; %m/s
P = Fe/c; %radiation pressure Pa
Ts = 365.2564*24*60*60; %sidereal year in seconds
ns = 2*pi / Ts;

N = [1 0 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 1; 0 0 -1; 1 0 0; -1 0 0; 1 0 0; -1 0 0];
A = [6*1e-2; 6*1e-2; 6*1e-2; 6*1e-2; 4*1e-2; 4*1e-2; 12*1e-2; 12*1e-2; 12*1e-2; 12*1e-2];
rhos = 0.5*ones(10,1);
rhos(7:10)=0.1;
rhos_=ones(10,1)-rhos;
rhod = 0.1*ones(10,1);
rf = 1e-2*[10 0 0; 0 10 0; -10 0 0; 0 -10 0; 0 0 15; 0 0 -15; 0 0 45; 0 0 45; 0 0 -45; 0 0 -45];
s0 = zeros(3,1);
A0 = eye(3);
g01 = -29404.8; g11 = -1450.9; h11 = 4652.5;
H0 = sqrt(g01^2+g11^2+h11^2)*1e-9;
incl = deg2rad(11.5);
omegaE = deg2rad(15)/3600;%angular velocity of earth in rad/s
Rt = astroConstants(23);
jb = [0.01; 0.05; 0.01];
%%
sim_options.SolverType='Fixed-Step';
sim_options.Solver='ode5';
sim_options.FixedStep='1';
sim_options.StartTime='t0';
sim_options.StopTime='tf';
out = sim("sat_dist_sim.slx", sim_options);
%%
Abn = out.Abn;
t = out.tout;
s = out.angles;
mag_tor = out.mag_tor;
Mg = out.Mg;
%%
omegab = NaN(3, length(t));
h_b = NaN(3, length(t));
h_n = NaN(3, length(t));
SRP = NaN(3,length(t));
dist = NaN(3,length(t));
for i = 1 : length(t)
    omegab(:,i) = out.omegab(:,1,i);
    h_b(:,i) = I * omegab(:,i);
    h_n(:,i) = Abn(:,:,i)'*h_b(:,i);
    SRP(:,i) = out.SRP(:,1,i);
    dist(:,i) = out.dist(:,1,i);
end
%%
figure
hold on
grid minor
plot(t/T,omegab(1,:))
plot(t/T,omegab(2,:))
plot(t/T,omegab(3,:))
lgd = legend('\omega_1', '\omega_2', '\omega_3');
lgd.FontSize=20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('omega_body [rad/s]');
yy.FontSize = 20;
%%
h_bn = vecnorm(h_b, 2, 1);
h_nn = vecnorm(h_n, 2, 1);
SRPn = vecnorm(SRP, 2, 1);
%%
figure
hold on
grid minor
plot(t/T, h_n(1,:))
plot(t/T, h_n(2,:))
plot(t/T, h_n(3,:))
plot(t/T, h_nn)
lgd = legend('h_1', 'h_2', 'h_3', 'norm');
lgd.FontSize = 20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('h_n  [kg*m^2/s]');
yy.FontSize = 20;
%%
figure
hold on
grid minor
plot(t/T,SRP(1,:))
plot(t/T,SRP(2,:))
plot(t/T,SRP(3,:))
plot(t/T, SRPn)
lgd = legend('SRP_1', 'SRP_2', 'SRP_3', '||SRP||');
lgd.FontSize = 20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('SRP torque [N*m]');
yy.FontSize = 20;
%%
figure
hold on
grid minor
plot(t/T, h_b(1,:))
plot(t/T, h_b(2,:))
plot(t/T, h_b(3,:))
plot(t/T, h_bn)
lgd = legend('h_1', 'h_2', 'h_3', 'norm');
lgd.FontSize = 20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('h_b  [kg*m^2/s]');
yy.FontSize = 20;
%%
figure
hold on
grid minor
plot(t/T, mag_tor(:,1))
plot(t/T, mag_tor(:,2))
plot(t/T, mag_tor(:,3))
lgd = legend('mag_1', 'mag_2', 'mag_3');
lgd.FontSize = 20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('magnetic torque [N*m]');
yy.FontSize = 20;
%%
figure
hold on
grid minor
plot(t/T, Mg(:,1))
plot(t/T, Mg(:,2))
plot(t/T, Mg(:,3))
lgd = legend('Mg_1', 'Mg_2', 'Mg_3');
lgd.FontSize = 20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('gravity gradient torque [N*m]');
yy.FontSize = 20;
%%
figure
hold on
grid minor
plot(t/T,dist(1,:))
plot(t/T,dist(2,:))
plot(t/T,dist(3,:))
lgd = legend('dist_1', 'dist_2', 'dist_3');
lgd.FontSize = 20;
xx = xlabel('t/T [-]');
xx.FontSize = 20;
yy = ylabel('torque [N*m]');
yy.FontSize = 20;