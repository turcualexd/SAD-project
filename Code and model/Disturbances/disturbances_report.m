clear; close all; clc;

%------------------Space-Craft Geometry and Properties---------------------


I1_ext = 22.97240268194;
I2_ext = 7.77699984285;
I3_ext = 23.18950342186;

I_ext = diag([I1_ext I2_ext I3_ext]);
Iinv_ext = I_ext\eye(3);

%-----------------------------Orbit properties-----------------------------

muE   = astroConstants(13);
a    = 6.851221970532768e+03;
e    = 0.001830122216180;
i    = 1.699980862034725;
om   = 1.772848103192913;
OM   = 0.554268509489784;
%th_0 = 3.109851878856139;
th_0 = 0;
nn    = sqrt(muE/a^3);
T    = 2*pi/nn;

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
n_sun = sqrt(mu_sun/a_E^3);
%------------Initial Conditions % constant values for Simulation-----------

%--------------------Dynamics and Kinematics subsystems--------------------

omega0  = zeros(3,1);
s0      = zeros(3,1);

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
load("IGRF13_2020.mat", "gh");

%----------------------------------SRP-------------------------------------

P      = astroConstants(31)/(astroConstants(5)*10^3); % solar pressure
N      = [0 1 0; 0 -1 0; sqrt(2)/2 0 sqrt(2)/2; -sqrt(2)/2 0 sqrt(2)/2; -sqrt(2)/2 0 -sqrt(2)/2; sqrt(2)/2 0 -sqrt(2)/2; 0 0 1; 0 0 -1; 0 0 1; 0 0 -1; 0 0 1; 0 0 -1; 0 0 1; 0 0 -1]; %normal vectors
rho_s  = [0.5 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]; 
c1     = [1 - rho_s; 1 - rho_s; 1 - rho_s];
rho_d  = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
A      = [0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49 0.49]; %area of panels [m^2]
[n, ~] = size(N);
onesM  = ones(3, n);
onesV  = ones(n,1);
onesVV = ones(3*n,1);
r_F = [0 0.35 0; 0 -0.35 0; 0.35/sqrt(2) 0 0.35/sqrt(2); -0.35/sqrt(2) 0 0.35/sqrt(2);-0.35/sqrt(2) 0 -0.35/sqrt(2); 0.35/sqrt(2) 0 -0.35/sqrt(2); 1.8 0 0; 1.8 0 0; 0.35*sqrt(2) + 0.48 0 0; 0.35*sqrt(2) + 0.48 0 0; -1.8 0 0; -1.8 0 0; -(0.35*sqrt(2) + 0.48) 0 0; -(0.35*sqrt(2) + 0.48) 0 0];
[n,~] = size(r_F);
r_Fskew = zeros(n*3, n*3);
j = 1;

for ii = 1:3:3*n-2
    r_Fskew(ii:ii+2,ii:ii+2) = skew(r_F(j,:));
    j = j+1;
end



%---------------------------------Air drag------------------------------

CD=2.2;              % drag coefficient, usually is between 1.5 and 2.6
d0=1.585*10^-12;     % stardart value of the density for the height h0 (values can be found in professor notes)
h0=450;              % standard height for the density equation model between 500 and 600 km of altitude from Earth
H=60.828; 


%---------------------------SIMULATION-------------------------------

t0 = 0;
tf = 3*T;
sim_options.SolverType='Fixed-Step';
sim_options.Solver='ode5';
sim_options.FixedStep='0.3';
sim_options.StartTime='t0';
sim_options.StopTime='tf';
out   = sim("disturbances_report_sim.slx", sim_options);
srpT  = out.SRP_torque.data;
magT  = squeeze(out.MAG_torque);
dragT = out.DRAG_torque.data;
ggT   = out.GRAV_torque.data;
t     = out.tout;


%%
figure;
subplot(2,2,1)
plot(t/T, magT, 'LineWidth', 2);
xlabel('Orbits','interpreter','latex'); ylabel('Magnetic Torque [Nm]','interpreter','latex')
title('Magnetic Torque','interpreter','latex');
lgd1 = legend('Mx', 'My', 'Mz','interpreter','latex','fontsize', 3);
grid on;

subplot(2,2,2)
plot(t/T, dragT, 'LineWidth',2);
xlabel('Orbits','interpreter','latex'); ylabel('Drag Torque [Nm]','interpreter','latex')
title('Drag Torque','interpreter','latex');
grid on;

subplot(2,2,3)
plot(t/T, ggT, 'LineWidth',2)
xlabel('Orbits','interpreter','latex'); ylabel('GG Torque [Nm]','interpreter','latex')
title('GG Torque','interpreter','latex');
grid on;

subplot(2,2,4)
plot(t/T, srpT, 'LineWidth',2)
xlabel('Orbits','interpreter','latex'); ylabel('SRP Torque [Nm]','interpreter','latex')
title('SRP Torque','interpreter','latex');
grid on;
