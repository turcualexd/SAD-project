clear; close all; clc;

%------------------Space-Craft Geometry and Properties---------------------

I1_retr = 49/6;
I2_retr = 49/6;
I3_retr = 49/6;

I_retr = diag([I1_retr I2_retr I3_retr]);
Iinv_retr = I_retr\eye(3);

I1_ext = 22.97240268194;
I2_ext = 7.77699984285;
I3_ext = 23.18950342186;

I_ext = diag([I1_ext I2_ext I3_ext]);
Iinv_ext = I_ext\eye(3);

%-----------------------------Orbit properties-----------------------------

mu   = astroConstants(13);
a    = 6.851221970532768e+03;
e    = 0.001830122216180;
i    = 1.699980862034725;
om   = 1.772848103192913;
OM = 0;
th_0 = 0;
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

omega0  = -0.06 + 0.12 * rand(3,1);
s0      = 2*pi*rand(3,1);

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

%--------------------------------Sensors-----------------------------------

% Magnetometer - AAC ClydeSpace MM200

freq_mag     = 30;         % Can be modulated up to 500Hz
noise_sd_mag = 1.18*1e-9;  % T/sqrt(Hz)
pwr_sd_mag   = noise_sd_mag^2;
bias_mag     = rand(3,1);
bias_mag     = bias_mag./norm(bias_mag);

% Sun Sensor   - Tensor Tech FSS-15 

freq_sun     = 30; 
bias_sun     = rand(3,1);
bias_sun     = bias_sun./norm(bias_sun);
variance_sun = deg2rad(0.3)^2;
pwr_sd_sun   = (1/freq_sun)*variance_sun;

% Earth Sensor - Meisei Electric
freq_earth     = 30; 
FOV_earth      = deg2rad(33);   
bias_earth     = rand(3,1);
bias_earth     = bias_earth./norm(bias_earth);
variance_earth = deg2rad(1)^2;
pwr_sd_earth   = (1/freq_earth)*variance_earth;

%--------------------------------------------------------------------------

freq_acds = 30;

%---------------------------Att. Determination-----------------------------

alpha1 = 0.05;
alpha2 = 0.15;
alpha3 = 0.8;

alpha12 = 0.2; 
alpha22 = 0.8;

%-------------------------------Actuators----------------------------------

max_dip   =  15; % [A m^2]
freq_act  = 30;
max_hdot = 0.1;
w_sat = (2 * pi * 6000 / 60);
I_wheel = 4 / w_sat;

%-----------------------------De-Tumbling----------------------------------

chsi = (i - deg2rad(11));
k_gain  = (4*pi/T) * (1 + sin(chsi)) * I1_retr;

out = sim("sat_model_delivery.slx");

%%
clc
close all;
t = out.tout;
omega = out.omega;
omega = squeeze(omega);
Mc = out.Mc;
alpha = out.alpha;
alpha = squeeze(alpha);
Abl = out.Abl;
hat_Bm_dot = out.hat_Bm_dot;
hat_Bm_dot = squeeze(hat_Bm_dot);
hat_Bm_dot_norm = vecnorm(hat_Bm_dot);

alpha_dot = nan(size(alpha,1), size(alpha,2));
for i = 1:size(alpha,1)
    alpha_dot(i,:) = omega(:,i) - Abl(:,:,i) * [0 0 n].';
end

alpha = rad2deg(alpha);
alpha_dot = rad2deg(alpha_dot);

% omega - true angular velocity
figure
hold on
grid on
plot(t/T, omega, 'LineWidth', 1.5)
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$\omega$ [rad/s]', 'Interpreter', 'latex')
legend('$\omega_x$', '$\omega_y$', '$\omega_z$', 'Interpreter', 'latex')

% Bm_dot - derivative of direction of measured magnetic field
figure
hold on
grid on
plot(t/T, hat_Bm_dot, 'LineWidth', 1.5)
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$\dot{B}_m$ [1/s]', 'Interpreter', 'latex')
legend('$\dot{B}_{m,x}$', '$\dot{B}_{m,y}$', '$\dot{B}_{m,z}$', 'Interpreter', 'latex')

% Bm_dot_norm - norm of derivative of magnetic field direction
figure
hold on
grid on
plot(t/T, hat_Bm_dot_norm, 'LineWidth', 1.5)
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$||\dot{B}_m||$ [1/s]', 'Interpreter', 'latex')

% alpha - attitude error angles
figure
hold on
grid on
plot(t/T, alpha, 'LineWidth', 1.5)
xlim([0.5 2])
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$\alpha [deg]$', 'Interpreter', 'latex')
legend('$\alpha_x$', '$\alpha_y$', '$\alpha_z$', 'Interpreter', 'latex')

% alpha - zoom of previous plot
figure
hold on
grid on
plot(t/T, alpha, 'LineWidth', 1.5)
xlim([1 2])
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$\alpha [deg]$', 'Interpreter', 'latex')
legend('$\alpha_x$', '$\alpha_y$', '$\alpha_z$', 'Interpreter', 'latex')

% alpha_dot - attitude error angle rate (omega_BL)
figure
hold on
grid on
plot(t/T, alpha_dot, 'LineWidth', 1.5)
xlim([0.5 2])
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$\dot{\alpha} [deg/s]$', 'Interpreter', 'latex')
legend('$\dot{\alpha}_x$', '$\dot{\alpha}_y$', '$\dot{\alpha}_z$', 'Interpreter', 'latex')

% alpha_dot - attitude error angle rate (omega_BL)
figure
hold on
grid on
plot(t/T, alpha_dot, 'LineWidth', 1.5)
xlim([1 2])
xlabel('Periods of orbit', 'Interpreter', 'latex')
ylabel('$\dot{\alpha} [deg/s]$', 'Interpreter', 'latex')
legend('$\dot{\alpha}_x$', '$\dot{\alpha}_y$', '$\dot{\alpha}_z$', 'Interpreter', 'latex')