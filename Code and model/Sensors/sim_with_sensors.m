clear; close all; clc;

I1 = 7776999842.85;

I2 = 22972402681.94;

I3 = 23189503421.86;

I = diag([I1 I2 I3])*1e-9;

Iinv = I\eye(3);

Aln0 = eye(3);

Abn0 = Aln0;

mu = astroConstants(13);%Km3/s2

e = 0.001830122216180;

a = 6.851221970532768e+03;

om = 1.772848103192913;

i = 1.699980862034725;

OM = 0.554268509489784;

theta0 = 3.109851878856139;

om_matrix = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];

OM_matrix = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];

i_matrix = [1 0 0;0 cos(i) sin(i); 0 -sin(i) cos(i)];

R_mat = om_matrix * i_matrix * OM_matrix;

n = sqrt(mu/a^3);

T = 2*pi/n;

t0 = 0;

tf = T;

omega0 = [0.01; 0.01; 0.01]*6;

s0 = zeros(3,1);

A0 = eye(3);

g01 = -29404.8; g11 = -1450.9; h11 = 4652.5;

H0 = sqrt(g01^2+g11^2+h11^2)*1e-9;

incl = deg2rad(11.5);

omegaE = deg2rad(15)/3600;%angular velocity of earth in rad/s

Rt = astroConstants(23);

jb = [0.01; 0.05; 0.01];

wE = [0; 0; omegaE];

incl_matrix = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];

wE = incl_matrix * wE;

% for B - dot control define gain k (reference: Markley Crassidirs)


chsi = (i - deg2rad(11)); % relative inclination of orbit wrt to equatorial magnetic plane

k_gain  = (4*pi/T) * (1 + sin(chsi)) * I1 *1e-9;
