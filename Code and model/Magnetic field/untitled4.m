%clear, clc, close all;

a=6371.2; % Reference radius used in IGRF
r = a+100;
theta = 30;
phi = 30;
N = 13;

[X1,Y1,Z1] = magnet(r,theta,phi)
[X2, Y2, Z2] = igrfs('01/01/2020', theta, phi, r, 'geocentric')