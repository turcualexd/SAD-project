clear, clc, close all;

a = 6371.2;

h = 200;
r = a + h;

theta = 30;
phi = 30;

B = calculateField(r, theta, phi)