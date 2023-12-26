clear, clc, close all;

N = 13;
a=6371.2;

h = 200;
r = a + h;

lat = -30;
lon = 30;

% XYZ = igrfmagm(h,lat,lon,decyear(2020,1,1),N).'
B = calculateField(r, lat, lon, N)