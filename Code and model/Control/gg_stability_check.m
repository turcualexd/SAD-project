clear; close all; clc;

I_x = 21.13;
I_y = 7.864;
I_z = 21.345;
I_vect = [I_x; I_y; I_z];

k_yaw  = (I_z - I_y) / I_x;
k_roll = (I_z - I_x) / I_y;
k_pitch= (I_y - I_x) / I_z;

figure;
plot(linspace(-1,1,100), ones(100,1),'Color', 'b'); 
hold on
plot(linspace(-1,1,100), -ones(100,1),'Color', 'b');
plot(ones(100,1), linspace(-1,1,100), 'Color', 'b');
plot(-ones(100,1),linspace(-1,1,100),'Color', 'b');

plot(linspace(-1,1,100),zeros(100,1),  'Color', 'k');
plot(zeros(100,1),linspace(-1,1,100),'Color', 'k');
plot(linspace(-1,1,100),linspace(-1,1,100), 'Color', [.7 .7 .7], 'LineWidth', 0.5)
axis equal;
plot(k_roll, k_pitch, 'Marker', 'o', 'Color', 'r');