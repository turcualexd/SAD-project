clear; close all; clc;

I_x = 21.13;
I_y = 7.864;
I_z = 21.345;
I_vect = [I_x; I_y; I_z];

k_yaw  = (I_z - I_y) / I_x;
k_roll = (I_z - I_x) / I_y;
k_pitch= (I_y - I_x) / I_z;

figure;
plot(linspace(-1,1,100), ones(100,1),'Color', 'b', 'LineWidth', 1.5); 
hold on
plot(linspace(-1,1,100), -ones(100,1),'Color', 'b', 'LineWidth', 1.5);
plot(ones(100,1), linspace(-1,1,100), 'Color', 'b', 'LineWidth', 1.5);
plot(-ones(100,1),linspace(-1,1,100),'Color', 'b', 'LineWidth', 1.5);

plot(linspace(-1,1,100),zeros(100,1),  'Color', 'k', 'LineWidth', 1.5);
plot(zeros(100,1),linspace(-1,1,100),'Color', 'k', 'LineWidth', 1.5);
plot(linspace(-1,1,100),linspace(-1,1,100), 'Color', 'r', 'LineWidth', 1.5)
axis equal;
plot(k_roll, k_pitch, 'Marker', 'x', 'Color', 'r', 'LineWidth', 1.5);

fimplicit(@(x,y) (1 + 3.*x + x.*y) - 4*sqrt(x.*y), [-1 0 -1 0],'Color','r', 'LineWidth', 1.5);
xlabel('$K_r$', 'interpreter', 'latex'); ylabel('$K_y$', 'interpreter', 'latex');


impl = @(x,y) (1 + 3.*x + x.*y) - 4*sqrt(x.*y);

k_yaw
k_roll
k_pitch