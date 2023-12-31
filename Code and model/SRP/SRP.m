clear; close all; clc;

I = diag([100.09 25.1 91.6]);
I_inv = inv(I);
w0 = [0.2; 0; 0.2];
s0 = [0; 0; 0];


OM = 0;
om = deg2rad(90);
i = deg2rad(98);

R_OM =  [ cos(OM),    sin(OM),    0;
          -sin(OM),   cos(OM),    0;
          0,          0,          1 ];

R_i =   [ 1,          0,          0;
          0,          cos(i),     sin(i);
          0,          -sin(i),    cos(i) ];

R_om =  [ cos(om),    sin(om),    0;
          -sin(om),   cos(om),    0;
          0,          0,          1 ];

A_313 = ( R_om * R_i * R_OM )'; %from perifocal to inertial

% Earth Sensor

optical = [1 0 0]';
A_MIS   = eye(3);
FOV = deg2rad(30);


%


t0 = 0;
tf = 100;
step = 0.1;
tol = 1e-1;

P = astroConstants(31)/(astroConstants(5)*10^3);
N = [1 0 0; 0 1 0; -1 0 0; 0 -1 0; 0 0 1; 0 0 -1; 1 0 0; -1 0 0; 1 0 0; -1 0 0];
rho_s = [0.5 0.5 0.5 0.5 0.5 0.5 0.1 0.1 0.1 0.1];
c1    = [1 - rho_s; 1 - rho_s; 1 - rho_s];
rho_d  = [0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1];
A = 1e-2*[6 6 6 6 4 4 12 12 12 12];
[n, ~] = size(N);
onesM = ones(3, n);
onesV = ones(n,1);
tol = 1e-3;
n_sun = 2*pi/(365*24*60*60);
eps = deg2rad(23.45);
e = 0.007;
a = 7500;
i = deg2rad(32.5934);


s0 = [0 0 0]';
 

r_F = 1e-2*[10 0 0; 0 10 0; -10 0 0; 0 -10 0; 0 0 15; 0 0 -15; 0 0 45; 0 0 45; 0 0 -45; 0 0 -45]+ 1e-2*[0 0 3];
[n,~] = size(r_F);
r_Fskew = zeros(n*3, n*3);
j = 1;

for i = 1:3:3*n-2
    r_Fskew(i:i+2,i:i+2) = skew(r_F(j,:));
    j = j+1;
end

out = sim("SRP.slx", "StartTime", "t0", "StopTime", "tf", "MaxStep", "step");

t = out.tout;
s = out.s.Data;
ss = s(end,:);
F = zeros(3, n);
T = zeros(3,n);
sss = zeros(3,1);
sumT = zeros(3,length(t));
AA = out.A.Data;
for j = 1:length(t)
    sss = AA(:,:,j)*[cos(n_sun*t(j)); sin(n_sun*t(j))*cos(eps); sin(n_sun*t(j))*sin(eps)];
    sss = sss';
    for i = 1:n 
        if dot(sss,N(i,:)) < 0
         F(:,i) = -P*A(i)*dot(sss,N(i,:)) * ((1 - rho_s(i))*sss + (((2*rho_s(i))*dot(sss,N(i,:))) + (2/3)*rho_d(i))*N(i,:));
         
        else
         F(:,i) = [0; 0; 0];
        end
        
        T(:,i) = cross(r_F(i,:), F(:,i));   
        sumT(:,j) = sumT(:,j) + T(:,i);
            
    end
    
 end
G = astroConstants(1);
Mt = astroConstants(13)/G;
R = astroConstants(23);

SRP_F = out.SRP_F.data;
SRP_T = out.SRP_T.data;


figure;
plot(t, SRP_T,'b')

figure;
plot(t, sumT, 'r')
