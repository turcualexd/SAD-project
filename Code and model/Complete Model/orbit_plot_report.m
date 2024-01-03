clear; close all; clc;

%-----------------------------Orbit properties-----------------------------

mu   = astroConstants(13);
a    = 6.851221970532768e+03;
e    = 0.001830122216180;
i    = 1.699980862034725;
om   = 1.772848103192913;
%OM   = 0.554268509489784;
OM = 0;
th_0 = 3.109851878856139;
n    = sqrt(mu/a^3);
T    = 2*pi/n;

%--------------------------Position of Earth wrt Sun-----------------------

% Sun Centred equatorial frame 
mu_sun = astroConstants(4);
e_E    = 1.743742542622805E-02;
i_E    = deg2rad(2.343620856296279E+01);
OM_E   = deg2rad(3.599997935134700E+02);
om_E   = deg2rad(1.030408771388982E+02);
th_E   = deg2rad(3.403771999857636E+02);
a_E    = 1.497081833259778E+08;

[rr, vv] = kep2car(a, e, i, OM, om, th_0, mu);
rs = -rr;
mjd2000date = date2mjd2000([2023, 12, 16, 00, 00, 00]);
[kep, ksun] = uplanet(mjd2000date,3);


%ECI
Terra_3D
[~,r,v]=plotOrbit(rr, vv,1e6, "color", 	'#D95319');
hold on;
quiver3(0,0,0,1.4*1e4,0,0,'Color','#0072BD','LineWidth',2);
quiver3(0,0,0,0,1.4*1e4,0,'Color','#0072BD','LineWidth',2);
quiver3(0,0,0,0,0,1.4*1e4,'Color','#0072BD','LineWidth',2);

[r_Earth, v_Earth] = kep2car(a_E, e_E, i_E, OM_E, om_E, th_E, mu_sun);

r_toSUN = -r_Earth./(1*1e4);
quiver3(0,0,0,r_toSUN(1),r_toSUN(2),r_toSUN(3),'Color',"#EDB120",'LineWidth',2)

legend('','Orbit','ECIEq frame', '','','Sun direction','fontsize',2,'interpreter','latex')


%Ecliptic
Terra_3D;
eps = deg2rad(23.44);
xp = [1 0 0];
yp = [0 cos(eps) -sin(eps)];
zp = [0 sin(eps) cos(eps)];
ROT = [xp; yp; zp]'; %from equatorial to ecliptic
for k = 1:length(r)
    r_ecl(:,k) = ROT*r(k,:)';
end

plot3(r_ecl(1,:), r_ecl(2,:),r_ecl(3,:))


[r_Earth, v_Earth] = kep2car(a_E, e_E, i_E, OM_E, om_E, th_E, mu_sun);

r_toSUN = -r_Earth./1e4;






