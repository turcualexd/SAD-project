w_E = 2*pi/(23*3600+56*60+47); % rad/sec Earth revolution velocity
R_E = 6371e3; % [m] Earth radius
a_E_S=1.495978707e8; % semimajor axis of Earth orbit to Sun
n_E=2*pi/(365*24*3600);
i_E=deg2rad(23.45);
e = 0;
th_0 = 0;
G=6.67e-11; 
mass_E=5.972e24; % [kg]
mu_E=G*mass_E;
R_E=6371e3; % [m]
a = R_E+500*1e3;
n_sc=sqrt(mu_E/(a^3)); % angular velocity of spacecraft orbit
am2 = 1; % we can change this value when we select the spacecraft and mission
Error_sun=0.5*pi/180; % 0.1 degrees of accuracy assumed; can be changed depending on sensor selection
Error_horizon=1.5*pi/180; % 1.5 degrees of accuracy assumed; can be changed depending on sensor selection
Error_mag_HYST=(0.6*60e-6); % taken from spacemaglite three axis magnometer document
Error_mag_ORTHO=deg2rad(0.5); % taken from spacemaglite three axis magnometer document 
T=1; % simulation step size time paramters will be determined later, this is a sample


% Coeficients for IGRF 2020

gnm = [-29404.8 -1450.9 0 0 0
       -2499.6 2982 1677 0 0
       1363.2 -2381.2 1236.2 525.7 0
       903.0 809.5 86.3 -309.4 48];
hnm = [0 4652.5 0 0 0
       0 -2991.6 -734.6 0 0
       0 -82.1 241.9 -543.4 0
       0 281.9 -158.4 199.7 -349.7];

% Earth's Magnetic Field

lenn = length(gnm(:,1));
lenm = length(gnm(1,:));
Snm = zeros(lenn,lenm);
Snm(1,1) = 1;

n = 0;
for j=2:lenn
    n = 1+n;
    Snm(j,1) = Snm(j-1,1)*(2*n-1)/n;
end 

k=1;
n=0;
for j=1:lenn
    n = n+1;
    for m=1:k
        if m==1
            kd = 1; %Kronecker delta
        else 
            kd = 0;
        end 
        Snm(j,m+1) = Snm(j,m)*( (kd+1)*(n-m+1)/(n+m))^(1/2);
        
    end 
    k = k+1;
end 

P = zeros(lenn,lenm);
dP = zeros(lenn,lenm);
gnm = Snm.*gnm;
hnm = Snm.*hnm;