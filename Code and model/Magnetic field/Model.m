clear, clc, close all;

%% Import IGRF 2020 coefficients from file

if ~isfile("IGRF13coeffs.mat")
     import_IGRF2020("IGRF13coeffs.xls")
end

load("IGRF13coeffs.mat")

%% Summatory to obtain V

R = 6371.2;     % Earth radius [km]
r = R + 100;
phi = 30;
theta = 30;
n_max = 13;      % order of the magnetic field

if (theta>-0.00000001 && theta<0.00000001)
    theta=0.00000001;
elseif(theta<180.00000001 && theta>179.99999999)
    theta=179.99999999;
end

theta= (90 - theta) * pi/180;
phi = phi * pi/180;

costheta = cos(theta);
sintheta = sin(theta);

m_ind = @(i) i+1;

K = @(n,m) ((n-1)^2 - m^2) / ((2*n - 1) * (2*n - 3));

P = zeros(m_ind(n_max), m_ind(n_max));
P(m_ind(0), m_ind(0)) = 1;
P(m_ind(1), m_ind(0)) = costheta;
P(m_ind(1), m_ind(1)) = sintheta;

dP = P;
dP(m_ind(0), m_ind(0)) = 0;
dP(m_ind(1), m_ind(0)) = -sintheta;
dP(m_ind(1), m_ind(1)) = costheta;

for n = 2:n_max
    for m = 0:n
        if n == m
            P(m_ind(n), m_ind(m)) = sintheta * P(m_ind(n-1), m_ind(m-1));
            dP(m_ind(n), m_ind(m)) = sintheta * dP(m_ind(n-1), m_ind(m-1)) + costheta * P(m_ind(n-1), m_ind(m-1));
        else
            P(m_ind(n), m_ind(m)) = costheta * P(m_ind(n-1), m_ind(m)) - K(n,m) * P(m_ind(n-2), m_ind(m));
            dP(m_ind(n), m_ind(m)) = costheta * dP(m_ind(n-1), m_ind(m)) - sintheta * P(m_ind(n-1), m_ind(m)) - K(n,m) * dP(m_ind(n-2), m_ind(m));
        end
    end
end

P(1,:) = [];
dP(1,:) = [];

Br = 0; Bt = 0; Bp = 0;
for n = 1:n_max
    for m = 0:n
        Br = Br + (R/r)^(n+2) * (n+1) * (g(n,m_ind(m)) * cos(m*phi) + h(n,m_ind(m)) * sin(m*phi)) * P(n, m_ind(m));
    end
end

Br