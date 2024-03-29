function B = calculateField(r, theta, phi)

N = 13;         % Order of the field
lat = pi/2 - theta;
lon = phi;

load("IGRF13_2020.mat", "gh");

R = 6371.2;     % Radius of Earth in IGRF
ratio = R / r;

sin_lat = sin(lat);
cos_lat = cos(lat);
sin_lon = sin(lon);
cos_lon = cos(lon);

npq = (N * (N+3)) / 2;

p = zeros(npq, 1);
p(1) = 2*sin_lat;
p(2) = 2*cos_lat;
p(3) = 4.5*sin_lat^2 - 1.5;
p(4) = 3*sqrt(3) * cos_lat * sin_lat;

q = zeros(npq, 1);
q(1) = -cos_lat;
q(2) = sin_lat;
q(3) = -3 * cos_lat * sin_lat;
q(4) = sqrt(3) * (sin_lat^2 - cos_lat^2);

sl = zeros(N+1, 1);
sl(1) = sin_lon;

cl = zeros(N+1, 1);
cl(1) = cos_lon;

Bx = 0; By = 0; Bz = 0;
l = 1; n = 0; m = 1;

for k = 1:npq

    if n < m
        m = 0;
        n = n + 1;
        rr = ratio^(n+2);
        fn = n;
    end

    fm = m;

    if k >= 5

        if m == n
            a = sqrt(1 - 0.5/fm);
            j = k - n - 1;
            p(k) = (1 + 1/fm) * a * cos_lat * p(j);
            q(k) = a * (cos_lat * q(j) + sin_lat/fm * p(j));
            sl(m) = sl(m-1)*cl(1) + cl(m-1)*sl(1);
            cl(m) = cl(m-1)*cl(1) - sl(m-1)*sl(1);
        else
            a = sqrt(fn^2 - fm^2);
            b = sqrt((fn-1)^2 - fm^2) / a;
            c = (2*fn - 1) / a;
            ii = k - n;
            j = k - 2*n + 1;
            p(k) = (fn+1) * (c * sin_lat/fn * p(ii) - b/(fn-1) * p(j));
            q(k) = c * (sin_lat * q(ii) - cos_lat/fn * p(ii)) - b * q(j);
        end

    end

    a = rr * gh(l);
    
    if m == 0

        Bx = Bx + a*q(k);
        Bz = Bz - a*p(k);
        
        l = l + 1;

    else
        
        b = rr * gh(l+1);
        c = a*cl(m) + b*sl(m);
        Bx = Bx + c * q(k);
        Bz = Bz - c * p(k);
        
        if cos_lat > 0
            By = By + (a*sl(m) - b*cl(m)) * fm * p(k)/((fn+1) * cos_lat);
        else
            By = By + (a*sl(m) - b*cl(m)) * q(k) * sin_lat;
        end
        
        l = l + 2;

    end

    m = m + 1;

end

B_NED = [Bx; By; Bz];

% Convert NED to ECEF
R = [-sin_lat * cos_lon,    -sin_lon,   -cos_lat * cos_lon;
     -sin_lat * sin_lon,    cos_lon,    -cos_lat * sin_lon;
     cos_lat,               0,          -sin_lat            ];

B = R * B_NED;