function B = calculateField(r, lat, lon, N)

% Ricorda: theta = 90 - lat, phi = lon

load("IGRF13_2020.mat", "gh");

R = 6371.2;     % Radius of Earth in IGRF
ratio = R / r;

sin_lat = sind(lat);
cos_lat = cosd(lat);

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
sl(1) = sind(lon);

cl = zeros(N+1, 1);
cl(1) = cosd(lon);

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
            aa = sqrt(1 - 0.5/fm);
            j = k - n - 1;
            p(k) = (1 + 1/fm) * aa * cos_lat * p(j);
            q(k) = aa * (cos_lat * q(j) + sin_lat/fm * p(j));
            sl(m) = sl(m-1)*cl(1) + cl(m-1)*sl(1);
            cl(m) = cl(m-1)*cl(1) - sl(m-1)*sl(1);
        else
            aa = sqrt(fn^2 - fm^2);
            bb = sqrt((fn-1)^2 - fm^2) / aa;
            cc = (2*fn - 1) / aa;
            ii = k - n;
            j = k - 2*n + 1;
            p(k) = (fn+1) * (cc * sin_lat/fn * p(ii) - bb/(fn-1) * p(j));
            q(k) = cc * (sin_lat * q(ii) - cos_lat/fn * p(ii)) - bb * q(j);
        end

    end

    aa = rr * gh(l);
    
    if m == 0

        Bx = Bx + aa*q(k);
        Bz = Bz - aa*p(k);
        
        l = l + 1;

    else
        
        bb = rr * gh(l+1);
        cc = aa*cl(m) + bb*sl(m);
        Bx = Bx + cc * q(k);
        Bz = Bz - cc * p(k);
        
        if cos_lat > 0
            By = By + (aa*sl(m) - bb*cl(m)) * fm * p(k)/((fn+1) * cos_lat);
        else
            By = By + (aa*sl(m) - bb*cl(m)) * q(k) * sin_lat;
        end
        
        l = l + 2;

    end

    m = m + 1;

end

B = [Bx; By; Bz];