function [a, e, i, OM, om, th] = car2par(rr, vv, unit, mu)

% Trasformation from cartesian coordinates to Keplerian parameters
% 
% [a, e, i, OM, om, th] = car2par(rr, vv, unit, mu)
% 
% -----------------------------------------------------------------------------------
% Input arguments:
% rr            [1x1]       position vector                                 [km]
% vv            [1x1]       velocity vector                                 [km/s]
% unit          [string]    unit of measure of angles ("rad" or "deg")      [-]
% mu            [1x1]       gravitational parameter                         [km^3/s^2]
% 
% -----------------------------------------------------------------------------------
% Output arguments:
% a             [1x1]       semi-major axis                                 [km]
% e             [1x1]       eccentricity                                    [-]
% i             [1x1]       inclination                                     [rad/deg]
% OM            [1x1]       RAAN (Right Ascension of the Ascending Node)    [rad/deg]
% om            [1x1]       pericenter anomaly                              [rad/deg]
% th            [1x1]       true anomaly                                    [rad/deg]
% 
% -----------------------------------------------------------------------------------

% If mu is not assigned, the default value is set to Earth
switch nargin
    case 2
        strcmp(unit,'deg');
        mu = 3.986 * 10^5;
    case 3
        mu = 3.986 * 10^5;
end

% -------------------------------------------------------------

% Finds a
r = norm(rr);
v = norm(vv);
a = 1 / ( (2/r) - (v^2 / mu) );

% Finds e
hh = cross(rr, vv);
h = norm(hh);
ee = cross(vv, hh)./mu - rr./r;
e = norm(ee);

% Finds i
i = acos(hh(3) / h);

% Calculates N vector (line of nodes)
if i >= 1e-5
    k = [0 0 1]';
    N = cross(k, hh) / norm( cross(k, hh) );
else
    N = [1 0 0]';
end


% Finds OM
OM = acos( N(1) );
if N(2) < 0
    OM = 2*pi - OM;
end

% Finds om
om = acos( dot(N, ee) / e );
if ee(3) < 0
    om = 2*pi - om;
end

% Finds th
vr = dot(vv, rr) / r;
th = acos( dot(rr, ee) / (r*e) );
if vr < 0 
    th = 2*pi - th;
end

% -------------------------------------------------------------

% If output angles have to be in degrees, does the conversion
if strcmp(unit,'deg')
    i = rad2deg(i);
    OM = rad2deg(OM);
    om = rad2deg(om);
    th = rad2deg(th);
end