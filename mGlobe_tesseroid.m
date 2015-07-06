function dg = mGlobe_tesseroid(point_of_calc,lon0,lat0,altitude0,height,density,delta_lonT,delta_latT,radius,sum_switch)
%MGLOBE_TESSEROID Function calculates the gravity effect of tesseroid
% Calculate the gravitational effect of a tesseroid of given centre (centre 
% on sphere, not in height). 
% Calculation is based on Heck and Seitz paper "A
% comparison of tesseroid, prism and point-mass approaches from mass
% reduction in gravity field modelling", 2007, J Geodesy, 81:121-136,
% (doi:10.1007/s00190-006-0094-0).
% It is a spherical approximation (Taylor expansion) not suitable for near 
% zone (see paper)! Use higher resolution for local zone. Results were 
% compared to spherical cup (analytical solution).
%
% Input:
%   point_of_calc   ...     computation point [longitude (deg), 
%                           latitude (deg), height (m)] (vector (3x1))
%   lon0            ...     matrix with longitude of tesseroids centre
%                           projected on reference sphere (deg)
%   lat0            ...     matrix with latitude of tesseroids centre
%                           projected on reference sphere (deg)
%   altitude0       ...     matrix with altitude of tesseroids bottom 
%                           sphere (m)
%   height          ...     matrix with tesseroids height used for the 
%                           calculation of upper sphere (m)
%   density         ...     matrix with tesseroids density (kg/m^3)
%   delta_lon       ...     scalar or matrix (same dim. as lon0,lat0)
%                           with tesseroids dimension in longitudal
%                           direction
%   delta_lat       ...     scalar or matrix (same dim. as lon0,lat0)
%                           with tesseroids dimension in latitude
%                           direction
%                           [delta_longitude (m), delta_latitude (m)]
%   radius          ...     radius (scalar) of reference sphere (m)
%   sum_switch      ...     = 1 for scalar, = 0 for matrix
% 
% Output:
%   dg              ...     resulting gravity effect (scalar) in m/s^2
% 
% Example:
%   dg = mGlobe_tesseroid([14,45,500],[10 11;10 11],...
%        [40 41;42 43],[1000 1100;900,10],[10 10;10 10],...
%        [2700 2700;2700 2700],5/60,5/60,6371000,1)
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

% Computation point
lon = point_of_calc(1)*pi/180;
lat = point_of_calc(2)*pi/180;
sinlat = sin(lat);
coslat = cos(lat);
r = radius + point_of_calc(3);
r2 = r^2;
% Running integration point
lon0 = lon0*pi/180;
lat0 = lat0*pi/180;
coslat0 = cos(lat0);
sinlat0 = sin(lat0);
sinlat02 = sinlat0.^2;
r0 = radius + altitude0 + height./2;
r02 = r0.^2;
r03 = r0.^3;
delta_r = height;
delta_lonT = delta_lonT*pi/180;
delta_latT = delta_latT*pi/180;
delta_lon = lon0-lon;
cosdelta_lon = cos(delta_lon);
sindelta_lon = sin(delta_lon);
coslat0xcosdelta_lon = coslat0.*cosdelta_lon;
coslatxsinlat0xcosdelta_lon = coslat.*sinlat0.*cosdelta_lon;
sinlatxsinlat0 = sinlat.*sinlat0;
sinlatxcoslat0 = sinlat.*coslat0;
% Constant
G = 6.674215*10^-11;

% Calc help variables
psi0 = acos(sinlatxsinlat0 + coslat.*coslat0.*cos(lon0 - lon));
cospsi0 = cos(psi0);
l0 = sqrt(r2 + r02 - 2*r.*r0.*cospsi0);
l02 = l0.^2;
l03 = l0.^3;
l04 = l02.*l02;
r0dl03 = (r0./l0).^3;

% Calc integral kernels
L000 = (r02.*(r - r0.*cospsi0).*coslat0)./l03;

L200 = ((r.*coslat0)./l03).*(2-((3.*r0)./l02).*(5*r0 - (2*r + 3*r0.*cospsi0).*cospsi0) + ...
        ((15*r03)./l04).*sin(psi0).^2.*(r0 - r.*cospsi0));

L020 = r0dl03.*coslat.*(1 - 2*sinlat02).*cosdelta_lon + ...
        (r02./l0.^5).*(-r.*(r2 + r02).*coslat0 + ...
        r0.*sinlat.*(-r.*r0.*(sinlatxcoslat0 - coslatxsinlat0xcosdelta_lon) + ...
        sinlat0.*coslat0.*(2*r2 + 4*r02 - 3*r.*r0.*sinlatxsinlat0)) + ...
        r02.*coslat.*cosdelta_lon.*(1 - 2*sinlat02).*...
        (r0 + r.*coslat.*coslat0xcosdelta_lon) + ...
        r.*r02.*coslat.*sinlat0.*coslat0xcosdelta_lon.*...
        (3*sinlatxcoslat0 - 4*coslatxsinlat0xcosdelta_lon)) + ...
        ((5*r.*r03)./l0.^7).*(-r.*(r2 + r02).*sinlat0 + ...
        r02.*coslat.*sinlat0.*coslat0xcosdelta_lon.*...
        (r0 + r.*coslat.*coslat0xcosdelta_lon) + ...
        r0.*sinlat.*(2*r2 - r02 - r.*r0.*cospsi0 + sinlat02.*...
        (r2 + 2*r02 - r.*r0.*sinlatxsinlat0))).*...
        (sinlatxcoslat0 - coslatxsinlat0xcosdelta_lon);

L002 = r0dl03.*coslat.*coslat0.^2.*...
        (cosdelta_lon - ((3*r)./l02).*(2*r0.*coslat.*coslat0.*sindelta_lon.^2 + ...
        (r - r0.*cospsi0).*cosdelta_lon) + ...
        ((15*r2.*r0)./l04).*coslat.*coslat0.*(r - r0.*cospsi0).*sindelta_lon.^2);

% Final calc
dg = G.*density.*delta_r.*delta_latT.*delta_lonT.*(L000 + (1/24)*(L200.*delta_r.^2 + L020.*delta_latT.^2 + L002.*delta_lonT.^2));
if sum_switch == 1
    dg = sum(sum(dg));
end

end





    