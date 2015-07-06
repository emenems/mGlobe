function [X,Y,Z] = mGlobe_elip2xyz(L,B,H)
%MGLOBE_ELIP2XYZ Transform ellipsoidal coordinates to 3D cartesian
%   The function transforms the ellipsoidal longitude, latitude and height
%   to 3D cartesian coordinates. The transformation assumes WGS84
%   ellipsoid.
%   Input:
%       L       ... longitude (rad)
%       B       ... latitude (rad)
%       H       ... height (m)
%   Output:
%       X,Y,Z   ... 3D coordinates (m)
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

if nargin == 2
    [r,s] = size(L);
    H(1:r,1:s) = 0;
end
a = 6378137;                                                                % WGS84 major axis (m)
b = 6356752.314245;
e = sqrt((a^2-b^2)/a^2);
N = a./(1-e^2*sin(B).^2).^0.5;

X = (N+H).*cos(B).*cos(L);
Y = (N+H).*cos(B).*sin(L);
Z = (N.*(1-e^2)+H).*sin(B);
end