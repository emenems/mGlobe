function [fiG] = mGlobe_elip2sphere(la,fi,h)
%MGLOBE_ELIP2SPHERE Transform ellipsoidal latitude to spherical
%   The function transforms the ellipsoidal latitude to spherical,
%   longitude stays the same. The transformation assumes the same 3D
%   cartesian coordinates (X,Y,Z) for both, ellipsoidal and spherical case.
%   Input:
%       la      ... longitude (rad)
%       fi      ... latitude (rad)
%       h       ... height (m)
%   Output:
%       fiG     ... spherical latitude (rad)
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

if nargin == 2
    [X,Y,Z] = mGlobe_elip2xyz(la,fi);
else
    [X,Y,Z] = mGlobe_elip2xyz(la,fi,h);
end
fiG = atan(Z./(X.^2+Y.^2).^0.5);
end

