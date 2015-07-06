function [dgE,dgP,la_out,fi_out,la,fi] = mGlobe_Global(lonD,latD,hranice,rozlisenie,tabulka,r,treshold)
%MGLOBE_GLOBAL calculate global correction
%   Function calculate the gravity response (loading and newton) to 1 mm 
%   of water for given grid using spherical approximation and equation
%   provided by: Farrell, W.E (1972). "Deformation of the Earth by
%   surface loads." Reviews of Geophysics 10(3): 761-797. 
% 
% Input:
%   lonD        ... longitude of point of observation (ellipsoidal, deg)
%   latD        ... latitude of point of observation (ellipsoidal, deg)
%   hranice     ... border of calculation grid (ellipsoidal deg [min(lon),
%                   min(lat);max(lon),max(lon)])
%   rozlisenie  ... grid resolution
%   tabulka     ... table with spherical distance and loading effect for
%                   1 kg of load [psi (deg), dg_loading (m/s^2 / kg)]
%   r           ... radius of replacement sphere (m)
%   treshold    ... threshold of calculation (psi in deg)
%
% Output:
%   dgE         ... gravity loading CORRECTION in m/s^2/ kg/m^2
%   dgP         ... gravity attraction CORRECTION in m/s^2/ kg/m^2
%   la_out      ... longitude: grid of omitted area (deg)
%   fi_out      ... latitude: grid of omitted area (deg)
%   la          ... longitude: created grid of calculation area (deg)
%   fi          ... latitude: created grid of calculation area (deg)
%
%
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

%% Grid preparation
num_step_la = length(hranice(1,1):rozlisenie:hranice(2,1));
num_step_fi = length(hranice(1,2):rozlisenie:hranice(2,2));
la = linspace(hranice(1,1),hranice(2,1),num_step_la);                       % uniformly spaced vector
fi = linspace(hranice(1,2),hranice(2,2),num_step_fi);
delta = [s2rad(abs(la(1)-la(2))) s2rad(abs(fi(1)-fi(2)))];                  % grid resolution
[la,fi] = meshgrid(la,fi);                                                  % meshgrid
fiG = mGlobe_elip2sphere(s2rad(la),s2rad(fi));                              % transform Grid coordinates to sphere
fiDG = mGlobe_elip2sphere(s2rad(lonD),s2rad(latD));                         % transform Point coordinates to sphere
deltaG = abs(fiG+delta(2)/2 - mGlobe_elip2sphere(s2rad(la),s2rad(fi)-delta(2)/2)); % calc new grid resolution in latitude direction
psi = psiSphere(s2rad(lonD),fiDG,s2rad(la),fiG);                            % calc spherical distance

%% Computation
surface = 2*r^2*delta(1).*cos(fiG).*sin(deltaG./2);                         % surface of each grid cell
if ~isempty(tabulka)
    psiint = interp1(tabulka(:,1),tabulka(:,2),psi,'pchip');                % interpolate loading correction
    dgE = psiint.*surface;                                                  % m/s-2
else
    dgE(size(psi)) = 0;
end
% Constant
GM = 3.986004418*10^14;
G = 6.674215*10^-11;
Me = GM/G;
dgP = (-9.80665./(4*Me.*sin(psi./2))).*surface;                             % calc. attraction correction (spherical approximation)

%% Remove area out of current zone
[riad,stlp] = find(psi < s2rad(treshold));
r1 = min([riad stlp]);
r2 = max([riad stlp]);
dgE(r1(1):r2(1),r1(2):r2(2)) = 0;                                           % squared area!! (not circle)
dgP(r1(1):r2(1),r1(2):r2(2)) = 0;
la_out = la(r1(1):r2(1),r1(2):r2(2));                                       % create omitted grid
fi_out = fi(r1(1):r2(1),r1(2):r2(2));

%% SubFunction for deg to rad transformation and spherical distance computation
    function rad = s2rad(s)
        rad = (s*2*pi)./360;
    end

    function psi = psiSphere(laD,fiDG,la,fiG)
        psi = acos(sin(fiG).*sin(fiDG) + cos(fiG).*cos(fiDG).*cos(laD-la));
    end
end

