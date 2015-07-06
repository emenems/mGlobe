function [dgE,dgP,la_out,fi_out,la,fi] = mGlobe_Local(lonD,latD,hD,dem,hranice,rozlisenie,tabulka,r,treshold_in,treshold_out)
%MGLOBE_LOCAL calculation of near zone gravity correction
%   Function calculates the gravity response (loading and newton) to 1 mm 
%   of water for given grid using spherical approximation + digital
%   elevation model. Loading effect is based on: Farrell, W.E (1972). 
%   "Deformation of the Earth by surface loads." Reviews of Geophysics 
%   10(3): 761-797. The attraction part is computed using point mass 
%   approximation. 
% 
% Input:
%   lonD        ... longitude of the point of observation (ellipsoidal,deg)
%   latD        ... latitude of the point of observation (ellipsoidal, deg)
%   hD          ... height of the point of observation (ellipsoidal = 
%                   spherical, m)
%   hranice     ... border/polygon of calculation grid (ellipsoidal deg 
%                   [min(lon),min(lat);max(lon),max(lon)])
%   rozlisenie  ... grid resolution (deg)
%   tabulka     ... table with spherical distance and loading effect for
%                   1 kg of load [psi (deg), dg_loading (m/s^2 / kg)]
%   r           ... radius of replacement sphere (m)
%   treshold_in ... threshold of calculation (psi in deg) for inner zone
%   treshold_out... threshold of calculation (psi in deg) for outer zone
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
num_step_la = length(hranice(1,1):rozlisenie:hranice(2,1));                 % determine the future grid dimensions
num_step_fi = length(hranice(1,2):rozlisenie:hranice(2,2));
la = linspace(hranice(1,1),hranice(2,1),num_step_la);                       % create new uniform vectors used for grid
fi = linspace(hranice(1,2),hranice(2,2),num_step_fi);
delta = [s2rad(abs(la(1)-la(2))) s2rad(abs(fi(1)-fi(2)))];                  % determine new grid resolution (ellipsoid, rad)
[la,fi] = meshgrid(la,fi);                                                  % create new grid on ellipsoid
fiG = mGlobe_elip2sphere(s2rad(la),s2rad(fi));                              % transform latitude to sphere
fiDG = mGlobe_elip2sphere(s2rad(lonD),s2rad(latD));                         % transform point of observation to sphere
deltaG = abs(fiG+delta(2)/2 - mGlobe_elip2sphere(s2rad(la),s2rad(fi)-delta(2)/2)); % calc. new grid resolution on sphere

%% Computation
psi = psiSphere(s2rad(lonD),fiDG,s2rad(la),fiG);clear fiDG                  % compute spherical distance between grid cells and point of observation
surface = 2*r^2*delta(1).*cos(fiG).*sin(deltaG./2); clear deltaG            % surface of each grid cell
% Indirect effect = loading
if ~isempty(tabulka)
    psiint = interp1(tabulka(:,1),tabulka(:,2),psi,'pchip');                % interpolate the loading effect (faster than Legendre polynomial evaluation)
    dgE = psiint.*surface;clear psiint                                      % 
else
    dgE(size(psi)) = 0;
end
% Direct effect = newton's attraction
G = 6.674215*10^-11;
if isempty(dem)
    hreg = 0;
else
    hreg = interp2(dem.lon,dem.lat,dem.height,la,fi,'linear');              % include digital elevation model
end
[Xr,Yr,Zr] = mGlobe_elip2xyz(s2rad(la),s2rad(fi),hreg);                     % transform to X,Y,Z coordinates (to include DEM)
[Xd,Yd,Zd] = mGlobe_elip2xyz(s2rad(lonD),s2rad(latD),hD);
dist_reg2 = (Xd-Xr).^2 + (Yd-Yr).^2 + (Zd-Zr).^2;
dist_reg = sqrt(dist_reg2);
clear Xr Yr Zr Xd Yd Zd
dgP0 = G*surface./dist_reg2;clear surface                                 % in the direction of the connection line
cos_alpha = (dist_reg2 + (r + hD).^2 - (r + hreg).^2)./(2*dist_reg.*(r + hD));clear dist_reg
dgP = -dgP0.*cos_alpha;                                                     % in the direction of plumb-line
clear dgP0 cos_alpha

%% Remove area out of current zone
if isempty(treshold_out)                                                    % squared area (not circle)
    [riad,stlp] = find(psi <= s2rad(treshold_in));
    r1 = min([riad stlp]);
    r2 = max([riad stlp]);
    if ~isempty(r1) && length(r1) > 1 && length(r2) > 1
        dgE(r1(1):r2(1),r1(2):r2(2)) = 0;
        dgP(r1(1):r2(1),r1(2):r2(2)) = 0;
        la_out = la(r1(1):r2(1),r1(2):r2(2));
        fi_out = fi(r1(1):r2(1),r1(2):r2(2));
    else
        la_out = [];
        fi_out = [];
    end
else                                                                        % circle area
    dgE(psi <= s2rad(treshold_in)) = 0;
    dgP(psi <= s2rad(treshold_in)) = 0;
    dgE(psi > s2rad(treshold_out)) = 0;
    dgP(psi > s2rad(treshold_out)) = 0;
    [riad,stlp] = find(psi <= s2rad(treshold_in));
    r1 = min([riad stlp]);
    r2 = max([riad stlp]);
    if ~isempty(r1) && length(r1) > 1 && length(r2) > 1
        la_out = la(r1(1):r2(1),r1(2):r2(2));
        fi_out = fi(r1(1):r2(1),r1(2):r2(2));
    else
        la_out = [];
        fi_out = [];
    end
end
%% SubFunction for deg to rad transformation and spherical distance computation
    function rad = s2rad(s)
        rad = (s*2*pi)./360;
    end

    function psi = psiSphere(laD,fiDG,la,fiG)
        psi = acos(sin(fiG).*sin(fiDG) + cos(fiG).*cos(fiDG).*cos(laD-la));
    end
end

