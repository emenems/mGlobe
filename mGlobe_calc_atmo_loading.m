function [dgE,la_out,fi_out,la,fi,zone_id] = mGlobe_calc_atmo_loading(lonD,latD,hranice,rozlisenie,tabulka,radius,treshold_in,treshold_out)
%MGLOBE_CALC_ATMO_LOADING Function calculates the atmospheric loading effect
%   Calculation based on: Merriam, J. B. (1992). 
%   "Atmospheric pressure and gravity." Geophysical Journal International 
%   109(3): 488-500.
% 
% INPUT:    
%   lonD          ...     longitude of the calculation point (deg)
%   latD          ...     latitude of the calculation point (deg)
%   hranice       ...     polygon of the computation zone (deg)
%   rozlisenie    ...     grid resolution (deg)
%   tabulka       ...     given table with spherical distance and loading
%                         effect
%   radius        ...     radius of the approximation sphere (m)
%   treshold_in   ...     minimal spherical distance to the point of
%                         calculation (deg)
%   treshold_out  ...     maximal spherical distance to the point of
%                         calculation (deg)
% 
% OUTPUT:  
%   dgE           ...     loading gravity effect (m/s^2/hPa)
%   la_out        ...     omitted output grid (deg)
%   fi_out        ...     omitted output grid (deg)
%   la            ...     used computation longitude grid (deg)
%   fi            ...     used computation latitude grid (deg)
%   zone_id       ...     identification matrix for used/omitted points of
%                         the la/fi grid (logical)
% 
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

%% Create new grid
num_step_la = length(hranice(1,1):rozlisenie:hranice(2,1));                 % determine the future grid dimensions (longitude)
num_step_fi = length(hranice(1,2):rozlisenie:hranice(2,2));                 % determine the future grid dimensions (latitude)
la = linspace(hranice(1,1),hranice(2,1),num_step_la);                       % create new uniform vectors used for grid
fi = linspace(hranice(1,2),hranice(2,2),num_step_fi);
delta = [s2rad(abs(la(1)-la(2))) s2rad(abs(fi(1)-fi(2)))];                  % ellipsoidal resolution (rad)
[la,fi] = meshgrid(la,fi);
fiG = mGlobe_elip2sphere(s2rad(la),s2rad(fi));                              % latitude on sphere (grid rad)
fiDG = mGlobe_elip2sphere(s2rad(lonD),s2rad(latD));                         % latitude on sphere (POINT rad)
% Calc new grid resolution
deltaG = fiG - circshift(fiG,[1 0]);                                        % approx. grid resolution in latitude direction (on sphere)
deltaG(1,:) = deltaG(2,:);
deltaG(end,:) = deltaG(end-1,:);
% Calc spherical distance
psi = psiSphere(s2rad(lonD),fiDG,s2rad(la),fiG);clear fiDG                  % spherical distance between grid cells and point of observation

%% Estimation of indirect == deformation part (dgE) FOR 1hPa on SURFACE
psiint = interp1(tabulka(:,1),tabulka(:,2),psi,'pchip');                    % tabled values in [rad, uGal/hPa]
area = 2*radius^2*delta(1).*cos(fiG).*sin(deltaG./2); clear deltaG          % surface area (on the sphere in m^2)
dgE = (psiint.*(area./radius^2))./(10^5*psi*2*pi*(1 - cosd(1)));            % uGal/hPa (surface in m^2 transformed to steradians)
dgE = dgE*10^-8;                                                            % m/s^2/hPa
clear psiint

%% Remove values out of current zone
zone_id = dgE;
zone_id(:,:) = 1;
if isempty(treshold_out)                                                    % for outer zones
    [riad stlp] = find(psi <= s2rad(treshold_in));
    r1 = min([riad stlp]);
    r2 = max([riad stlp]);
    if ~isempty(r1) && length(r1)+length(r2) > 2
        zone_id(r1(1):r2(1),r1(2):r2(2)) = 0;
        la_out = la(r1(1):r2(1),r1(2):r2(2));
        fi_out = fi(r1(1):r2(1),r1(2):r2(2));
    else
        la_out = [];
        fi_out = [];
    end
else
    zone_id(psi <= s2rad(treshold_in)) = 0;                                 % for inner zones
    zone_id(psi > s2rad(treshold_out)) = 0;
    [riad stlp] = find(psi <= s2rad(treshold_in));
    r1 = min([riad stlp]);
    r2 = max([riad stlp]);
    if ~isempty(r1) && length(r1)+length(r2) > 2
        la_out = la(r1(1):r2(1),r1(2):r2(2));
        fi_out = fi(r1(1):r2(1),r1(2):r2(2));
    else
        la_out = [];
        fi_out = [];
    end
end
dgE(zone_id == 0) = 0;

%% SubFunction for deg to rad transformation and spherical distance computation
    function rad = s2rad(s)
        rad = (s*2*pi)./360;
    end

    function psi = psiSphere(laD,fiDG,la,fiG)
        psi = acos(sin(fiG).*sin(fiDG) + cos(fiG).*cos(fiDG).*cos(laD-la));
    end
end

