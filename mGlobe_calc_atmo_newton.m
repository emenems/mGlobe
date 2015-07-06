function dgP = mGlobe_calc_atmo_newton(lonD,latD,hd,lon0,lat0,surface,height,density,zone_id,radius,approx,max_level,h0_surface)
%MGLOBE_CALC_ATMO_NEWTON Function calculates the atmospheric attraction effect
%   Calculation based on:
%   Tesseroid approx.: Heck, B. and K. Seitz (2007). 
%   "A comparison of the tesseroid, prism and point-mass approaches for 
%   mass reductions in gravity field modelling." Journal of Geodesy 81(2): 
%   121-136.
%   Point mass approx.: Farrell, W.E (1972). "Deformation of the Earth by
%   surface loads." Reviews of Geophysics 10(3): 761-797.
%
% INPUT:    
%   lonD          ...     longitude of the calculation point (deg)
%   latD          ...     latitude of the calculation point (deg)
%   hd            ...     height of the calculation point (m)
%   lon0          ...     longitude of the computation grid (deg)
%   lat0          ...     latitude of the computation grid (deg)
%   surface       ...     structure variable with surface pressure,
%                         temperature, humidity and density
%   height        ...     model pressure level heights 
%   density       ...     model density (kg/m^3)
%   zone_id       ...     identification matrix for used/omitted points of
%                         the lon0/lat0 grid (logical)
%   radius        ...     radius of the approximation sphere (m)
%   approx        ...     approximation switch (2D|3D|3Dlocal)
%   max_level     ...     maximal level of computation (1-37)
%   h0_surface    ...     orography (m)
%   h0_dem        ...     digital elevation model used for 3Dlocal only
%                         (h0_dem.height m) + surface air density 
%                         (h0_dem.density kg/m^3) - only for beta version
% 
% OUTPUT:  
%   dgP           ...     attraction gravity effect (m/s^2)
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                06.04.2015
%                                                                      v1.0


%% Prepare grid
fiG = mGlobe_elip2sphere(s2rad(lon0),s2rad(lat0))*180/pi;                   % latitude on sphere (grid deg)
fiDG = mGlobe_elip2sphere(s2rad(lonD),s2rad(latD))*180/pi;                  % latitude on sphere (POINT deg)
% Calc new grid resolution
delta_latT = fiG - circshift(fiG,[1 0]);                                    % approx. grid resolution in latitude direction (on sphere)
delta_latT(1,:) = delta_latT(2,:);
delta_latT(end,:) = delta_latT(end-1,:);
delta_lonT = lon0(1,2)-lon0(1,1);

switch approx
    case '2D'                                                               % Compute the atmospheric effect for spherical approximation
%% Estimation of direct == newtonian attraction (dgP) for thin layer
GM = 3.986004418*10^14;                                                     % set constants
G = 6.674215*10^-11;
Me = GM/G;
psi = psiSphere(s2rad(lonD),s2rad(fiDG),s2rad(lon0),s2rad(fiG));clear fiDG  % spherical distance to grid cells
area = 2*radius^2*s2rad(delta_lonT).*cos(s2rad(fiG)).*sin(s2rad(delta_latT)./2); % surface (on the sphere in m^2)
dgP = (9.80665./(4*Me.*sin(psi./2)));                                       % effect in m/s^2 for 1 kg
dgP = (dgP.*area)/9.80665;                                                  % for 1 Pascal                               
dgP(zone_id == 0) = 0;                                                      % remove areas out of current zone

    case '3D'                                                               % compute the 3D tesseroid approximation (no local DEM used)
%% Estimation of direct == newtonian attraction (dgP) for 3D layers
i = max_level;                                                              % for the first layer = surface layer
h_lower = h0_surface;                                                       % height of the lower tesseroid boundary
h_upper = mGlobe_interpolation(surface.lon,surface.lat,height.data(:,:,i),lon0,lat0,0); % height of the upper boundary
ro_lower = mGlobe_interpolation(surface.lon,surface.lat,surface.density,lon0,lat0,0); % density at the lower boundary
ro_upper = mGlobe_interpolation(surface.lon,surface.lat,density.data(:,:,i),lon0,lat0,0); % density at the upper boundary
tess_density = (ro_lower+ro_upper)/2;                                       % mean tesseroid density 
tess_density(h_upper<h_lower) = 0;                                          % set density to zero if the tesseroid is below the orography
tess_density(zone_id == 0) = 0;                                             % remove areas out of the current zone                                  
ro_upper(h_upper<h_lower) = ro_lower(h_upper<h_lower);                      % shift the boundary to the next layer
h_upper(h_upper<h_lower) = h_lower(h_upper<h_lower);                        % shift the boundary to the next layer                                 
h_upper(isnan(tess_density)) = h_lower(isnan(tess_density));                % shift the boundary to the next layer (for points with NaN -> MERRA only). Shift only upper boundary, not density as that might be NaN for pressure levels, but not for orography.
tess_density(isnan(tess_density)) = 0;                                      % remove all points with NaNs = MERRA only
tess_height = h_upper - h_lower;                                            % tesseroid height (keep in mind that the density is set to zero where h_lower > h_upper)
dgP(i) = mGlobe_tesseroid([lonD,fiDG,hd],lon0,fiG,h_lower,tess_height,tess_density,abs(delta_lonT),abs(delta_latT),radius,1); % gravitational effect of the first layer
while i >= 2                                                                % Loop for other layers (same procedure as for the first layer)
    i = i - 1;                                                              % same steps as for the first layer (i.e. between orography - 1000 hPa)
    h_lower = h_upper;                                                      % move to next level
    ro_lower = ro_upper;
    h_upper = mGlobe_interpolation(surface.lon,surface.lat,height.data(:,:,i),lon0,lat0,0);
    ro_upper = mGlobe_interpolation(surface.lon,surface.lat,density.data(:,:,i),lon0,lat0,0);
    tess_density = (ro_lower+ro_upper)/2;
    tess_density(h_upper<h_lower) = 0;
    tess_density(zone_id == 0) = 0;
    ro_upper(h_upper<h_lower) = ro_lower(h_upper<h_lower);
    h_upper(h_upper<h_lower) = h_lower(h_upper<h_lower);
    h_upper(isnan(tess_density)) = h_lower(isnan(tess_density)); 
    tess_density(isnan(tess_density)) = 0;
    tess_height = h_upper - h_lower;
    dgP(i) = mGlobe_tesseroid([lonD,fiDG,hd],lon0,fiG,h_lower,tess_height,tess_density,abs(delta_lonT),abs(delta_latT),radius,1);
end

end
%% SubFunctions: degrees to radians + spherical distance
    function rad = s2rad(s)
        rad = (s*2*pi)./360;
    end

    function psi = psiSphere(laD,fiDG,la,fiG)
        psi = acos(sin(fiG).*sin(fiDG) + cos(fiG).*cos(fiDG).*cos(laD-la));
    end
end

