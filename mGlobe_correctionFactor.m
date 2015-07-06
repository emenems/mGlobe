function mGlobe_correctionFactor(Input,dem_file,orography_file,mGlobe_output,local_pressure_tsf,local_pressure_channel)
% MGLOBE_CORRECTIONFACTOR computes the correction factor
% Function serves for the computation of site- and model-dependent
% correction/admittance factor. This factor should be use for the multiplication
% of pressure residuals (in situ - ERA Interim). The corresponding effect must
% be then added to the mGlobe atmospheric effect in order to reduce the
% effect of deficient spatial and temporal resolution.
% This function requires Matlab R2012a or higher including Mapping and
% Statistics toolboxes + 8 GB of RAM. Additionally, this function uses 
% several functions of the mGlobe toolbox and therefore should be stored in 
% the same folder. To reduce the resolution of the computation grid, modify
% row 56.
% 
% Inputs:
%   Input           ...     Coordinates of the gravimeter (WGS84 + height)
%                           [Latitude (deg), Longitude (deg), Height (m)]
%   dem_file        ...     Matlab structure array file containing digital
%                           elevation model (dem.lon,dem.lat,dem.height).
%                           This file can by generated using the mGlobe
%                           build-in conversion tool. Loading of a DEM with
%                           insufficient extension will result in NaNs or
%                           correction factor equal to zero. Maximal
%                           integration area = 0.1 deg (spherical distance)
%                           The difference between the Input Height and
%                           interpolated DEM height is used to shift the
%                           whole DEM "towards" the gravimeter. If you do
%                           not want to shift the DEM (e.g. underground 
%                           gravimeter) set Input(3) = NaN.
%   orography_file  ...     NetCDF (ERA) or HDF (MERRA) file containing the 
%                           orography. Same as for the computation of 
%                           mGlobe atmospheric effect.
%   mGlobe_output   ...     txt mGlobe output file (mGlobe-Atmo).
%   local_pressure_tsf.     TSF file containing local (in situ) pressure
%                           variations in hPa.
%   local_pressure_channel  TSF channel containing the pressure variations
%                           in hPa (channel not the column).
% 
% Outputs:
%   The result will be printed to the command line. Additional plots will
%   be used to visualize the results.
% 
% Example:
%   Input = [47.9283 15.8598 1044];
%   dem_file = 'CO_DEM_upTo_01deg_WGS84.mat'; 
%   orography_file = '\ERAinterim\Invariant\ERA_INVARIANT_ALL.nc';
%   mGlobe_output = 'F:\mikolaj\Documents\MYpapers\mGlobe\Results\ATMO\mGlobE\CO\CO_MGLOBE_ATMO_ERA_012010_122013_6h_dem0.txt'; 
%   local_pressure_tsf = 'CO_SG_g_and_p.tsf'; 
%   local_pressure_channel = 2;
%   mGlobe_correctionFactor(Input,dem_file,orography_file,mGlobe_output,local_pressure_tsf,local_pressure_channel);
%   
%                                                   M. Mikolaj, 18.11.2014
%                                                   mikolaj@gfz-potsdam.de

%% Main settings
psi = 0.10;                                                                 % fixed integration radius
resol = 0.0002;                                                             % fixed spatial resolution

%% Load input DEM and mGlobe output
fprintf('Reading DEM...\n');
dem = importdata(dem_file);
dem.height = dem.height;
Hdem = interp2(dem.lon,dem.lat,dem.height,Input(2),Input(1));               % height of the compuation point with respect to orography
if ~isnan(Input(3))
    dem.height = dem.height + (Input(3)-Hdem);                                  % shift the DEM towards the actual height.
end
meteo = load(mGlobe_output);    
time = meteo(:,1);                                                          % time
p_era = meteo(:,9);                                                         % pressure in Pa
t = meteo(:,10);                                                            % temperature in K
q = meteo(:,11);                                                            % specific humidity in kg/kg
Input(4) = mGlobe_elip2sphere(Input(2)*pi/180,Input(1)*pi/180,Hdem)*180/pi; % spherical latitude (deg)

%% Load local pressure data
fprintf('Reading TSF...\n');
fid = fopen(local_pressure_tsf);
tline = fgetl(fid);header_lines = 0;
if ~isempty(tline)
while isnan(str2double(tline(1)))                                           % determine number of header lines
    tline = fgetl(fid);
    try
        if strcmp(char(tline(2:9)),'UNDETVAL')                              % try to find flagged values
            flagg = tline(11:end);
            flagg = str2double(flagg);
        end
    end
    if isempty(tline)
        tline = 'mGlobe_RULEZ';
    end
    header_lines = header_lines + 1;                                        % count header lines
end
fclose(fid);                                                                % close the file to open it with known number of header lines
grav = dlmread(local_pressure_tsf,'',header_lines,0);                       % warning no footer info are allowed
grav(grav(:,local_pressure_channel+6)==flagg,local_pressure_channel+6) = NaN; % remove flagged values
p_obs = interp1(datenum(grav(:,1:6)),grav(:,local_pressure_channel+6)*100,time); % interpolate the local pressure+convert to Pa

%% Load Orography
fprintf('Reading Orography...\n');
switch orography_file(end-1:end)
    case 'nc'
        ncid_ref = netcdf.open(orography_file,'NC_NOWRITE');                        % open given netCDF files
        [ndims,nvars] = netcdf.inq(ncid_ref);                                       % get variable names
        for nc = 1:nvars
            varname = netcdf.inqVar(ncid_ref,nc-1);   
            switch varname
                % MERRA
                case 'XDim'
                    orography.lon = netcdf.getVar(ncid_ref,nc-1,'double');          % get longitude
                case 'YDim'
                    orography.lat = netcdf.getVar(ncid_ref,nc-1,'double');          % get latitude
                case 'PHIS'
                    orography.data = netcdf.getVar(ncid_ref,nc-1,'double');         % get orography
                    try
                        scale_factor = netcdf.getAtt(ncid_ref,nc-1,'scale_factor'); % check if scaling factor does exist
                    catch
                        scale_factor = 1;
                    end
                    try
                        add_offset = netcdf.getAtt(ncid_ref,nc-1,'add_offset');     % check if offset does exist
                    catch
                        add_offset = 0;
                    end
                    orography.data = orography.data'*scale_factor + add_offset;     % scale and transpose the orography matrix
                    orography.data = double(orography.data)./9.80665;               % transform to meters using CONSTANT gravity (no latitude dependency!)
                    [orography.lon,orography.lat] = meshgrid(double(orography.lon),double(orography.lat)); % vectors to matrices
                % ERA
                case 'longitude'
                    orography.lon = netcdf.getVar(ncid_ref,nc-1,'double');          % get longitude
                case 'latitude'
                    orography.lat = netcdf.getVar(ncid_ref,nc-1,'double');          % get latitude
                case 'z'
                    orography.data = netcdf.getVar(ncid_ref,nc-1,'double');         % get orography
                    try
                        scale_factor = netcdf.getAtt(ncid_ref,nc-1,'scale_factor'); % check if scaling factor does exist
                    catch
                        scale_factor = 1;
                    end
                    try
                        add_offset = netcdf.getAtt(ncid_ref,nc-1,'add_offset');     % check if offset does exist
                    catch
                        add_offset = 0;
                    end
                    orography.data = orography.data'*scale_factor + add_offset;     % scale and transpose the orography matrix
                    orography.data = orography.data./9.80665;                       % transform to meters using CONSTANT gravity (no latitude dependency!)
                    [orography.lon,orography.lat] = meshgrid(orography.lon,orography.lat);      % create lon/lat matrices
                    clear ndims nvars nc add_offset scale_factor                                % remove used variables
                    orography.lon(orography.lon>=180) = orography.lon(orography.lon>=180) - 360; % transform coordinates to (-180,180) system  
                    ri = find(abs(diff(orography.lon(1,:)))==max(abs(diff(orography.lon(1,:)))));
                    orography.lon = horzcat(orography.lon(:,ri+1:end),orography.lon(:,1:ri));   % Connect matrices to remove discontinuity
                    orography.lat = horzcat(orography.lat(:,ri+1:end),orography.lat(:,1:ri));
                    orography.data = horzcat(orography.data(:,ri+1:end,:),orography.data(:,1:ri,:)); 
            end
        end
        netcdf.close(ncid_ref);                                                       % geopotential to height
    case 'df'                                                               % MERRA only
        orography.data = hdfread(orography_file,'PHIS');                    % read hdf file containing MERRA orography
        orography.data(orography.data>1000000) = NaN;                       % set flagged values to NaN (no flagged values are expected though)
        orography.data = double(orography.data)./9.80665;                   % convert to m  
        orography.lon = hdfread(orography_file,'XDim:EOSGRID');             % longitude vector
        orography.lat = hdfread(orography_file,'YDim:EOSGRID');             % latitude vector
        [orography.lon,orography.lat] = meshgrid(double(orography.lon),double(orography.lat)); % vectors to matrices
end

%% Create grid
fprintf('Preparing computation...\n');
psi_lon = 0.05;delta_lon = 0.05;
while psi_lon<psi
   delta_lon = delta_lon + 0.02;
   psi_lon = acos(sin(Input(4)*pi/180).*sin(Input(4)*pi/180) + cos(Input(4)*pi/180).*cos(Input(4)*pi/180).*cos((delta_lon)*pi/180))*180/pi; 
end
psi_lat = 0.05;delta_lat = 0.05;
while psi_lat<psi
   delta_lat = delta_lat + 0.02;
   psi_lat = acos(sin(Input(4)*pi/180).*sin((Input(4)-delta_lat)*pi/180) + cos(Input(4)*pi/180).*cos((Input(4)-delta_lat)*pi/180).*cos(0))*180/pi; 
end
[comp.lon,comp.lat] = meshgrid(Input(2)-delta_lon-resol/2:resol:Input(2)+delta_lon+resol/2,... % computation grid for spherical coordinated
                               Input(4)-delta_lat-resol/2:resol:Input(4)+delta_lat+resol/2); 
% Transform to sphere
a = 6378137;                                                                % ellipsoidal major axis (m)
b = 6356752.314245;
e = sqrt((a^2-b^2)/a^2);e2c = (a^2-b^2)/b^2;
R = sqrt(1+2/3*e^2+3/5*e^4+4/7*e^6+5/9*e^8+6/11*e^10+7/13*e^12)*b;          % Radius of the replacement sphere -> equal surface
X = (R).*cos(comp.lat*pi/180).*cos(comp.lon*pi/180);                        % approximated solution (zero heights)                        
Y = (R).*cos(comp.lat*pi/180).*sin(comp.lon*pi/180);
Z = (R).*sin(comp.lat*pi/180);
p = sqrt(Y.^2+X.^2);
tet = atan((Z.*a)./(p.*b));
comp.lat_ell = atan((Z+e2c*b*(sin(tet)).^3)./(p-e^2*a*(cos(tet)).^3));
clear tet p X Y Z e a b e2c                                               % remove used variables
comp.lat_ell = comp.lat_ell*180/pi;
comp.dem = interp2(dem.lon,dem.lat,dem.height,comp.lon,comp.lat_ell);       % get topography heights
comp.oro = interp2(orography.lon,orography.lat,orography.data,comp.lon,comp.lat_ell);          % get orography heights
comp.psi = acos(sin(comp.lat*pi/180).*sin(Input(4)*pi/180) + cos(comp.lat*pi/180).*cos(Input(4)*pi/180).*cos((Input(2)-comp.lon)*pi/180)); % spherical distance
Hdem = interp2(dem.lon,dem.lat,dem.height,Input(2),Input(1));               % height of the compuation point with respect to orography
Horo = interp2(orography.lon,orography.lat,orography.data,Input(2),Input(1)); % height of the compuation point with respect to orography
height = comp.dem-comp.oro;                                                 % tesseroid heights
clear dem orography;                                                        % remove used variables

%% Determine the density of each tesseroid
lower_boundary = comp.oro;                                                  % create lower boundary variable
lower_boundary(height<0) = comp.dem(height<0);                              % new lower boundary for tesseroids where orography is above topography
dens = height; dens(:,:) = 1;                                               % create density variable
dens(height<0) = -1;                                                        % negative density for tesseroids below the orography
dens(comp.psi*180/pi>psi) = 0;                                              % remove tesseroids outside the local zone
comp.dem = [];

%% Show altitude differences
% figure
% mesh(comp.lon,comp.lat,comp.oro);view(0,90);colorbar;
% title('orography');
% 
% figure
% mesh(comp.lon,comp.lat,comp.dem);view(0,90);colorbar;
% title('DEM');
% 
% figure
% mesh(comp.lon,comp.lat,lower_boundary);view(0,90);colorbar;
% title('Lower boundary');
% 
% figure
% mesh(comp.lon,comp.lat,dens);view(0,90);colorbar;
% title('density matrix');
% comp.psi = []; comp.dem = [];
% 
% figure
% mesh(comp.lon,comp.lat,height);view(0,90);colorbar;
% title('DEM-orography');

%% Compute gravitational effect
fprintf('Computing graivty effect (DEM-Orography)...\n');
dg = mGlobe_tesseroid([Input(2),Input(4),Hdem],comp.lon,comp.lat,lower_boundary,abs(height),dens,resol,resol,R,0);
dg = dg*1e+9;                                                               % convert to nm/s^2
dg(comp.psi*180/pi>psi) = 0;                                                % set points outside the integration radius to zero
clear dens;
%% Gravitational effect for each time epoch
effect(1:length(time),1) = NaN;
p_rep(1:length(time),1) = NaN;
clc
for i = 1:1:length(time)
    fprintf('Computing correction factor... (%3.0f%%)\n',round((i/length(time)*100)));
    p_down = p_era(i)*exp((9.81*(comp.oro-lower_boundary))./(287*t(i)*(1+0.608*q(i)))); % pressure for lower boundary
    p_up =   p_era(i)*exp((9.81*(comp.oro-(lower_boundary+abs(height))))./(287*t(i)*(1+0.608*q(i)))); % pressure for upper boundary
    p_rep(i) = p_era(i)*exp((9.81*(Horo-Hdem))./(287*t(i)*(1+0.608*q(i)))); % pressure estimation for gravimeter position
    t_down = t(i)+(comp.oro-lower_boundary)*-0.0065;                        % temperature at the lower boundary
    t_up = t(i)+(comp.oro-(lower_boundary+abs(height)))*-0.0065;            % temperature at the upper boundary
    density_up = p_up./(287.*t_up.*(1-q(i)+q(i)/0.62197));                  % density at the upper boundary
    density_down = p_down./(287.*t_down.*(1-q(i)+q(i)/0.62197));            % density at the lower boundary
    density_up(comp.psi*180/pi>psi) = 0;                                    % set points outside the integration radius to zero
    density_down(comp.psi*180/pi>psi) = 0;                                  % set points outside the integration radius to zero
    clear p_up p_down
    effect(i,1) = sum(sum(dg.*abs((density_up+density_down)/2)));           % total gravity effect for the mean density
    clc
end

%% Estimate the correction factor
delta_p = (p_obs-p_era)/100;
x = delta_p - mean(delta_p(~isnan(delta_p)));                             % x values = observed pressure - orography
if Hdem > Horo
    y = effect - mean(effect(~isnan(effect)));                                           % y = function of x
else
    y = - effect + mean(effect(~isnan(effect)));                                         % y = - effect because of the sign switch related to missing+above lying air
end
r = find(isnan(x)|isnan(y));
if ~isempty(r)
    x(r) = [];y(r) = [];
end
poly = polyfit(x,y,1);                                                      % compute slope
correction_factor = poly(1);                                                  
out_fit = polyval(poly,x);                                                  % fitted dp vs dg

%% Plot
figure
subplot(14,1,1:4);
plot(x,y,'k.',x,out_fit,'r');
xlabel('pressure residuals (hPa)');
ylabel('gravity effect (nm s^{-2})');
legend('Results','Fitted');
title('Pressure residuals vs. gravitational effect');

subplot(14,1,7:9)
plot(time,p_obs/100,'k',time,p_era/100,'b',time,p_rep/100,'r');
xlabel('time (matlab)');
ylabel('pressure (hPa)');
legend('Observed','Orography','Estimated');
set(gca,'XTickLabel',[]);
grid on

subplot(14,1,11:14)
plotyy(time,p_rep/100-p_era/100,time,p_rep/100-p_obs/100);
% yl1 = get(h(1),'YLim');
% delta_p0 = p_rep/100-p_obs/100;
% set(h(2),'YLim',[round(mean(delta_p(~isnan(delta_p0)))-(range(yl1))/2),round(mean(delta_p(~isnan(delta_p0)))+(range(yl1))/2)]);
xlabel('time (matlab)');
ylabel('pressure (hPa)');
legend('Estimated-Orography','Estimated-Observed');
grid on

%% Output
if correction_factor < -2 && correction_factor > -4
    fprintf('Estimated correction factor = %5.2f nm/s^2/hPa\n',correction_factor);
else
    t = meteo(:,4); t = t - mean(t(~isnan(t)));
    f = meteo(:,9)/100; f = f - mean(f(~isnan(f)));
    c2 = polyfit(f,t,1);
    fprintf('The estimated correction factor (=%5.2f) lies outside the expected range, i.e. <-4,-2> nm/s^2/hPa.\nIt is not recommended to use such correction factor. Please check the plotted results.\nThe correct computation requires good performance of the estimated pressure, i.e. the estimated pressure must be similar to the observed air pressure.\nAdditionally, the height interpolated from the loaded DEM must be very similar to the real gravimeter altitude.\nMake sure that the input time series (mGlobe output) are at least 2 years long and sampled every 6 hours.\nThis algorithm fails for areas where the ERA Interim orography is to close to the actual topography (see paper for details).\nValues around 0 nm/s^2/hPa may indicate that the volume above the orography cancels out the effect of the volume below the orography.\nWe recommend to use correction factor %5.2f /s^2/hPa (estimated from least square adjustment of mGlobe total effect and ERA air pressure).\nAlternatively, use site-specific correction factor determined for frequencies higher than 4 cpd.\nOtherwise, the neglected deficient temporal resolution may introduce an error.\n',correction_factor,c2(1));
end

end