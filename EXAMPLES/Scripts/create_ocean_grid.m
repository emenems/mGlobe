%% Create reference grid used for Land/Ocean identification
% This grid will be used for computation of GHE, NTOL and ATMO effects
% WARNING: the computation may take several hours depending on the output
% resolution!! 
% Octave version requires IO and MAPPING packages
close all
clear
clc

%% Settings
% Set output grid. Output grid will show 1 in area of oceans/lakes and 0
% over continents
output_file = 'f:\mikolaj\code\libraries\mGlobe\mGlobe_DATA_OceanGrid.mat';
output_resolution = 0.2; % in degrees!

% Input shape with continents (must be a Shapefile!).
% Go to http://www.naturalearthdata.com/downloads/ and download 'land' and
% 'lakes' shapefiles in desired resolution. It is not recommended to use
% resolution higher than 1:50m
% Set lakes to [] for no lake identification. 
input_continent = 'f:\mikolaj\data\global_model\boarders\naturalearthdata\50m_physical\ne_50m_land.shp';
input_lakes = 'f:\mikolaj\data\global_model\boarders\naturalearthdata\50m_physical\ne_50m_lakes.shp';

% Show output (1=yes)
show_grid = 1;

%% Load data
v = version;
if ~strcmp(v(end),')')
    pkg load io
    pkg load mapping
end
boarders = shaperead(input_continent,'UseGeoCoords',true);
if ~isempty(input_lakes)
    lakes = shaperead(input_lakes,'UseGeoCoords',true);
end

%% Create grid
mz = 0;mz2=1;
cc = 0;cc2=0;
[oceans.lon,oceans.lat] = meshgrid(-180+output_resolution/2:output_resolution:180-output_resolution/2,...
                                   -90+output_resolution/2:output_resolution:90-output_resolution/2);
oceans.id = oceans.lon.*0+1;
% Go through all continents/islands and find point within the polygon given
% by input shapefile
for i = 1:length(boarders)
    if strcmp(v(end),')')
        fprintf('Land ID remaining: %d\n',length(boarders)-i);
    end
    % Check if the polygon contains NaN and remove it (required for octave)
    temp_lon = [boarders(i).Lon];
    temp_lat = [boarders(i).Lat];
    temp = find(isnan(temp_lon+temp_lat));
    if length(temp)>1% && ~strcmp(v(end),')')
        % First part = main polygon, other parts after NaN = remove
        % from inside of polygon => revert sign
        oceans.id = oceans.id - double(inpolygon(oceans.lon,oceans.lat,temp_lon(1:temp(1)-1),temp_lat(1:temp(1)-1)));
        for j = 2:length(temp)
            c = temp(j-1)+1;
            oceans.id = oceans.id + double(inpolygon(oceans.lon,oceans.lat,temp_lon(c:temp(j)-1),temp_lat(c:temp(j)-1)));
        end
        clear c j
    elseif length(temp) == 1 && temp~= length(temp)
        % Special case for Euro-asia + Caspian Sea
        oceans.id = oceans.id - double(inpolygon(oceans.lon,oceans.lat,temp_lon(1:temp-1),temp_lat(1:temp-1)));
        oceans.id = oceans.id + double(inpolygon(oceans.lon,oceans.lat,temp_lon(temp+1:end),temp_lat(temp+1:end)));
    else
        oceans.id = oceans.id - double(inpolygon(oceans.lon,oceans.lat,temp_lon,temp_lat));
    end
    clear temp temp_lon temp_lat
    clc
end

% Do the same for lakes
if ~isempty(input_lakes)
    for i = 1:length(lakes)
        if strcmp(v(end),')')
            fprintf('Lakes ID remaining: %d\n',length(lakes)-i);
        end
        % Check if the polygon contains NaN and remove it (required for octave)
        temp_lon = [lakes(i).Lon];
        temp_lat = [lakes(i).Lat];
        temp = find(isnan(temp_lon+temp_lat));
        if length(temp)>1% && ~strcmp(v(end),')')
            % First part = main polygon, other parts after NaN = remove
            % from inside of polygon => revert sign
            oceans.id = oceans.id + double(inpolygon(oceans.lon,oceans.lat,temp_lon(1:temp(1)-1),temp_lat(1:temp(1)-1)));
            for j = 2:length(temp)
                c = temp(j-1)+1;
                oceans.id = oceans.id - double(inpolygon(oceans.lon,oceans.lat,temp_lon(c:temp(j)-1),temp_lat(c:temp(j)-1)));
            end
            clear c j
        elseif length(temp) == 1 && temp~= length(temp)
            oceans.id = oceans.id + double(inpolygon(oceans.lon,oceans.lat,temp_lon(1:temp-1),temp_lat(1:temp-1)));
            oceans.id = oceans.id - double(inpolygon(oceans.lon,oceans.lat,temp_lon(temp+1:end),temp_lat(temp+1:end)));
        else
            oceans.id = oceans.id + double(inpolygon(oceans.lon,oceans.lat,temp_lon,temp_lat));
        end
        clc
    end
end
oceans.id(oceans.id<1) = 0;
oceans.id(oceans.id>1) = 1;

%% Save results
if show_grid == 1
    figure
    mesh(oceans.lon,oceans.lat,oceans.id);view(0,90);
end
% Reduce size
oceans.lon = oceans.lon(1,1:1:end);
oceans.lat = oceans.lat(1:1:end,1);
% Save
if ~strcmp(v(end),')')
    save(output_file,'oceans','-mat7-binary');
else
    save(output_file,'oceans');
end