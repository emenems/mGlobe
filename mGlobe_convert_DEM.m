function mGlobe_convert_DEM(DEM_input,DEM_output,DEM_type)
%MGLOBE_CONVERT_DEM Function serves for the conversion of DEMs
%   Function is used for the conversion of DEM in geographic coordinate
%   system (latitude, longitude (DEG)). Input DEM can by transformed to 
%   matlab *.mat format used by the mGlobe. 
%   The mGlobe tool requires fix format:
%   dem.lon,dem.lat,dem.height.
% 
% Input:
%   DEM_input       ...     string with full path/name of the input DEM
%                           Example: 'VI_DEM_arc.ascii';
%   DEM_output      ...     string with full path/name of the output DEM
%                           Example: 'VI_DEM_arc.mat';
%   DEM_type        ...     number (1,3) for the identification of the
%                           input file format:  1 ... free (lon,lat,H)
%                                               2 ... free (lat,lon,H)
%                                               3 ... arc ASCII
%                                               4 ... grd 6 text
%                                               5 ... netCDF
%                           Example: 3
% Output:
%	dem.lon 	    ...		longitude (in input units)
%	dem.lat 	    ...		latitude (in input units)
%	dem.height 	    ...		longitude (in input units)
%	dem.input_file 	...		input file name
%	dem.units 	    ...		dem.height units
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

%% DEM conversion
set(findobj('Tag','text_status'),'String','Models: loading DEM ...');drawnow % status message
check_out = 0;                                                              % control value
try
    switch DEM_type
        case 1                                                              % Lon,Lat,Height txt format
            dem_original = load(DEM_input);                                 % load the input file
            dem.lat = linspace(min(min(dem_original(:,2))),max(max(dem_original(:,2))),length(unique(dem_original(:,2)))); % latitude vector
            dem.lon = linspace(min(min(dem_original(:,1))),max(max(dem_original(:,1))),length(unique(dem_original(:,1)))); % longitude vector
            [dem.lon,dem.lat] = meshgrid(dem.lon,dem.lat);                  % create new meshgrid lon/lat matrices
            dem.height = griddata(dem_original(:,1),dem_original(:,2),dem_original(:,3),dem.lon,dem.lat); % interpolate to new grid (linear interpolation)
			dem.input_file = DEM_input;                                     % store input model file name
			dem.units = 'm/see input model';                                % store model units
        case 2                                                              % Lon,Lat,Height txt format
            dem_original = load(DEM_input);
            dem.lon = linspace(min(min(dem_original(:,2))),max(max(dem_original(:,2))),length(unique(dem_original(:,2))));
            dem.lat = linspace(min(min(dem_original(:,1))),max(max(dem_original(:,1))),length(unique(dem_original(:,1))));
            [dem.lon,dem.lat] = meshgrid(dem.lon,dem.lat);
            dem.height = griddata(dem_original(:,2),dem_original(:,1),dem_original(:,3),dem.lon,dem.lat);
			dem.input_file = DEM_input;
			dem.units = 'm/see input model';
        case 3                                                              % arc ascii file format
            [dem.height,Ref] = arcgridread(DEM_input);                      % read file
            row = dem.height;col = dem.height;row(:,:) = 0;col(:,:) = 0;    % prepare variables
            for i = 1:size(dem.height,1)                                    % prepare matrix dimensions
                row(i,:) = ones(1,size(dem.height,2))*i;
                col(i,:) = 1:1:size(dem.height,2);
            end
            [dem.lon,dem.lat] = pix2map(Ref,row,col);                       % stack read height values
			dem.input_file = DEM_input;                                     % store input file name
			dem.units = 'm/see input model';                                % store units
            dem.lon = flipud(dem.lon);
            dem.lat = flipud(dem.lat);
            dem.height = flipud(dem.height);
        case 5                                                              % NetCDF format
            ncid = netcdf.open(DEM_input,'NC_NOWRITE');                     % open netcdf file
            [ndims,nvars] = netcdf.inq(ncid);                               % get variables
            for i = 1:nvars
                varname(i) = {netcdf.inqVar(ncid,i-1)};                     % sort variables/layers
            end
            % Longitude
            [selection,confirm] = listdlg('ListString',char(varname),'Name','Longitude (deg)','ListSize',[round(160*2),round(300*1.1)]); % let user to select the longitude layer
            if confirm == 1 && length(selection) == 1
                dem.lon = netcdf.getVar(ncid,selection-1,'double');
                try
                    scale_factor = netcdf.getAtt(ncid,selection-1,'scale_factor'); % check if scale factor does exist
                catch
                    scale_factor = 1;
                end
                try
                    add_offset = netcdf.getAtt(ncid,selection-1,'add_offset'); % check if add constant does exist
                catch
                    add_offset = 0;
                end
                dem.lon = dem.lon'*scale_factor + add_offset;               % scale and transpose the loaded height matrix
            end
            % Latitude
            [selection,confirm] = listdlg('ListString',char(varname),'Name','Latitude (deg)'); % same procedure as for longitude layer
            if confirm == 1 && length(selection) == 1
                dem.lat = netcdf.getVar(ncid,selection-1,'double');
                try
                    scale_factor = netcdf.getAtt(ncid,selection-1,'scale_factor');
                catch
                    scale_factor = 1;
                end
                try
                    add_offset = netcdf.getAtt(ncid,selection-1,'add_offset');
                catch
                    add_offset = 0;
                end
                dem.lat = dem.lat'*scale_factor + add_offset;
            if size(dem.lon,1) == 1 || size(dem.lon,2) == 1
                [dem.lon,dem.lat] = meshgrid(dem.lon,dem.lat);
            end
            end
            % Height
            [selection,confirm] = listdlg('ListString',char(varname),'Name','Height (m)'); % same procedure as for longitude layer
            if confirm == 1 && length(selection) == 1
                dem.height = netcdf.getVar(ncid,selection-1,'double');
                try
                    scale_factor = netcdf.getAtt(ncid,selection-1,'scale_factor');
                catch
                    scale_factor = 1;
                end
                try
                    add_offset = netcdf.getAtt(ncid,selection-1,'add_offset');
                catch
                    add_offset = 0;
                end
                dem.height = dem.height'*scale_factor + add_offset;
            end
            dem.input_file = DEM_input;
			dem.units = 'm/see input model';
            netcdf.close(ncid);
        case 4                                                              % grd 6 txt (grapher) file format
            fid = fopen(DEM_input);                                         % open file
            dsaa = fscanf(fid,'%s', 1);                                     % read header id (not used)
            dimen = fscanf(fid,'%d %d',2);                                  % get dimensions
            xlim = fscanf(fid,'%f %f',2);                                   % get x limit
            step = abs(diff(xlim)/(dimen(1)-1));                            % set resolution
            ylim = fscanf(fid,'%f %f',2);                                   % get x limit
            zlim = fscanf(fid,'%f %f',2);                                   % get z limit (not used)
            dem.height = fscanf(fid,'%f',dimen');
            x = xlim(1):step(1):xlim(2);                                    % create longitude vector
            if step ~= abs(diff(ylim)/(dimen(2)-1))
                   step(2) = abs(diff(ylim)/(dimen(2)-1));
                   y = ylim(1):step(2):ylim(2);                             % create latitude vector
            else
                y = ylim(1):step:ylim(2);
            end
            dem.height = dem.height';                                       % transpose height matrix
            if x(end) ~=xlim(2)
                x(end+1) = xlim(2);
            end
            if y(end) ~= ylim(2)
                y(end+1) = ylim(2);
            end
            dem.input_file = DEM_input;                                     % store input file name
			dem.units = 'm/see input model';                                % store units
            [dem.lon,dem.lat] = meshgrid(x,y);                              % create lon/lat matrices
    end
        set(findobj('Tag','text_status'),'String','Models: DEM converted...');drawnow % write status message
catch exception
    set(findobj('Tag','text_status'),'String','Models: Could not load the input DEM file (check format)'); drawnow
    check_out = 1;
    return
end

if check_out == 0;
    save(DEM_output,'dem');                                                 % save the transformed file
end 

end

