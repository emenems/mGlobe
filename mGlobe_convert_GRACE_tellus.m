function mGlobe_convert_GRACE_tellus(start_calc,end_calc,file_path,output_name,ghc_path,ocean_land,scale_file)
%MGLOBE_CONVERT_GRACE Read and convert GRACE TELLUS netcdf data 
% Extract required GRACE (TELLUS) data stored in netcdf format
% Data download: ftp://podaac-ftp.jpl.nasa.gov/allData/tellus/L3/ocean_mass/RL05/netcdf/
%                ftp://podaac-ftp.jpl.nasa.gov/allData/tellus/L3/land_mass/RL05/netcdf/
% 
% % ASSUMTPION:  ... layer 0 == longitude
%                    layer 1 == latitude
%                    layer 2 == time
%                    layer 3 or 4 equivalent water thickness (cm)
%                    flagged value == 32767 or -9999         
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%					   Example: datenum(2012,1,1,12,0,0);
%   end_calc       ... finish time in matlab format (days)
%					   Example: datenum(2013,1,1,12,0,0);
%   file_path      ... full file name with GRACE data in netCDF format (string)
%					   Example: fullfile('E','Models','GRACE','GRCTellus.GFZ.200208_201501.OCN.RL05.DSTvDPC1409.nc');
%   file_name      ... file name with GRACE data in netCDF format (used for
%                       output name)
%                      Example: 'GRC_GFZ_RL05_CONv1409s'
%   ghc_path       ... path used for output
%                      Example: fullfile('GRACE','LAND');
%   ocean_land     ... Ocean or Land grid (2 => Ocean, 1 => Land)
%                      Example: 1
%   scale_file     ... path to scale file (used for multiplication with 
%                       given GRACE grid to minimize the effect of 
%                       filtering). [] if ocean_land == 2;
%						Example: fullfile('E','Models','GRACE','CLM4.SCALE_FACTOR.DS.G300KM.RL05.DSTvSCS1409.nc');
% 
% OUTPUT (automatically saved):
%   grace         ... structure array (several matrices) containing:
%   grace.lon     ... longitude (degrees)
%   grace.lat     ... latitude  (degrees)
%   grace.time    ... [start,mid,end]
%   grace.total   ... equivalent water thickness in cm (flagged values=NaN)
%   grace.units   ... grace.total units
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

%% LOAD data
try
    set(findobj('Tag','text_status'),'String','Models: Loading GRACE...'); drawnow  % write status message
    if ~isempty(scale_file)                                                 % load scale matrix-if given
        scale_ncid = netcdf.open(scale_file,'NC_NOWRITE');
        [scale_ndims,scale_nvars] = netcdf.inq(scale_ncid);
        for i = 1:scale_nvars                                               % let user pick the scale matrix
            varname(i) = {netcdf.inqVar(scale_ncid,i-1)};
        end
        [scale_selection,scale_confirm] = listdlg('ListString',char(varname),'Name','Pick SCALE matrix','ListSize',[round(160*2),round(300*1.1)]); % pick scale layer
        if scale_confirm == 1 && length(scale_selection) == 1
            temp_scale = netcdf.getVar(scale_ncid,scale_selection-1,'double'); % load scale matrix (transpose later)
            temp_scale(temp_scale == 32767 | temp_scale == -9999) = NaN;    % remove flagged values (assumption flag = 32767|-9999)
        else
            temp_scale = 1;                                                 % one if user does not pick the scale matrix
        end
    else
        scale_file = 'No scaling';
        temp_scale = 1;
    end
    ncid = netcdf.open(file_path,'NC_NOWRITE');                             % open GRACE file
    latitude = netcdf.getVar(ncid,1,'double');                              % read latitude
    longitude = netcdf.getVar(ncid,0,'double');                             % read longitude
    time = netcdf.getVar(ncid,2,'double');                                  % read time info (in days)
    try_var = netcdf.getVar(ncid,3,'double');  
    if size(try_var,1) > 3
        time_bounds = netcdf.getVar(ncid,4,'double');                       % bounding time
        lwe_thickness = netcdf.getVar(ncid,3,'double');                     % water layer in cm
    else
        time_bounds = netcdf.getVar(ncid,3,'double');                       % bounding time
        lwe_thickness = netcdf.getVar(ncid,4,'double');                     % water layer in cm
    end
    time_civil = datenum(2002,1,1)+time;                                    % transform to civil date
    rid = find(time_civil>=start_calc & time_civil<=end_calc);              % use only data specified by user
    if ~isempty(rid)
        set(findobj('Tag','text_status'),'String','Models: Converting GRACE...'); drawnow
        for i = 1:length(rid)
            temp_var = lwe_thickness(:,:,rid(i));
            switch ocean_land                                               % remove flagged value for Ocean or Land
                case 1 
                    temp_var(temp_var==32767) = NaN;
                case 2
                    temp_var(temp_var==-9999) = NaN;
            end
            grace.total = temp_scale'.*temp_var';                           % transpose matrix so it corresponds to grace.lon, grace.lat and multiply with scale matrix
            grace.time = [time_bounds(1,rid(i)),time(rid(i)),time_bounds(2,rid(i))];       % output time
            [grace.lon,grace.lat] = meshgrid(longitude,latitude);           % longitude and latitude matrix (in 0 to 360 system)
            grace.lon(grace.lon>=180) = grace.lon(grace.lon>=180) - 360;    % transform coordinates to (-180,180) system  
            ri = find(abs(diff(grace.lon(1,:)))==max(abs(diff(grace.lon(1,:)))));
            grace.lon = horzcat(grace.lon(:,ri+1:end),grace.lon(:,1:ri));   % Connect matrices to remove discontinuity
            grace.lat = horzcat(grace.lat(:,ri+1:end),grace.lat(:,1:ri));
            grace.total = horzcat(grace.total(:,ri+1:end),grace.total(:,1:ri));
            [year_temp,month_temp,day_temp,hour_temp,minute_temp,sec_temp] = datevec(time_civil(rid(i)));
            grace.input_file = file_path;
            grace.scale_input = scale_file;
            grace.units = 'cm';
            switch ocean_land
                case 1
                    save(fullfile(ghc_path,'LAND',sprintf('%s_%04d%02d%02d_%02d.mat',output_name,year_temp,month_temp,day_temp,hour_temp)),'grace');
                case 2
                    save(fullfile(ghc_path,'OCEAN',sprintf('%s_%04d%02d%02d_%02d.mat',output_name,year_temp,month_temp,day_temp,hour_temp)),'grace');
            end
            clear ri temp_var grace year_temp month_temp day_temp hour_temp minute_temp sec_temp
        end
    end
    netcdf.close(ncid);
catch
    set(findobj('Tag','text_status'),'String','Models: GRACE: Load valid netCDF file (with lon,lat,time,time_bounds and lwe_thickness Layers)'); drawnow
end

end
