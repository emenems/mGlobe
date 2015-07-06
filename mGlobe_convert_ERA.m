function mGlobe_convert_ERA(start_calc,end_calc,time_resol,file,ghc_path)
%MGLOBE_CONVERT_ERA Read and convert ECMWF ERA Interim netcdf data
% Extract required ERA Interim (surface) model data stored in netcdf
% format. This function transforms all layers stored in NetCDF format to 
% mGlobe supported file format, i.e. structure array.
% Data download: http://apps.ecmwf.int/datasets/data/interim_full_daily/
% 
% ASSUMTPION:      ... layer 0 == longitude
%                      layer 1 == latitude
%                      layer 2 == time count
%                      layer 3...end == svwl1,svwl2,...sd (arbitrary order)
%                      no negative values for layers 3,4,...
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%   end_calc       ... finish time in matlab format (days)
%   time_resol     ... time resolution switcher (not in time units)
%   file           ... full path to the file with ERA data in netCDF format
%   ghc_path       ... path used for output
% 
% OUTPUT (automatically saved):
%   out_mat        ... structure array (several matrices) containing:
%   out_mat.lon    ... longitude (degrees)
%   out_mat.lat    ... latitude  (degrees)
%   out_mat.time   ... ERA time (begins on 01/01/1900 in hours)
%   out_mat.svwl   ... soil moisture for layer X (m3/m3)
%   out_mat.sd     ... snow water equivalent  (m)
%   out_mat.input_file ... input file name
%   out_mat.units  ... units of svwl and sd
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

%% Time setting
[year_s,month_s] = datevec(start_calc);                                     % transform matlab time to civil date
[year_e,month_e] = datevec(end_calc);
if time_resol == 6                                                          % create time for MONTHly data
    j = 1;
    for year = year_s:year_e
        if j == 1
            mz = month_s;
        else
            mz = 1;
        end
        if year == year_e
            mk = month_e;
        else
            mk = 12;
        end
        for m = mz:mk
            time(j,1) = year;
            time(j,2) = m;
            j = j + 1;
        end
    end
    time(:,3) = 1;
    time(:,7) = datenum(time(:,1),time(:,2),time(:,3));
else                                                                        % create time for other resolutions
    switch time_resol
        case 1
            time_resol_in_days = 3/24;
        case 2
            time_resol_in_days = 6/24;
        case 3 
            time_resol_in_days = 12/24;
        case 4
            time_resol_in_days = 1;
        case 5
            time_resol_in_days = 2;
    end
    days = start_calc:time_resol_in_days:end_calc;
    time = datevec(days);
    time(:,7) = days;
    clear days
end

%% Download data
try
    ncid = netcdf.open(file,'NC_NOWRITE');                                  % open netcdf file (ECMWF standard)
    latitude = netcdf.getVar(ncid,1,'double');                              % get latitude                      
    longitude = netcdf.getVar(ncid,0,'double');                             % get longitude
    time_count = netcdf.getVar(ncid,2,'double');                            % get ERA Interim time
    [numdims,numvars] = netcdf.inq(ncid);                                   % get ERA variables
catch
    set(findobj('Tag','text_status'),'String','Models: Load valid ERA Interim (netCDF) file'); drawnow 
    return
end

for i = 1:size(time,1);                                                     % create new file for each time epoch
    r = find(time_count == (time(i,7)-datenum(1900,1,1,0,0,0))*24);         % find corresponding time epoch
    if ~isempty(r)                                                          % continue only if such time epoch does exist
    for j = 3:numvars-1                                                     % transform all layers (not only svwlX and sd)!!
        name = netcdf.inqVar(ncid,j);                                       % get variable name
        temp_var = netcdf.getVar(ncid,j,[0 0 r-1],[length(longitude) length(latitude) 1],'double'); % temporary variable
        temp_var = temp_var';                                               % transpose
        scale_factor = netcdf.getAtt(ncid,j,'scale_factor');                % get scaling factor
        add_offset = netcdf.getAtt(ncid,j,'add_offset');                    % get offset
        temp_var = temp_var*scale_factor + add_offset;                      % transform to final units
        temp_var(temp_var<0) = 0;                                           % remove negative values (no negative values are expected)
        [out_mat.(name)] = temp_var;                                        % add new layer to the structure area field
    end
    [out_mat.lon,out_mat.lat] = meshgrid(longitude,latitude);               % meshgrid lon/lat matrices
    out_mat.time = time_count(r);                                           % store ERA interim time
	out_mat.input_file = file;                                              % store used input file
	out_mat.units = 'svwl = m3/m3; sd = m';                                 % store units
    if size(time,1) > 2
        out_message = sprintf('Models: converting ERA model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
    else
        out_message = sprintf('Models: converting ERA model ...'); % create status message
    end
    set(findobj('Tag','text_status'),'String',out_message); drawnow         % write status message 
    if time_resol == 6                                                      % create new output file name (monthly or hourly data)
            nazov = fullfile(ghc_path,sprintf('ERA_INTERIM_M_%4d%02d.mat',time(i,1),time(i,2))); 
    else
            nazov = fullfile(ghc_path,sprintf('ERA_INTERIM_6H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4)));
    end
    save(nazov,'out_mat');                                                  % save create matrix
    else
        out_mat = [];
        if time_resol == 6
            out_message = sprintf('Models: Warning: data for ERA_INTERIM_M_%4d%02d.mat not found!',time(i,1),time(i,2)); % create warning message
        else
            out_message = sprintf('Models: Warning: data for ERA_INTERIM_6H_%4d%02d%02d_%02d.mat not found!',time(i,1),time(i,2),time(i,3),time(i,4));
        end
        set(findobj('Tag','text_status'),'String',out_message); drawnow     % write warning message
        fprintf(out_message);fprintf('\n');
    end
    clear out_mat nazov out_time name temp_var scale_factor add_offset temp_var out_message
end
try
    netcdf.close(ncid);                                                     % close opened netcdf file
    set(findobj('Tag','text_status'),'String','Models: Conversion completed'); % final status message
end
end

