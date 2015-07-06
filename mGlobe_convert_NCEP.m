function mGlobe_convert_NCEP(start_calc,end_calc,time_resol,input_path,ghc_path)
%MGLOBE_CONVERT_NCEP Read and convert NCEP Reanalysis-2 netcdf data
% Extract required NCEP (surface) model data stored in netcdf
% format. This function transforms all layers and thus require that all
% layers (soil moisture 0-10,10-200 and snow) are stored in the same 
% folder, i.e input_path.
% Data download:
% http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis2.gaussian.html
% 4-times Daily Values: Soil moisture (0-10cm)
%                       Soil moisture (10-200cm)
%                       Water equiv. of snow depth
% Monthly values are stored in one file for each compartment.
% Do not reaname the downloaded files!
% 
% ASSUMTPION:      ... layer 1 == longitude
%                      layer 0 == latitude
%                      layer 2 == time count
%                      layer 3(4) == snow (monthly data)
%                      layer 5(6) == soil moisture (monthly data)
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%					   Example: datenum(2012,1,1,12,0,0);
%   end_calc       ... finish time in matlab format (days)
%					   Example: datenum(2013,1,1,12,0,0);
%   time_resol     ... time resolution switcher: 1 == 3 hours, 2 == 6 hours,
%                      3 == 12 hours, 4 == 24 hours, 5 == 48 hours, 6 == month.
%                      example: 4
%   input_path     ... one of the NCEP NetCDF files
%                      Example: fullfile('E','models','NCEP','soilw.0-10cm.gauss.2012.nc');
%   ghc_path       ... path used for output
%                      Example: fullfile('GHM','NCEP');
% 
% OUTPUT (automatically saved):
%   out_mat        ... structure array (several matrices) containing:
%   out_mat.lon    ... longitude (degrees)
%   out_mat.lat    ... latitude  (degrees)
%   out_mat.time   ... NCEP time (begins on 01/01/1800 in hours)
%   out_mat.soilw1 ... soil moisture for layer 1 (0-10 cm, VolSM)
%   out_mat.soilw2 ... soil moisture for layer 2 (10-200cm, VolSM)
%   out_mat.weasd  ... snow water equivalent  (m)
%   out_mat.input_file ... input file name
%   out_mat.units  ... units of soil moisture and snow
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                27.11.2014
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
for i = 1:size(time,1);                                                     % create new file for each time epoch
    try
    % Soil moisture 0-10
    if time_resol == 6
        file = fullfile(input_path,'soilw.0-10cm.mon.mean.nc');
    else
        file = fullfile(input_path,sprintf('soilw.0-10cm.gauss.%04d.nc',time(i,1)));
    end
    ncid = netcdf.open(file,'NC_NOWRITE');                                  % open netcdf file
    longitude = netcdf.getVar(ncid,1,'double');                             % get latitude                      
    latitude = netcdf.getVar(ncid,0,'double');                              % get longitude
    time_count = netcdf.getVar(ncid,4,'double');                            % get NCEP time
    r = find(time_count == (time(i,7)-datenum(1800,1,1,0,0,0))*24);         % find corresponding time epoch
    if ~isempty(r)                                                          % continue only if such time epoch does exist
        if time_resol == 6
            channel = 6;
        else
            channel = 5;
        end
        temp_var = netcdf.getVar(ncid,channel,[0 0 0 r-1],[length(longitude) length(latitude) 1 1],'double'); % temporary variable
        temp_var(temp_var==32766) = NaN;                                    % remove negative values (no negative values are expected)
        scale_factor = netcdf.getAtt(ncid,channel,'scale_factor','double'); % get scaling factor
        add_offset = netcdf.getAtt(ncid,channel,'add_offset','double');     % get offset
        temp_var = temp_var'*scale_factor + add_offset;                     % scale
        temp_var(temp_var<0|isnan(temp_var)) = 0;
        out_mat.soilw1 = temp_var;
        [out_mat.lon,out_mat.lat] = meshgrid(longitude,latitude);           % meshgrid lon/lat matrices
        out_mat.time = time_count(r);                                       % store time
        out_mat.input_fileSoilw1 = file;                                    % store used input file
        out_mat.units = 'VolSM and kg/m2';                                  % store units
    end
    clear longitude latitude time_count r file temp_var channel
    netcdf.close(ncid);
    
    % Soil moisture 10-200
    if time_resol == 6
        file = fullfile(input_path,'soilw.10-200cm.mon.mean.nc');
    else
        file = fullfile(input_path,sprintf('soilw.10-200cm.gauss.%04d.nc',time(i,1)));
    end
    ncid = netcdf.open(file,'NC_NOWRITE');                                  % open netcdf file
    longitude = netcdf.getVar(ncid,1,'double');                             % get latitude                      
    latitude = netcdf.getVar(ncid,0,'double');                              % get longitude
    time_count = netcdf.getVar(ncid,4,'double');                            % get NCEP time
    r = find(time_count == (time(i,7)-datenum(1800,1,1,0,0,0))*24);         % find corresponding time epoch
    if ~isempty(r)                                                          % continue only if such time epoch does exist
        if time_resol == 6
            channel = 6;
        else
            channel = 5;
        end                                                                 % continue only if such time epoch does exist
        temp_var = netcdf.getVar(ncid,channel,[0 0 0 r-1],[length(longitude) length(latitude) 1 1],'double'); % temporary variable
        temp_var(temp_var==32766) = NaN;                                    % remove negative values (no negative values are expected)
        scale_factor = netcdf.getAtt(ncid,channel,'scale_factor','double'); % get scaling factor
        add_offset = netcdf.getAtt(ncid,channel,'add_offset','double');     % get offset
        temp_var = temp_var'*scale_factor + add_offset;                     % scale
        temp_var(temp_var<0|isnan(temp_var)) = 0;
        out_mat.soilw2 = temp_var;
        out_mat.input_fileSoilw2 = file;                                    % store used input file
    end
    clear longitude latitude time_count r file temp_var channel
    netcdf.close(ncid);
    
    % Snow water equivalent
    if time_resol == 6
        file = fullfile(input_path,'weasd.sfc.mon.mean.nc');
    else
        file = fullfile(input_path,sprintf('weasd.sfc.gauss.%04d.nc',time(i,1)));
    end
    ncid = netcdf.open(file,'NC_NOWRITE');                                  % open netcdf file
    longitude = netcdf.getVar(ncid,1,'double');                             % get latitude                      
    latitude = netcdf.getVar(ncid,0,'double');                              % get longitude
    time_count = netcdf.getVar(ncid,2,'double');                            % get NCEP time
    r = find(time_count == (time(i,7)-datenum(1800,1,1,0,0,0))*24);         % find corresponding time epoch
    if ~isempty(r)                                                          % continue only if such time epoch does exist
        if time_resol == 6
            channel = 4;
        else
            channel = 3;
        end  
        temp_var = netcdf.getVar(ncid,channel,[0 0 r-1],[length(longitude) length(latitude) 1],'double'); % temporary variable
        temp_var(temp_var==32766) = NaN;                                    % remove negative values (no negative values are expected)
        scale_factor = netcdf.getAtt(ncid,channel,'scale_factor','double'); % get scaling factor
        add_offset = netcdf.getAtt(ncid,channel,'add_offset','double');     % get offset
        temp_var = temp_var'*scale_factor + add_offset;                     % scale
        temp_var(temp_var<0|isnan(temp_var)) = 0;
        out_mat.weasd = temp_var;
        out_mat.input_fileweasd = file;                                    % store used input file
    end
    clear longitude latitude time_count r file temp_var channel
    netcdf.close(ncid);
    
    
    if size(time,1) > 2
        out_message = sprintf('Models: converting NCEP model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
    else
        out_message = sprintf('Models: converting NCEP model ...'); % create status message
    end
    set(findobj('Tag','text_status'),'String',out_message); drawnow         % write status message 
    if time_resol == 6                                                      % create new output file name (monthly or hourly data)
            nazov = fullfile(ghc_path,sprintf('NCEP_REAN2_M_%4d%02d.mat',time(i,1),time(i,2))); 
    else
            nazov = fullfile(ghc_path,sprintf('NCEP_REAN2_6H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4)));
    end
    save(nazov,'out_mat');                                                  % save create matrix
    
    catch
        if time_resol == 6
            out_message = sprintf('Models: Warning: could not convert data for NCEP_REAN2_M_%4d%02d.mat',time(i,1),time(i,2)); 
        else
            out_message = sprintf('Models: Warning: could not convert data for NCEP_REAN2_6H_%4d%02d%02d_%02d.mat!',time(i,1),time(i,2),time(i,3),time(i,4));
        end
        set(findobj('Tag','text_status'),'String',out_message); drawnow     % write warning message
        fprintf(out_message);fprintf('\n');
    end
end
set(findobj('Tag','text_status'),'String','Models: Conversion completed'); % final status message

end

