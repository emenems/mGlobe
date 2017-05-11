function mGlobe_convert_ECCO(start_calc,end_calc,time_resol,file,ghc_path,input_path,model_version)
%MGLOBE_CONVERT_ECCO Read and covert ECCO data
% Extract ocean bottom pressure data and save them to matlab format.  
% 
% ASSUMTPION:     
%   ECCO-JPL (month).. file name example:
%                      ECCO_kf080_2009335_2009365_AveRmvd_OBP.txt
%                  ... header lines == 14
%                  ... flagged values are omitted
%                  ... Ocean bottom pressure in hPa
%                  ... WARNING: only monthly JPL data supported
%                  ... ftp://podaac-ftp.jpl.nasa.gov/allData/tellus/L3/ecco_obp/ (unpacked!!)
% 
%   ECCO-JPL (12H) ... layer 0 == latitude (deg)
%                  ... layer 2 == longitude (deg)
%                  ... layer 5 == time
%                  ... layer 7 == Ocean Bottom Pressure Pot. Anomaly
%                  ... the conversion requires the same folder structure as:
%                      ftp://snowwhite.jpl.nasa.gov/data4/KalmanFilter/
% 
%   ECCO2          ... layer 1 == latitude (deg)
%                  ... layer 2 == longitude (deg)
%                  ... layer 0 == time
%                  ... layer 3 == Bottom Pressure Pot. Anomaly
%                  ... WARNING: ECCO2 in beta version
%                  ... ftp://ecco2.jpl.nasa.gov/data1/cube/cube92/lat_lon/quart_90S_90N/PHIBOT.nc/
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%					   Example: datenum(2012,1,1,12,0,0);
%   end_calc       ... finish time in matlab format (days)
%					   Example: datenum(2013,1,1,12,0,0);
%   time_resol     ... time resolution switcher (not in time units)
%						1 = 3 hours, 2 = 6 hours, 3 = 12 hours, 
%						4 = one day, 5 = two days
%						6 = month
%						Example: 4
%   file           ... input file name (string)
%						Example: 'OBPano_08_08.00001_02160_012.cdf'
%   ghc_path       ... path used for output (string)
%                      Example: fullfile('OBPM','ECCO1');
%   input_path     ... full path to one of the input files (string)
%                      Example: 'E:\Models\data4\KalmanFilter\kf080_2012\n10day_01_09\'
%   model_version  ... model ID: 1 = ECCO1 (JPL), 4 = ECCO2
%                      Example: 1
% 
% OUTPUT (automatically saved):
%   out_mat        ... structure array (several matrices) containing:
%   out_mat.lon    ... longitude (degrees)
%   out_mat.lat    ... latitude  (degrees)
%   out_mat.time   ... ECCO time (begins 2002/01/01 in days)
%   out_mat.obp    ... ECCO OBP (100*(newton/m2)=mbar~=cmH2O)
%   out_mat.input_file ... input file name
%   out_mat.units  ... out_mat.obp units
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

%% Convert data
switch model_version
    case 1
        if time_resol < 6
            for i = 1:length(time(:,7))
                try
                    if strcmp(input_path(end),'\')
                        fix_path = input_path(1:end-18);
                    else
                        fix_path = input_path(1:end-17);
                    end
                    year_folder = [fix_path,sprintf('%04d',time(i,1))];
                    day_of_year = fix(time(i,7)-datenum(time(i,1),1,1,0,0,0))+1;
                    if day_of_year<=90
                        % Get the folder name (depends on the total number
                        % of days => not constant)
                        temp = dir(fullfile(year_folder,'n10day_01_*'));
                        % Get file name (just like folder name)
                        loc_file = dir(fullfile(year_folder,temp(1).name,'OBPano*'));
                        % Combine folder and file names for reading
                        file_name = fullfile(year_folder,temp(1).name,loc_file(1).name);
                    elseif day_of_year<=180
                        temp = dir(fullfile(year_folder,'n10day_10_*'));
                        loc_file = dir(fullfile(year_folder,temp(1).name,'OBPano*'));
                        file_name = fullfile(year_folder,temp(1).name,loc_file(1).name);
                    elseif day_of_year<=270
                        temp = dir(fullfile(year_folder,'n10day_19_*'));
                        loc_file = dir(fullfile(year_folder,temp(1).name,'OBPano*'));
                        file_name = fullfile(year_folder,temp(1).name,loc_file(1).name);
                    else
                        temp = dir(fullfile(year_folder,'n10day_28_*'));
                        loc_file = dir(fullfile(year_folder,temp(1).name,'OBPano*'));
                        file_name = fullfile(year_folder,temp(1).name,loc_file(1).name);
                    end
					clear temp;
                    ncid = netcdf.open(file_name,'NC_NOWRITE');                     % open NetCDF file
                    out_mat.lat = netcdf.getVar(ncid,0,'double');                   % get latitude
                    out_mat.lon = netcdf.getVar(ncid,2,'double');                   % get longitude
                    time_all = netcdf.getVar(ncid,5,'double');                      % get time
                    r = find(datenum(1970,1,1,0,0,0)+time_all/24==time(i,7));
                    if ~isempty(r)
                        out_mat.time = time_all(r);
                        obp_all = netcdf.getVar(ncid,7,'double');
                        out_mat.obp = obp_all(:,:,r)';
                        out_mat.obp(out_mat.obp<=-10000000000.0000) = NaN;
                        out_mat.obp = out_mat.obp./100;
        %                 out_mat.lon(out_mat.lon>=180) = out_mat.lon(out_mat.lon>=180)-360;
                        [out_mat.lon,out_mat.lat] = meshgrid(out_mat.lon,out_mat.lat);  % crate new meshgrid matrices
                        out_mat.lon(out_mat.lon>=180) = out_mat.lon(out_mat.lon>=180) - 360; % transform coordinates to (-180,180) system  
                        ri = find(abs(diff(out_mat.lon(1,:)))==max(abs(diff(out_mat.lon(1,:)))));
                        out_mat.lon = horzcat(out_mat.lon(:,ri+1:end),out_mat.lon(:,1:ri)); % Connect matrices to remove discontinuity
                        out_mat.lat = horzcat(out_mat.lat(:,ri+1:end),out_mat.lat(:,1:ri));
                        out_mat.obp = horzcat(out_mat.obp(:,ri+1:end),out_mat.obp(:,1:ri));clear ri;
                        out_mat.input_file = file_name;                     % store input file name
                        out_mat.units = 'mbar';                             % store units
                        netcdf.close(ncid)                                  % close the netcdf file
                        save(fullfile(ghc_path,sprintf('ECCO1_12H_%04d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4))),'out_mat');
                        if size(time,1) > 2
                            out_message = sprintf('Models: converting ECCO model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
                        else
                            out_message = sprintf('Models: converting ECCO model ...'); % create status message
                        end
                    set(findobj('Tag','text_status'),'String',out_message); drawnow % write status message
                        clear out_mat ncit file_name fix_path year_folder day_of_year quarter_folder r
                    else
                        set(findobj('Tag','text_status'),'String',['Models: Could not find model values for : %04d/%02d/%02d %02d (in %s)',time(i,1),time(i,2),time(i,3),time(i,4),file_name]); drawnow
                        fprintf('Models: Could not find model values for : %04d/%02d/%02d %02d (in %s)\n',time(i,1),time(i,2),time(i,3),time(i,4),file_name);
                    end
                catch
                    set(findobj('Tag','text_status'),'String',['Models: Could not load: ',file_name]); drawnow
                    fprintf('Models: Could not load: %s\n',file_name);
                end
            end
        else                                                                  % switch to ECCO1 (JPL) 
            for i = 1:length(time(:,7))
                if time(i,2) == 1                                               % create input file name (using day of year)
                    file_name = sprintf('%s%04d%03d_%04d%03d%s',[input_path,file(1:11)],time(i,1),1,time(i,1),fix(time(i,2)*(365/12)),file(end-15:end));
                else
                    file_name = sprintf('%s%04d%03d_%04d%03d%s',[input_path,file(1:11)],time(i,1),round(time(i,2)*(365/12))-30,time(i,1),fix(time(i,2)*(365/12)),file(end-15:end));
                end
                try
                    data_original = dlmread(file_name,'',14,0);                 % load txt data
                    out_mat.lat = -89.5:1:89.5;                                 % crate output latitude vector
                    out_mat.lon = 0.5:1:359.5;
                    [out_mat.lon,out_mat.lat] = meshgrid(out_mat.lon,out_mat.lat); % create output meshgrid matrices
                    lon_lat = horzcat(out_mat.lon(:),out_mat.lat(:));
                    lon_lat_data(1:length(lon_lat)) = NaN;
                    for ir = 1:length(lon_lat)                                  % create ECCO-JPL OBP model
                        rid = find(lon_lat(ir,1) == data_original(:,1) & lon_lat(ir,2) == data_original(:,2)); % finding correct value for each grid cell (omitted values = NaN) 
                        if ~isempty(rid)
                            lon_lat_data(ir) = data_original(rid,3);
                        end
                        clear rid
                    end
                    out_mat.obp = reshape(lon_lat_data,size(out_mat.lon));      % reshape the created matrix
                    out_mat.lon(out_mat.lon>=180) = out_mat.lon(out_mat.lon>=180) - 360; % transform coordinates to (-180,180) system  
                    ri = find(abs(diff(out_mat.lon(1,:)))==max(abs(diff(out_mat.lon(1,:)))));
                    out_mat.lon = horzcat(out_mat.lon(:,ri+1:end),out_mat.lon(:,1:ri)); % Connect matrices to remove discontinuity
                    out_mat.lat = horzcat(out_mat.lat(:,ri+1:end),out_mat.lat(:,1:ri));
                    out_mat.obp = horzcat(out_mat.obp(:,ri+1:end),out_mat.obp(:,1:ri));clear ri;
                    out_mat.input_file = file_name;                             % store input file name
                    out_mat.units = 'mbar/see input model';                     % store output file name
                    save(fullfile(ghc_path,sprintf('ECCO1_M_%04d%02d.mat',time(i,1),time(i,2))),'out_mat'); % save model values
                    if size(time,1) > 2
                        out_message = sprintf('Models: converting ECCO model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
                    else
                        out_message = sprintf('Models: converting ECCO model ...'); % create status message
                    end
                    set(findobj('Tag','text_status'),'String',out_message); drawnow % write status message
                catch
                    set(findobj('Tag','text_status'),'String',['Models: Could not load: ',file_name]); drawnow
                    fprintf('Models: Could not load: %s\n',file_name);
                end
            end
        end
    case 4
        for i = 1:length(time(:,7))
            try
                file_name = sprintf('%s%04d%02d%02d.nc',[input_path,file(1:end-11)],time(i,1),time(i,2),time(i,3)); % create new file name (input)
                ncid = netcdf.open(file_name,'NC_NOWRITE');                     % open NetCDF file
                out_mat.lat = netcdf.getVar(ncid,1,'double');                   % get latitude
                out_mat.time = netcdf.getVar(ncid,0,'double');                  % get time
                out_mat.lon = netcdf.getVar(ncid,2,'double');                   % get longitude
                out_mat.obp = netcdf.getVar(ncid,3,'double');                   % for ECCO2  Bottom Pressure Pot. Anomaly 
                out_mat.obp(out_mat.obp < -10e+21) = NaN;                       % Remove flagged values
                out_mat.obp = out_mat.obp'*1027.5/100;                          % absolute bottom pressure in Pa is: Depth*rhonil*g + PHIBOT*rhonil (rhonil=1027.5 kg/m^3), /100 = to mbar
                [out_mat.lon,out_mat.lat] = meshgrid(out_mat.lon,out_mat.lat);  % crate new meshgrid matrices
                out_mat.lon(out_mat.lon>=180) = out_mat.lon(out_mat.lon>=180) - 360; % transform coordinates to (-180,180) system  
                ri = find(abs(diff(out_mat.lon(1,:)))==max(abs(diff(out_mat.lon(1,:)))));
                out_mat.lon = horzcat(out_mat.lon(:,ri+1:end),out_mat.lon(:,1:ri)); % Connect matrices to remove discontinuity
                out_mat.lat = horzcat(out_mat.lat(:,ri+1:end),out_mat.lat(:,1:ri));
                out_mat.obp = horzcat(out_mat.obp(:,ri+1:end),out_mat.obp(:,1:ri));clear ri;
                out_mat.input_file = file_name;                             % store input file name
                out_mat.units = 'mbar';                                     % store units
                if time_resol == 6                                          % save model values 
                    save(fullfile(ghc_path,sprintf('ECCO2_M_%04d%02d.mat',time(i,1),time(i,2))),'out_mat');
                else
                    save(fullfile(ghc_path,sprintf('ECCO2_D_%04d%02d%02d_12.mat',time(i,1),time(i,2),time(i,3))),'out_mat');
                end
                netcdf.close(ncid)                                              % close the netcdf file
                if size(time,1) > 2
                    out_message = sprintf('Models: converting ECCO model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
                else
                    out_message = sprintf('Models: converting ECCO model ...'); % create status message
                end
                set(findobj('Tag','text_status'),'String',out_message); drawnow % write status message
            catch exception
                set(findobj('Tag','text_status'),'String',['Models: Could not load: ',file_name]); drawnow
                fprintf('Models: Could not load: %s\n',file_name);
            end
        end
end

end

