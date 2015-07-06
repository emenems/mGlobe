function mGlobe_download_GLDAS(start_calc,end_calc,model,time_resol,ghc_path)
%MGLOBE_DOWNLOAD_GLDAS Download GLDAS/MERRA data using OPeNDAP server architecture
% This function works only with Matlab r2012a or higher!
% Download required hydrological GLDAS/MERRA data, i.e. soil moisture and snow.
% Water store in vegetation is not downloaded.
% 
% ASSUMTPION:      
%   NOAH           ... layer 21-24 == soilm1-soilm4
%                  ... layer 17 == swe
%                  ... [longitude,latitude] = -179.875:0.25:179.875,
%                                             -59.875:0.25:89.875
%   CLM            ... layer 21-30 == soilm1-soilm10
%                  ... layer 17 == swe
%                  ... [longitude,latitude] = -179.5:1:179.5,-59.5:1:89.5
%   MOSAIC         ... layer 22-24 == soilm1-soilm3
%                  ... layer 17 == swe
%                  ... [longitude,latitude] = -179.5:1:179.5,-59.5:1:89.5
%   VIC            ... layer 15-17 == soilm1-soilm3
%                  ... layer 13 == swe
%                  ... [longitude,latitude] = -179.5:1:179.5,-59.5:1:89.5
%   MERRA          ... MAT1NXLND layer 34 = total water storage (TWLAND)
%                      size(lon,lat) = 541x360
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%                      Example: datenum([2012,1,1,12,0,0]);
%   end_calc       ... finish time in matlab format (days)
%                      Example: datenum([2013,1,1,12,0,0]);
%   time_resol     ... time resolution switcher (not in time units)
%					   1 == 3 hours, 2 == 6 hours,
%                      3 == 12 hours, 4 == 24 hours, 5 == 48 hours, 6 == month.
%                      Example: 4
%   model          ... GLDAS model identification: 
%         3        ... NOAH with 0.25 deg spatial resolution
%         4        ... NOAH with 1.0 deg spatial resolution
%         1        ... CLM model (1.0 deg spatial resolution)
%         2        ... MOSAIC model (1.0 deg spatial resolution)
%         5        ... VIC model (1.0 deg spatial resolution)
%         6        ... MERRA model (0.667x0.5 deg spatial resolution)
%                      Example: 1
%   ghc_path 	   ... output path (string)
%                      Example: fullfile('GHM','CLM');
% 
% OUTPUT (saved automatically):
%   out_mat        ... structure array (several matrices) containing:
%   out_mat.lon    ... longitude (degrees)
%   out_mat.lat    ... latitude  (degrees)
%   out_mat.time   ... GLDAS time (begins one year earlier then MATLAB 
%                    time, i.e. 01/01/0001 at 00:00 (matlab 01/01/0000 !!) 
%                    add 365 to use with datevec function
%   out_mat.soilmX ... soil moisture for layer X (kg/m2)
%   out_mat.swe    ... snow water equivalent  (kg/m2)
%   out_mat.twland ... Total Land Water Storage (kg/m2)
%   out_mat.units  ... swe and soilmX units
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
    if model==6
        time(:,5) = 30;
    end
    time(:,7) = datenum(time(:,1:6));
    clear days
end
progres_perc = linspace(1,99,size(time,1));
for i = 1:size(time,1);
    check_out = 0;
try
switch model
        %% N O A H   0.25 deg
    case 3
        switch time_resol
            case 6
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_NOAH025_M'; % MONTHLY
                ncid = netcdf.open(modis);                                  % open netCDF (crate connection)
                time_count = round((time(i,7) - (730181.0+365))/30.4167);   % determine the time (for given reference value)
                nazov = fullfile(ghc_path,sprintf('GLDAS_NOAH025_M_%4d%02d.mat',time(i,1),time(i,2))); % MONTHly data
            otherwise
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_NOAH025SUBP_3H'; % 3-Hourly data
                ncid = netcdf.open(modis);                                  % open netCDF (crate connection)
                time_count = round((time(i,7)- (730175.0+365))*24/3);       % determine the time (for given reference value)
                nazov = fullfile(ghc_path,sprintf('GLDAS_NOAH025SUBP_3H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4))); % Daily/3-hourly data   
        end
            lat_vec = netcdf.getVar(ncid,1,'double');
            lon_vec = netcdf.getVar(ncid,2,'double');
            lonsize = length(lon_vec);
            latsize = length(lat_vec);
            [out_mat.lon,out_mat.lat] = meshgrid(lon_vec,lat_vec);
            % Download data 
            soilm1 = netcdf.getVar(ncid,21,[0 0 time_count],[lonsize latsize 1],'double'); % get the soil moisture data
            soilm2 = netcdf.getVar(ncid,22,[0 0 time_count],[lonsize latsize 1],'double');
            soilm3 = netcdf.getVar(ncid,23,[0 0 time_count],[lonsize latsize 1],'double');
            soilm4 = netcdf.getVar(ncid,24,[0 0 time_count],[lonsize latsize 1],'double');
            swe = netcdf.getVar(ncid,17,[0 0 time_count],[lonsize latsize 1],'double'); % snow (water stored in vegetation neglected)
            out_time_gldas = netcdf.getVar(ncid,0,time_count);              % get the default GLDAS time
            % Create output matrix
            soilm1(soilm1>9.999e+19 | soilm1<0) = 0;out_mat.soilm1 = soilm1'; % set flagged values to zero + transpose the matrix (standard for NetCDF files)
            soilm2(soilm2>9.999e+19 | soilm2<0) = 0;out_mat.soilm2 = soilm2';
            soilm3(soilm3>9.999e+19 | soilm3<0) = 0;out_mat.soilm3 = soilm3';
            soilm4(soilm4>9.999e+19 | soilm4<0) = 0;out_mat.soilm4 = soilm4';
            swe(swe>9.999e+19 | swe<0) = 0;out_mat.swe = swe';
            out_mat.time = out_time_gldas;                                  % store the GLDAS time
            out_mat.units = 'mm';                                           % store the GLDAS units
            netcdf.close(ncid);                                             % close the netcdf file
            save(nazov,'out_mat');                                          % save the created matrix to *.mat file
            clear out_mat nazov out_time gldas soilm1 soilm2 soilm3 soilm4 swe
        %% N O A H   1 deg
    case 4
        switch time_resol
            case 6
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_NOAH10_M'; % server for monthly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - (722451.0+365))/30.4167);
                nazov = fullfile(ghc_path,sprintf('GLDAS_NOAH10_M_%4d%02d.mat',time(i,1),time(i,2))); 
            otherwise
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_NOAH10SUBP_3H'; % server for 3-hourly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - (722452.0+365))*24/3);      % OK
                nazov = fullfile(ghc_path,sprintf('GLDAS_NOAH10SUBP_3H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4)));
        end
            lat_vec = netcdf.getVar(ncid,1,'double');
            lon_vec = netcdf.getVar(ncid,2,'double');
            lonsize = length(lon_vec);
            latsize = length(lat_vec);
            [out_mat.lon,out_mat.lat] = meshgrid(lon_vec,lat_vec);
            % Download data 
            soilm1 = netcdf.getVar(ncid,21,[0 0 time_count],[lonsize latsize 1],'double');
            soilm2 = netcdf.getVar(ncid,22,[0 0 time_count],[lonsize latsize 1],'double');
            soilm3 = netcdf.getVar(ncid,23,[0 0 time_count],[lonsize latsize 1],'double');
            soilm4 = netcdf.getVar(ncid,24,[0 0 time_count],[lonsize latsize 1],'double');
            swe = netcdf.getVar(ncid,17,[0 0 time_count],[lonsize latsize 1],'double');
            out_time_gldas = netcdf.getVar(ncid,0,time_count);
            % Create output matrix
            soilm1(soilm1>9.999e+19 | soilm1<0) = 0;out_mat.soilm1 = soilm1';
            soilm2(soilm2>9.999e+19 | soilm2<0) = 0;out_mat.soilm2 = soilm2';
            soilm3(soilm3>9.999e+19 | soilm3<0) = 0;out_mat.soilm3 = soilm3';
            soilm4(soilm4>9.999e+19 | soilm4<0) = 0;out_mat.soilm4 = soilm4';
            swe(swe>9.999e+19 | swe<0) = 0;out_mat.swe = swe';
            out_mat.time = out_time_gldas;
            out_mat.units = 'mm';
            netcdf.close(ncid);
            save(nazov,'out_mat');
            clear out_mat nazov out_time gldas soilm1 soilm2 soilm3 soilm4 swe
        %% C L M
    case 1
        switch time_resol
            case 6
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_CLM10_M'; % server for monthly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - (722451.0+365))/30.4167);
                nazov = fullfile(ghc_path,sprintf('GLDAS_CLM10_M_%4d%02d.mat',time(i,1),time(i,2))); % MONTHly data
            otherwise
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_CLM10SUBP_3H'; % server 3-hourly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - (722452.0+365))*24/3);      % OK
                nazov = fullfile(ghc_path,sprintf('GLDAS_CLM10SUBP_3H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4))); % Daily data
                
        end
            lat_vec = netcdf.getVar(ncid,1,'double');
            lon_vec = netcdf.getVar(ncid,2,'double');
            lonsize = length(lon_vec);
            latsize = length(lat_vec);
            [out_mat.lon,out_mat.lat] = meshgrid(lon_vec,lat_vec);
            % Download data 
            soilm1 = netcdf.getVar(ncid,21,[0 0 time_count],[lonsize latsize 1],'double');
            soilm2 = netcdf.getVar(ncid,22,[0 0 time_count],[lonsize latsize 1],'double');
            soilm3 = netcdf.getVar(ncid,23,[0 0 time_count],[lonsize latsize 1],'double');
            soilm4 = netcdf.getVar(ncid,24,[0 0 time_count],[lonsize latsize 1],'double');
            soilm5 = netcdf.getVar(ncid,25,[0 0 time_count],[lonsize latsize 1],'double');
            soilm6 = netcdf.getVar(ncid,26,[0 0 time_count],[lonsize latsize 1],'double');
            soilm7 = netcdf.getVar(ncid,27,[0 0 time_count],[lonsize latsize 1],'double');
            soilm8 = netcdf.getVar(ncid,28,[0 0 time_count],[lonsize latsize 1],'double');
            soilm9 = netcdf.getVar(ncid,29,[0 0 time_count],[lonsize latsize 1],'double');
            soilm10 = netcdf.getVar(ncid,30,[0 0 time_count],[lonsize latsize 1],'double');
            swe = netcdf.getVar(ncid,17,[0 0 time_count],[lonsize latsize 1],'double');
            out_time_gldas = netcdf.getVar(ncid,0,time_count);
            % Create output matrix
            soilm1(soilm1>9.999e+19 | soilm1<0) = 0;out_mat.soilm1 = soilm1';
            soilm2(soilm2>9.999e+19 | soilm2<0) = 0;out_mat.soilm2 = soilm2';
            soilm3(soilm3>9.999e+19 | soilm3<0) = 0;out_mat.soilm3 = soilm3';
            soilm4(soilm4>9.999e+19 | soilm4<0) = 0;out_mat.soilm4 = soilm4';
            soilm5(soilm5>9.999e+19 | soilm5<0) = 0;out_mat.soilm5 = soilm5';
            soilm6(soilm6>9.999e+19 | soilm6<0) = 0;out_mat.soilm6 = soilm6';
            soilm7(soilm7>9.999e+19 | soilm7<0) = 0;out_mat.soilm7 = soilm7';
            soilm8(soilm8>9.999e+19 | soilm8<0) = 0;out_mat.soilm8 = soilm8';
            soilm9(soilm9>9.999e+19 | soilm9<0) = 0;out_mat.soilm9 = soilm9';
            soilm10(soilm10>9.999e+19 | soilm10<0) = 0;out_mat.soilm10 = soilm10';
            swe(swe>9.999e+19 | swe<0) = 0;out_mat.swe = swe';
            out_mat.time = out_time_gldas;
            out_mat.units = 'mm';
            netcdf.close(ncid);
            save(nazov,'out_mat');
            clear out_mat nazov out_time gldas soilm1 soilm2 soilm3 soilm4 soilm5 soilm 6 soilm7 soilm8 soilm10 swe
        %% M O S A I C
    case 2
        switch time_resol
            case 6
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_MOS10_M'; % server for monthly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - (722451.0+365))/30.4167);
                nazov = fullfile(ghc_path,sprintf('GLDAS_MOS10_M_%4d%02d.mat',time(i,1),time(i,2))); 
            otherwise
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_MOS10SUBP_3H'; % server for 3-hourly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - datenum(1979,1,2,0,0,0))*24/3); % OK 
                nazov = fullfile(ghc_path,sprintf('GLDAS_MOS10SUBP_3H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4)));
        end
            lat_vec = netcdf.getVar(ncid,1,'double');
            lon_vec = netcdf.getVar(ncid,2,'double');
            lonsize = length(lon_vec);
            latsize = length(lat_vec);
            [out_mat.lon,out_mat.lat] = meshgrid(lon_vec,lat_vec);
            % Download data 
            soilm1 = netcdf.getVar(ncid,22,[0 0 time_count],[lonsize latsize 1],'double');
            soilm2 = netcdf.getVar(ncid,23,[0 0 time_count],[lonsize latsize 1],'double');
            soilm3 = netcdf.getVar(ncid,24,[0 0 time_count],[lonsize latsize 1],'double');
            swe = netcdf.getVar(ncid,17,[0 0 time_count],[lonsize latsize 1],'double');
            out_time_gldas = netcdf.getVar(ncid,0,time_count);
            % Create output matrix
            soilm1(soilm1>9.999e+19 | soilm1<0) = 0;out_mat.soilm1 = soilm1';
            soilm2(soilm2>9.999e+19 | soilm2<0) = 0;out_mat.soilm2 = soilm2';
            soilm3(soilm3>9.999e+19 | soilm3<0) = 0;out_mat.soilm3 = soilm3';
            swe(swe>9.999e+19 | swe<0) = 0;out_mat.swe = swe';
            out_mat.time = out_time_gldas;
            out_mat.units = 'mm';
            netcdf.close(ncid);
            save(nazov,'out_mat');
            clear out_mat nazov out_time gldas soilm1 soilm2 soilm3 soilm4 swe
        %% V I C
    case 5
        switch time_resol
            case 6
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_VIC10_M'; % server for monthly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - (722451.0+365))/30.4167);
                nazov = fullfile(ghc_path,sprintf('GLDAS_VIC10_M_%4d%02d.mat',time(i,1),time(i,2))); 
            otherwise
                modis = 'http://hydro1.sci.gsfc.nasa.gov/dods/GLDAS_VIC10_3H'; % server for 3-hourly data
                ncid = netcdf.open(modis);
                time_count = round((time(i,7) - datenum(1979,1,1,0,0,0))*24/3-1); % OK 
                nazov = fullfile(ghc_path,sprintf('GLDAS_VIC10_3H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4)));
        end
            lat_vec = netcdf.getVar(ncid,1,'double');
            lon_vec = netcdf.getVar(ncid,2,'double');
            lonsize = length(lon_vec);
            latsize = length(lat_vec);
            [out_mat.lon,out_mat.lat] = meshgrid(lon_vec,lat_vec);
            % Download data 
            soilm1 = netcdf.getVar(ncid,15,[0 0 time_count],[lonsize latsize 1],'double');
            soilm2 = netcdf.getVar(ncid,16,[0 0 time_count],[lonsize latsize 1],'double');
            soilm3 = netcdf.getVar(ncid,17,[0 0 time_count],[lonsize latsize 1],'double');
            swe = netcdf.getVar(ncid,13,[0 0 time_count],[lonsize latsize 1],'double');
            out_time_gldas = netcdf.getVar(ncid,0,time_count);
            % Create output matrix
            soilm1(soilm1>9.999e+19 | soilm1<0) = 0;out_mat.soilm1 = soilm1';
            soilm2(soilm2>9.999e+19 | soilm2<0) = 0;out_mat.soilm2 = soilm2';
            soilm3(soilm3>9.999e+19 | soilm3<0) = 0;out_mat.soilm3 = soilm3';
            swe(swe>9.999e+19 | swe<0) = 0;out_mat.swe = swe';
            out_mat.time = out_time_gldas;
            out_mat.units = 'mm';
            netcdf.close(ncid);
            save(nazov,'out_mat');
            clear out_mat nazov out_time gldas soilm1 soilm2 soilm3 soilm4 swe
    case 6
        %% M E R R A
        switch time_resol
            case 6
                modis = 'http://goldsmr2.sci.gsfc.nasa.gov/dods/MATMNXLND/'; 
                ncid = netcdf.open(modis,'NC_NOWRITE');                     % open netCDF (crate connection)
                time_count = round((time(i,7) - datenum(1979,1,1,0,30,0))/30.4167); % OK, checked agains downloaded NetCDF
                nazov = fullfile(ghc_path,sprintf('MERRA_M_%4d%02d.mat',time(i,1),time(i,2)));
            otherwise
                modis = 'http://goldsmr2.sci.gsfc.nasa.gov/dods/MAT1NXLND/'; % Hourly
                ncid = netcdf.open(modis,'NC_NOWRITE');                                  % open netCDF (crate connection)
                time_count = round((time(i,7) - datenum(1979,1,1,0,30,0))*24); % OK, checked agains downloaded NetCDF, 00:30 01/01/1979 => ref_time = datenum(1979,1,1,0,30,0)-365 = 722451.020833333
                nazov = fullfile(ghc_path,sprintf('MERRA_1H_%4d%02d%02d_%02d.mat',time(i,1),time(i,2),time(i,3),time(i,4)));
        end
            lat_vec = netcdf.getVar(ncid,1,'double');
            lon_vec = netcdf.getVar(ncid,2,'double');
            lonsize = length(lon_vec);
            latsize = length(lat_vec);
            [out_mat.lon,out_mat.lat] = meshgrid(lon_vec,lat_vec);
            % Download
            out_mat.twland = netcdf.getVar(ncid,34,[0 0 time_count],[lonsize,latsize,1],'double');%TWLAND = MAXSOILWAT – CATDEF + SRFEXC + RZEXC + SNOMAS + CAPAC
            out_mat.twland(out_mat.twland>9.999e+19 | out_mat.twland<0) = 0;
            if sum(sum(out_mat.twland)) ~= 0                                % write only if downloaded model is not empty
                out_mat.twland = out_mat.twland';
                out_mat.time = time(i,7);
                out_mat.units = 'mm';
                netcdf.close(ncid);
                save(nazov,'out_mat');
            else
                fpritnf('Models: Empty MERRA model output: %s\n',nazov);
            end
            clear out_mat nazov
end
catch
    check_out = 1;
end
switch check_out
    case 0
        if size(time,1) > 2
            out_message = sprintf('Models: downloading GLDAS/MERRA model ... (%3.0f%%)',progres_perc(i)); % create message
        else
            out_message = sprintf('Models: downloading GLDAS/MERRA model ...'); % create message
        end
    case 1
        out_message= sprintf('Models: could not download data for: %04d%02d%02d_%02d',time(i,1),time(i,2),time(i,3),time(i,4));
        fprintf('%s\n',out_message);
end
set(findobj('Tag','text_status'),'String',out_message); drawnow  
clear nazov
end

