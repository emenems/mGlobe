%% Compute mGlobe ATMO effect
% Using this script, user can call the mGlobe atmospheric effect function
% and compute the gravity effect for ERA Interim model.
clear
clc
% pkg load netcdf

%% Main settings
% Location
location = 'GGP_COORDINATES.txt'; 											% this variable points the the file with coordinates
output_file_prefix = fullfile('EXAMPLES','Scripts','Results');              % folder for results
select_site = 'CO';                                                         % set GGP abbreviation (e,g., 'PE'), or 'All' for all sites (see GGP_COORDINATES.txt file)
% Time
start_calc = datenum(2012,1,1,12,0,0);										% starting time
end_calc = datenum(2013,1,1,12,0,0);										% last epoch
step_calc = 2;                                                              % see  'switch step_calc' section (2 => 6 hours, 3=> 12 hours,...)
% Model
file_ref = 'I:\GlobalModel\ERAinterim\Invariant\ERA_INVARIANT_ALL.nc';      % orography file
file_temp = 'I:\GlobalModel\ERAinterim\Temperature\ERA_TEMP_24hStep_00_2000.nc'; % temperature file (for arbitrary year and hour)
file_humid = 'I:\GlobalModel\ERAinterim\SpecificHumidity\ERA_SHUMID_24hStep_00_2000.nc'; % spec. humidity file
file_height = 'I:\GlobalModel\ERAinterim\Geopotential\ERA_GEOPOT_24hStep_00_2000.nc'; % geopotential height file
file_sp = 'I:\GlobalModel\ERAinterim\Surface\ERA_SP_2T_2D_ALLh_06hStep_D_2000.nc'; % surface pressure, temperature, humidity file
% mGlobe
mglobe_folder = fullfile('..','..');                                        % mGlobe folder (all mGlobe functions are in this folder)

%% Additional settings
output_file_type = [0 1 0];                                                 % [xls, txt, tsf], 0 - off, 1 - on
subtract_average = 0;                                                       % subtract average values from all results, 0- off, 1 - on

currentFolder = pwd;                                                        % get current folder (do not change)
close_after_comutation = 0;                                                 % 0 = no, 1 = yes,i.e. close Matlab/Octave after computation

%% Load SG locations
if ~isempty(location) 														% if location variable points to a file
    fid = fopen(location);													% open the file
    temp = fgetl(fid);clear temp;                                           % read header
    sites = textscan(fid,'%s %f %f %f %s');									% all files
    fclose(fid);															% close the file
end

%% Compute
for i = 1:length(sites{1,2});											% loop for each site
    if strcmp(sites{1,1}(i),select_site) || strcmp(select_site,'All')   % compute only for selected sites
        % Name - time
        [cyear1,cmonth1] = datevec(start_calc);
        [cyear2,cmonth2] = datevec(end_calc);
        switch step_calc												% time resolution switch
            case 1
                time_resol = '3h';
            case 2
                time_resol = '6h';
            case 3 
                time_resol = '12h';
            case 4
                time_resol = '24h';
            case 5
                time_resol = '48h';
        end
        output_file = fullfile(output_file_prefix,sprintf('%s_%02d%04d_%02d%04d_%s_ERA_ATMO.txt',... % create output file name
            char(sites{1,1}(i)),cmonth1,cyear1,cmonth2,cyear2,time_resol));
        % Input
        Input = [sites{1,3}(i),sites{1,2}(i),sites{1,4}(i)];            % input coordinates
        % Compute
        cd(mglobe_folder)                                               % change folder to mGlobe
        mGlobe_calc_Atmo_ERA(Input,output_file,output_file_type,file_ref,file_temp,file_humid,file_height,file_sp,start_calc,end_calc,step_calc,subtract_average)
        clear Input output_file
        cd(currentFolder);                                              % change folder back
    end
end
if close_after_comutation == 1
    quit
end
