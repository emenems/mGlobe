%% Prepare ECCO data for mGlobe
% Use this script to download and convert ECCO1/2 (net)cdf files (no monthly
% data conversion). Make sure the connection is not blocked by firewall.
clear
clc

%% Settings
% Set time interval. Set carefully, this may overwrite existing files
time_start = [2016,12,31];
time_stop  = [2017,3,1];

% Set folder for data download ECCO1 & ECCO2
path_download1 = 'd:\GlobalModel\ECCO1\kfh_080\'; % ECCO1
path_download2 = 'd:\GlobalModel\ECCO2\PHIBOT\'; % ECCO2
% Set path with mGlobe OBPM (without 'ECCO1' & 'ECCO2' subfolder)
path_mglobe_obpm = 'F:\mikolaj\data\global_model\mglobe\OBPM\';
% Set path with mGlobe functions
path_mglobe = 'F:\mikolaj\code\libraries\mGlobe';

% Set what should be done
download_data = 0; % 0 = Off, 1 == ECCO1, 2 == ECCO2, 3 == ECCO1+ECCO2
convert_data = 1; % just like 'download_data'

% Close matlab/octave after downloading/conversion
close_after_comutation = 0; 

%% ECCO FTP server settings
% set FTP connection to ECCO1
ftp_ecco1 = 'snowwhite.jpl.nasa.gov';
remote_dir1 = 'NearRealTime/KalmanFilter';
username1 = 'anonymous';
password1 = 'anonymous';
ecco1_subfolder = 'kf080h_'; % kalman filter, version dependent. tested for *80h only

% set FTP connection to ECCO2
ftp_ecco2 = 'ecco2.jpl.nasa.gov';
remote_dir2 = 'data1/cube/cube92/lat_lon/quart_90S_90N/';
username2 = 'anonymous';
password2 = 'anonymous';
ecco2_subfolder = 'PHIBOT.nc';

%% Download ECCO1
script_folder = pwd;
if download_data == 1 || download_data == 3
    try
        % Establishe connection
        ftp_connect = ftp(ftp_ecco1,username1,password1);
        cd(ftp_connect,remote_dir1);
        % Run loopf for all years
        for i = time_start(1):time_stop(1)
            % switch to current year and create folder/subfolder names
            cur_folder = sprintf('%s%04d',ecco1_subfolder,i);
            cd(ftp_connect,cur_folder);
            % check if the download folder contains required structure
            if exist(fullfile(path_download1,cur_folder),'dir') ~= 7
                cd(path_download1)
                mkdir(cur_folder);
                cd(script_folder);
            end
            % Try downloading all 4 quarters
            ftp_dir = dir(ftp_connect);
            for j = 1:length(ftp_dir)
                if strcmp(ftp_dir(j).name(1:4),'n10d')
                    % Check if the file to be downloaded is within the
                    % requested time interval
                    cur_time(1) = datenum(i,1,1) + str2double(ftp_dir(j).name(end-1:end))*10;
                    if strcmp(ftp_dir(j).name(end-4:end-3),'01')
                        cur_time(2) = datenum(i,1,1) + str2double(ftp_dir(j).name(end-4:end-3));
                    else
                        cur_time(2) = datenum(i,1,1) + str2double(ftp_dir(j).name(end-4:end-3))*10;
                    end
                    if (cur_time(1)>=datenum(time_start) && cur_time(1)<=datenum(time_stop)) || ...
                       (cur_time(2)>=datenum(time_start) && cur_time(2)<=datenum(time_stop))
                        cd(ftp_connect,ftp_dir(j).name);
                        if exist(fullfile(path_download1,cur_folder,ftp_dir(j).name),'dir') ~= 7
                            cd(fullfile(path_download1,cur_folder));
                            mkdir(ftp_dir(j).name);
                            cd(script_folder);
                        end
                        % Go to download folder for downloading
                        cd(fullfile(path_download1,cur_folder,ftp_dir(j).name));
                        % download OBP anomaly only
                        disp(['Downloading ',cur_folder,'/',ftp_dir(j).name]);
                        mget(ftp_connect,'OBPano_*');
                        cd(ftp_connect,'..');
                        cd(fullfile(path_download1,cur_folder));
                        clc
                    end
                end
            end
            % Back to default ftp folder
            cd(ftp_connect,'..');
            cd(script_folder);
            clear cur_folder j ftp_dir
        end
        close(ftp_connect);
        clear i
    catch exception
        fprintf('An error occurred during ECCO1 downloading:\n%s\n',exception.message);
        try 
            cd(script_folder);
            close(ftp_connect);
        end
    end
end

%% Download ECCO2
if download_data == 2 || download_data == 3
    try
        % Establish connection
        ftp_connect2 = ftp(ftp_ecco2,username2,password2);
        cd(ftp_connect2,remote_dir2);
        cd(ftp_connect2,ecco2_subfolder);
        % Look up all available data
        ftp_dir = dir(ftp_connect2);
        cd(path_download2);
        for i = 1:length(ftp_dir)
            if length(ftp_dir(i).name)>26
                % Get the file data and compare it with set time interval
                file_date = str2double(ftp_dir(i).name(end-10:end-3));
                file_year = floor(file_date/10000);
                file_month= floor((file_date - file_year*10000)/100);
                file_day  = floor((file_date - file_year*10000 - file_month*100));
                if datenum(file_year,file_month,file_day)>=datenum(time_start) && ...
                   datenum(file_year,file_month,file_day)<=datenum(time_stop)
                    mget(ftp_connect2,ftp_dir(i).name);
                end
                clear file_date file_year file_month file_day
            end
        end
        cd(script_folder);
        close(ftp_connect2);
        clear i ftp_dir
    catch exception
        fprintf('An error occurred during ECCO2 downloading:\n%s\n',exception.message);
        try 
            cd(script_folder);
            close(ftp_connect2);
        end
    end
end

%% Convert data ECCO1
if convert_data == 1 || convert_data == 3
    try
        cd(path_download1);
        % First, get one file with ECCO1 (required for mGlobe_convert_ECCO)
        loc_dir = dir(fullfile(path_download1,...
                    sprintf('%s%04d',ecco1_subfolder,time_start(1)),'n10d*'));
        if length(loc_dir)>1
            % Look for one file
            for j = 1:length(loc_dir)
                loc_file = dir(fullfile(path_download1,...
                        sprintf('%s%04d',ecco1_subfolder,time_start(1)),...
                        loc_dir(j).name,'OBPano*'));
                if length(loc_file) >= 1
                    input_path = fullfile(path_download1,...
                        sprintf('%s%04d',ecco1_subfolder,time_start(1)),...
                        loc_dir(j).name);
                    break;
                end
            end
            if length(loc_file) >= 1
                cd(path_mglobe);
                mGlobe_convert_ECCO(datenum([time_start,6,0,0]),datenum([time_stop,18,0,0]),...
                    3,loc_file(1).name,fullfile(path_mglobe_obpm,'ECCO1'),...
                    input_path,1)
                cd(script_folder);
            else
                disp('Convert data ECCO1: set starting time to date with downloaded file');
            end
        else
            disp('Convert data ECCO1: set starting time to date with downloaded file');
        end
    catch exception
        cd(script_folder);
        fprintf('An error occurred during ECCO1 conversion:\n%s\n',exception.message);
    end
end

%% Convert data ECCO2
if convert_data == 2 || convert_data == 3
    try
        cd(path_download2);
        % First, get one file with ECCO2 (required for mGlobe_convert_ECCO)
        loc_dir = dir('PHIBOT.*.nc');
        if length(loc_dir)>1
            input_file = loc_dir(1).name;
            cd(path_mglobe);
            mGlobe_convert_ECCO(datenum(time_start),datenum(time_stop),...
                    4,input_file,fullfile(path_mglobe_obpm,'ECCO2'),...
                    path_download2,4)
        else
            disp('Convert data ECCO2: set starting time to date with downloaded file');
        end
    catch exception
        cd(script_folder);
        fprintf('An error occurred during ECCO2 conversion:\n%s\n',exception.message);
    end
end

if close_after_comutation == 1
    quit
end