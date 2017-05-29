%% Use this script to 
% 1. For Windows users: Update list for wget/DownThemAll data downloading.
%    The aim is to remove all files that are not of interest. The NASA SSW 
%    portal generates file containing all available hours. If the user, 
%    however, wishes to compute GHE only for daily data than the NASA SSW 
%	 URL list needs to  be modified so only files of interest are 
%	 downloaded (e.g., at 12:00 each day).
%	 Does not work with GLDASv2.1!
%	 In Linux/Unix use 'grep' for list reduction 
%	 Example for 12:00 only: grep '1200' input_list.txt > output.txt)
%
% 2. Convert the downloaded files:
%    Use this part to convert the downloaded data without using GUI.
%
clear 
close all
clc

%% Main Settings
% Set which part of the script should be carried out, i.e. update URL list
% (1) or convert downloaded files (2)
task_switch = 2;

% Close Matlab/Octave after completing the script?
close_matlab = 0; % 1 = yes, 0 = no;

%% Settings: Task 1
% Set the input file with URLs downloaded from: http://disc.gsfc.nasa.gov/SSW
input_urls = 'F:\Downloads\SSW_download_2016-08-18T15_51_40_35173_fdLkMV3Q.inp.txt'; 
% Set the output file name with modified URL list
output_urls = 'F:\Downloads\SSW_download_2016-08-18T15_51_40_35173_fdLkMV3Q.html'; % use html to allow opening in Mozilla.
% Set which time stamp should be kept (all other will be removed from the
% list). Use this script only for hourly not monthly data. Does not use with
% MERRA as all hours are in one file.
use_hour = {'0000','1200'};

%% Settings: Task 2
% Set the folder containing mGlobe
mglobe_folder = fullfile('..','..');
% Starting date
date_start = [2015 08 01 12 00 00]; % yyyy mm dd HH MM SS
% Last date
date_stop = [2015 08 08 12 00 00]; % yyyy mm dd HH MM SS
% Time resolution switch
step_calc = 4; % 2 => 6 hours, 3=> 12 hours, ...
% Set folder name containing all NetCDF files. Use the 
% model_ver for GLDAS NOAH to switch between 1 or 2(.1)!
input_path = 'd:\GlobalModel\MERRA2\';
model_ver = 1; % 1 or 2. will be used only if GLDAS NOAH model is converted
% Set ONE file name (just to generate correct input file name for all
% required time steps)
input_file = 'MERRA2_400.tavg1_2d_lnd_Nx.20161218.nc4.nc';
% Output folder. Warning, the settings in 'mGlobe_PATH_Settings.txt' are
% in this case irrelevant!! The correct sub-folder name (e.g. NOAH025) will
% be appended automatically.
output_folder = 'f:\mikolaj\code\libraries\mGlobe\GHM\';
% NO progress-bar will be shown!

%% Main code
% Switch between tasks (not possible to carry out both at once)
switch task_switch
    %% TAKS 1
    case 1
        % Open the input file for reading
        fid = fopen(input_urls,'r');
        % Open the output file for writing
        fid_out = fopen(output_urls,'w');
        % Read first line of the input file
        row = fgetl(fid);
        
        % Get the model version to use correct index of time
        model_index = strfind(row,'CLM10SUBP_3H');
        if ~isempty(model_index)
            time_index = model_index(2) + [22:25];
        end
        model_index = strfind(row,'MOS10SUBP_3H');
        if ~isempty(model_index)
            time_index = model_index(2) + [22:25];
        end
        model_index = strfind(row,'NOAH025SUBP_3H');
        if ~isempty(model_index)
            time_index = model_index(2) + [24:27];
        end
        model_index = strfind(row,'NOAH10SUBP_3H');
        if ~isempty(model_index)
            time_index = model_index(2) + [23:26];
        end
        model_index = strfind(row,'VIC10_3H');
        if ~isempty(model_index)
            time_index = model_index(2) + [18:21];
        end
        
        if exist('time_index','var')
            %% Run loop checking each line
            while ischar(row)
                % Run loop for all hours that should not be deleted
                for i = 1:length(use_hour)
                    if strcmp(row(time_index),use_hour{i})
                        fprintf(fid_out,'%s\n',row);
                    end
                end
                row = fgetl(fid);
            end
        end
        fclose(fid);
        fclose(fid_out);
        
    %% TASK 2
    case 2
        % Get current folder to switch back after conversion
        currentFolder = pwd;
        % Switch folders to mGlobe
        cd(mglobe_folder) % change folder to mGlobe
        % Get input model number (switch) and set model output folder
        switch input_file(1:11)
            case 'GLDAS_CLM10'
                model = 1;
                ghc_path = fullfile(output_folder,'CLM');
            case 'GLDAS_MOS10'
                model = 2;
                ghc_path = fullfile(output_folder,'MOS'); 
            case 'GLDAS_NOAH0'
				if model_ver == 1
					model = 3;
					ghc_path = fullfile(output_folder,'NOAH025');
				elseif model_ver == 2
					model = 8;
					ghc_path = fullfile(output_folder,'NOAH025v21');
				end
            case 'GLDAS_NOAH1'
                model = 4;
                ghc_path = fullfile(output_folder,'NOAH10');
            case 'GLDAS_VIC10'
                model = 5;
                ghc_path = fullfile(output_folder,'VIC');
            case 'MERRA300.pr'
                model = 6;
                ghc_path = fullfile(output_folder,'MERRA');
            case 'MERRA2_200.'
                model = 7;
                ghc_path = fullfile(output_folder,'MERRA2');
            case 'MERRA2_300.'
                model = 7;
                ghc_path = fullfile(output_folder,'MERRA2');
            case 'MERRA2_400.'
                model = 7;
                ghc_path = fullfile(output_folder,'MERRA2');
            otherwise
                model = 0;
                ghc_path = fullfile(output_folder,'OTHER');
        end
        % Call conversion script
        mGlobe_convert_GLDAS(datenum(date_start),datenum(date_stop),...
							model,step_calc,ghc_path,input_path,input_file);
        % change folder back
        cd(currentFolder);           
end

if close_matlab == 1
    quit
end