%% Compute mGlobe NTOL effect
% Using this script, user can call the mGlobe non-tidal ocean loading
% effect function and compute the gravity effect for supported OBPM
% models and selected sites.
clear
clc

%% Main settings
% Location
location = 'GGP_COORDINATES.txt'; 											% this variable points the the file with coordinates
output_file_prefix = 'F:\Documents\mGlobe\EXAMPLES\Ocean';                  % folder for results. Use relative or absolute Path
select_site = 'All';                                                        % set GGP abbreviation (e,g., 'PE'), or 'All' for all sites (see GGP_COORDINATES.txt file)
% Time
start_calc = datenum(2014,1,1,6,0,0);										% starting time
end_calc = datenum(2014,1,5,6,0,0);                                         % last epoch
step_calc = 3;                                                              % see  'switch step_calc' section (4 => daily, 3=> 12 hours,...)
% Model
select_model = 1;                                                           % use either number of the OBPM (1 = ECCO1) or 'All' for all models (see 'switch model_calc(j)' section)
ghc_treshold = 0.10;                                                        % in degree (spherical distance)
% mGlobe
mglobe_folder = fullfile('..','..');                                        % mGlobe folder (all mGlobe functions are in this folder). Relative or absolute path
mean_field = 2;                                                             % Subtract area mean (pressure) for each time epoch
pressure_time_series = {[],[],[]};                                          % non empty only if mean_field==3 (see mGlobe_calc_Ocean.m). Example: {'output_file.txt',{1},{9}};
% Model folder. Warning, the settings in 'mGlobe_PATH_Settings.txt' are
% in this case irrelevant!! The correct sub-folder name (e.g. NOAH025) will
% be appended automatically.
model_folder = 'F:\Documents\mGlobe\OBPM';

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
if ischar(select_model)														% use 'All' models if select_model = 'All'
    model_calc = [1 4 5 6 7];                                               % only daily models are supported
else
    model_calc = select_model;												% use only selected model
end

for j = 1:length(model_calc);												% loop for each model
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
                otherwise
                    time_resol = 'M';
            end
            % Name - model
            switch model_calc(j)                                            % switch between models
                case 1  
                   model_name = 'ECCO1';  
                   ghc_path = fullfile(model_folder,'ECCO1');
                case 4 
                   model_name = 'ECCO2';  
                   ghc_path = fullfile(model_folder,'ECCO2');
                case 5 
                   model_name = 'OMCT_oba'; 
                   ghc_path = fullfile(model_folder,'OMCT');
                case 6 
                   model_name = 'OMCT_ocn'; 
                   ghc_path = fullfile(model_folder,'OMCT');
                case 7 
                   model_name = 'OMCT_atm'; 
                   ghc_path = fullfile(model_folder,'OMCT');
            end
            output_file = fullfile(output_file_prefix,sprintf('%s_%02d%04d_%02d%04d_%s_%s_mf%1d_tre%03d.txt',... % create output file name
                char(sites{1,1}(i)),cmonth1,cyear1,cmonth2,cyear2,time_resol,model_name,mean_field,ghc_treshold*100));
            % Input
            Input = [sites{1,3}(i),sites{1,2}(i),sites{1,4}(i)];            % input coordinates

            % Compute
            cd(mglobe_folder);                                              % change folder to mGlobe
            mGlobe_calc_Ocean(Input,output_file,output_file_type,start_calc,end_calc,step_calc,ghc_treshold,ghc_path,model_calc(j),subtract_average,mean_field,pressure_time_series);
            clear Input output_file                                			% remove used variables
            cd(currentFolder);                                              % change folder back
        end
    end
end
if close_after_comutation == 1
    quit
end
