%% Compute mGlobe CWS effect
% Using this script, user can call the mGlobe continental water storage
% effect function and compute the gravity effect for supported hydrological
% models and selected sites.
clear
clc

%% Main settings
% Location
location = 'GGP_COORDINATES.txt'; 											% this variable points the the file with coordinates
output_file_prefix = fullfile('EXAMPLES','Scripts','Results');              % folder for results. Path relative to mGlobe folder
select_site = 'All';                                                        % set GGP abbreviation (e,g., 'PE'), or 'All' for all sites (see GGP_COORDINATES.txt file)
% Time
start_calc = datenum(2012,1,1,12,0,0);										% starting time
end_calc = datenum(2014,1,1,12,0,0);										% last epoch
step_calc = 4;                                                              % see  'switch step_calc' section (4 => daily, 3=> 12 hours,...)
% Model
select_model = 4;                                                           % use either number of the GHM (4 = GLDAS/NOAH10) or 'All' for all models (see 'switch model_calc(j)' section)
mass_conserv = 2;                                                           % 1 - off, 2 - ocean layer (from mass excess), 3 - as given by input model
ghc_treshold = 0.10;                                                        % in degree (spherical distance)
% mGlobe
mglobe_folder = fullfile('..','..');                                        % mGlobe folder (all mGlobe functions are in this folder)
DEM_file = fullfile('EXAMPLES','Models','DEM_ETOPO2_example.mat');          % *.mat DEM (matlab array containing *.lon,*.lat,*.height layers), set to [] if not required. Path relative to mGlobe folder
INCLUDE_file = [];                                                           % [] for no inclusion polygon = all grid cells are used, otherwise, txt file for "inlude  only" area (longitude,latidude (deg))

%% Additional settings
exclude_calc = [1,1];                                                       % [greenland, antarctica], 0 - off, 1 - exclude
output_file_type = [0 1 0];                                                 % [xls, txt, tsf], 0 - off, 1 - on
model_layer = 1;                                                            % model layer, 1 - total water storage, 2 soil moisture 1,...
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
    model_calc = [1 2 3 5 6 7 10];                                          % only daily models are supported
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
                   model_name = 'CLM';  
                   ghc_path = 'GHM\CLM';
                case 2  
                   model_name = 'MOS'; 
                   ghc_path = 'GHM\MOS'; 
                case 3  
                   model_name = 'NOAH025';  
                   ghc_path = 'GHM\NOAH025';
                case 4 
                   model_name = 'NOAH10';  
                   ghc_path = 'GHM\NOAH10';
                case 5 
                   model_name = 'VIC'; 
                   ghc_path = 'GHM\VIC';
                case 6 
                   model_name = 'ERA'; 
                   ghc_path = 'GHM\ERA';
                case 7 
                   model_name = 'MERRA'; 
                   ghc_path = 'GHM\MERRA';
                case 8 
                   model_name = 'Other'; 
                   ghc_path = 'GHM\OTHER';
                case 9 
                   model_name = 'GRACE'; 
                   ghc_path = 'GRACE\LAND';
                case 10 
                   model_name = 'NCEP'; 
                   ghc_path = 'GHM\NCEP';
            end
            if isempty(INCLUDE_file)										% inclusion polygon
                inc = 'All';
            else
                inc = 'Polygon';
            end
            output_file = fullfile(output_file_prefix,sprintf('%s_%02d%04d_%02d%04d_%s_%s_total_gre%1d_ant%1d_inc%s_mc%1d_dem1_tre%03d.txt',... % create output file name
                char(sites{1,1}(i)),cmonth1,cyear1,cmonth2,cyear2,time_resol,model_name,exclude_calc(1),exclude_calc(2),inc,mass_conserv-1,ghc_treshold*100));
            % Input
            Input = [sites{1,3}(i),sites{1,2}(i),sites{1,4}(i)];            % input coordinates

            % change folder to mGlobe
            cd(mglobe_folder)                                               
			% load DEM to check if it contains point of computation
            if ~isempty(DEM_file)	
                if ~exist('dem','var')
                    dem = importdata(DEM_file);
                end
                Hint = interp2(dem.lon,dem.lat,dem.height,Input(2),Input(1)); % check DEM height
                if isnan(Hint)
                    curr_DEM = [];                                          % set to [], if out of DEM area
                    Input(3) = 0;
                else
                    curr_DEM = DEM_file;                                    % otherwise use given DEM
                end
            else
				curr_DEM = [];
                Input(3) = 0;                                               % set height to zero if DEM not loaded
            end
			% Compute
            mGlobe_calc_Hydro(Input,output_file,output_file_type,curr_DEM,start_calc,end_calc,step_calc,exclude_calc,model_calc(j),model_layer,mass_conserv,ghc_treshold,ghc_path,subtract_average,INCLUDE_file)
            clear Input output_file curr_DEM                                        % remove used variables
            cd(currentFolder);                                              % change folder back
        end
    end
end
if close_after_comutation == 1
    quit
end
