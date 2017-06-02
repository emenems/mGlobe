function mGlobe_convert_OMCT(start_calc,end_calc,time_resol,file,ghc_path,input_path,model_version)
%MGLOBE_CONVERT_OMCT Read and covert OMCT spherical harmonics 
% Extract ocean bottom pressure data and save them to matlab format.  
% 
% ASSUMTPION:     
%   The AOD1B data has been downloaded (and unpacked) from:
%           http://isdc.gfz-potsdam.de    OR
%           ftp://podaac-ftp.jpl.nasa.gov/allData/grace/L1B/GFZ/AOD1B/
%           RL05 max degree = 100 + 6H time resolution max
%           RL05 max degree = 180 + 3H time resolution max
%           + see section Constants
%   The computation is based on the paper by Wahr et al., 1998. Time 
%   variability of the Earth's gravity field: Hydrological and oceanic 
%   effects and their possible detection using GRACE. Journal of 
%   Geophysical Research, vol. 103, no. B12, p30205-30229. 
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%                      Example: datenum([2012,1,1,12,0,0]);
%   end_calc       ... finish time in matlab format (days)
%                      Example: datenum([2013,1,1,12,0,0]);
%   time_resol     ... time resolution switcher (not in time units): 1 == 3 hours, 2 == 6 hours,
%                      3 == 12 hours, 4 == 24 hours, 5 == 48 hours.
%                      Example: 4
%   file           ... one of the OMCT input files (string)
%                      Example: 'AOD1B_2012-01_05.asc';
%   ghc_path       ... path used for output (string)
%                      Example: fullfile('OBPM','OMCT');
%   input_path     ... full path to one of the input files
%                      Example: fullfile('E','models','OMCT');
%   model_version  ... Product switch: 
%                      RL05: 5 = oba, 6 = ocn, 7 = atm
%                      RL06: 8 = oba, 9 = ocn, 10 = atm
%                      Example: 5
% 
% OUTPUT (automatically saved):
%   omct           ... structure array (several matrices) containing:
%   omct.lon       ... longitude (degrees)
%   omct.lat       ... latitude  (degrees)
%   omct.time      ... matlab time
%   omct.obp       ... OMCT water column in mm
%   omct.input_file ... input file name
%   omct.units     ... omct.obp units
%   omct.spherHarmDegree ... max. degree of spherical harmonics
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%   
%% Time setting
% transform matlab time to civil date 
[year_s,month_s] = datevec(start_calc); 
[year_e,month_e] = datevec(end_calc);
% create time for MONTHly data
if time_resol == 6 
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
else % create time for other resolutions
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

%% Constants
% spatial resolution of the new grid (deg)
resolution = 1.0;                                                           
% maximal degree of development + ellipsoid definition
if model_version < 8 % RL05
    max_deg = 100;
    a = 0.63781364600000E+07;% major axis (m)
else % RL06
    max_deg = 180;
    a = 0.63781366000000E+07;% major axis (m)
end
ro_ave = 5517; % mean density of the Earth
ro_water = 1000;  % mean water density
% longitude and latitude grid (deg) of the output file
[omct.lon,omct.lat] = meshgrid(0+resolution/2:resolution:360-resolution/2,... 
                           -90+resolution/2:resolution:90-resolution/2);
omct.units = 'mm';
% Constant term used for spherical synthesis
term0 = ((a*ro_ave)/(3*ro_water));
% Switch between OMCT releases and layers to crate output file name
% (date+time will be appended at the end)
switch model_version
    case 5
        aod1b_type = 'oba';
        output_prefix = 'OMCT_oba_6H_';
    case 6
        aod1b_type = 'ocn';
        output_prefix = 'OMCT_ocn_6H_';
    case 7
        aod1b_type = 'atm';
        output_prefix = 'OMCT_atm_6H_';
    case 8
        aod1b_type = 'oba';
        output_prefix = 'OMCT6_oba_3H_';
    case 9
        aod1b_type = 'ocn';
        output_prefix = 'OMCT6_ocn_3H_';
    case 10
        aod1b_type = 'atm';
        output_prefix = 'OMCT6_atm_3H_';
end
% file prefix and suffix of the input (ascii) file
file_prefix = 'AOD1B_';                                                     
file_suffix = file(end-8:end);
% load love numbers input file (part of mGlobe toolbox)
file_love = 'mGlobe_DATA_Load_degree_k.txt'; 
% auxiliary grid for Continents/Oceans (part of mGlobe toolbox)
file_oceans = 'mGlobe_DATA_OceanGrid.mat';

%% Prepare variables
% unique co-latitude (to reduce computation time, compute only for unique
% latitude)
lat_uniq = 90 - unique(omct.lat(:));
% compute cosine for faster computation (constant for all time epochs)
cosinus = cos(lat_uniq*pi/180);     
% all longitude grid cells = to use unique latitude and all possible
% longitude combination
lon_vec = omct.lon(:)'*pi/180;
% Prepare 'unique' indices for reshaping the compute vector (lon/lat
% combination back to full matrix)
in_uniq(1:length(omct.lat(:))) = 0;
for i = 1:length(omct.lat(:))
    in_uniq(i) = find(90-omct.lat(i)==lat_uniq);
end
% remove used variables
clear i lat_uniq

%% Load required files (aux) files and prepare grid
% load love numbers
if ~isempty(file_love)                                                      
    love.data = load(file_love);
    % degree in first column
    love.degree = love.data(:,1);
    % k in second column !!!!
    love.k = love.data(:,2);
end
% load auxiliary grid used for continents/oceans
if ~isempty(file_oceans)
    oceans = importdata(file_oceans);
    if size(oceans.lon,1) == 1 || size(oceans.lon,2) == 1
       [oceans.lon,oceans.lat] = meshgrid(oceans.lon,oceans.lat); 
    end
    new = omct;
    % transform coordinates/longitude to (-180,180) system  
    new.lon(new.lon>=180) = new.lon(new.lon>=180) - 360;                    
    omct_id = interp2(oceans.lon,oceans.lat,oceans.id,new.lon,omct.lat);clear new;
    omct_id(omct_id>0) = 1;
    % Additional ID/Mask polygon (inland lakes)
    polyx = [30.75 27.25 132.3 100.8 38.75 42.75 40.25 30.75]; % Euro-Asia
    polyy = [59.75 63.25 63.25 28.75 36.25 42.75 47.75 59.75];
    ida = inpolygon(omct.lon,omct.lat,polyx,polyy);
    omct_id(ida) = 0;clear polyx polyy ida
    polyx = [0 28.25 44.25 36.25 36.25 24.25 12.25 0 0]; % Africa
    polyy = [29.75 29.75 5.75 -2.25 -14.25 -30.25 5.75 9.75 29.75];
    ida = inpolygon(omct.lon,omct.lat,polyx,polyy);
    omct_id(ida) = 0;clear polyx polyy ida
    polyx = [-75.75 -55.75 -57.75 -39.75 -59.75 -67.75 -67.75 -75.5 -75.75]+360; % South America
    polyy = [5.75 1.45 -6.25 -10.25 -34.25 -42.25 -18.25 -14.25 5.75];
    ida = inpolygon(omct.lon,omct.lat,polyx,polyy);
    omct_id(ida) = 0;clear polyx polyy ida
    polyx = [-155.8 -131.8 -123.8 -119.8 -111.8 -95.75 -95.45 -79.75, ... % North America
            -75.75 -67.75 -61.5 -62.75 -71.75 -74.75 -80.75 -110.8,...
            -119.8 -119.8 -134.8 -155.8]+360; 
    polyy = [ 69.75  65.75  67.75  67.75  65.75  65.75  57.75  49.75,...
             52.0   57.75  53.9   50.75 50.75  41.75  30.75  29.75,...
             38.75  50.75  59.75  69.75];
    ida = inpolygon(omct.lon,omct.lat,polyx,polyy);
    omct_id(ida) = 0;clear polyx polyy ida
else
    % If no grid loaded, compute for whole globe
    omct_id = omct.lon.*0+1;
end
% Compute only for ocean grid cells
ocean_index = 1:length(lon_vec);
ocean_id = omct_id(:)';
in_uniq_ocean = in_uniq(ocean_id==1);
lon_vec_ocean = lon_vec(ocean_id==1);
ocean_index = ocean_index(ocean_id==1);

%% Conversion
% Check how much RAM has the computer (will be used to switch between
% faster and slower conversion): MATLAB ONLY (no 'memory' function)!
v = version; % matlab version returns string ending with ')'
if strcmp(v(end),')')
    [~,mem2] = memory; 
else
    mem2.PhysicalMemory.Available = 0;
end
for i = 1:length(time(:,7))
    % get current year, month and day to read the input file
    [cyear,cmonth,cday,chour] = datevec(time(i,7)); 
    % create file name that will be loaded
    file_name = fullfile(input_path,sprintf('%s%04d-%02d-%02d%s',file_prefix,cyear,cmonth,cday,file_suffix)); 
    try
        % load AOD1B data
        [degree,order,Cnm,Snm] = mGlobe_readAOD1B(file_name,aod1b_type,time(i,7));  
        % Check the max degree of the loaded file and warn user if
        % discrepancy
        if max(degree) < max_deg
           cmax_deg = max(degree);
           fprintf('Oceans: Conversion OMCT model: %s max degree only %3.0f\n',file_name,cmax_deg);
        else
           cmax_deg = max_deg;
        end
        % Get k values for loaded degree
        k = interp1(love.degree,love.k,degree);                             
        % prepare variable
        delta_wc = in_uniq_ocean.*0; 
        % Switch between faster but RAM demanding approach and slower but
        % RAM sufficient approach. The main difference is in the
        % computation of Legendre polynomials that are either computed
        % once and stored (faster) or computed each time (always in octave)
        if mem2.PhysicalMemory.Available/1e+9 > 3 && length(time(:,7)) > 4
            % compute Legendre polynomials for given degree and unique
            % co-latitude (do this only first time)
            if ~exist('Pnm','var')
                for n = 0:cmax_deg
                   Pnm_u = legendre(n,cosinus,'norm');       
                   % Apply full geodetic norm
                   Pnm_u(1,:) = Pnm_u(1,:)*sqrt(2);
                   Pnm_u(2:end,:) = Pnm_u(2:end,:)*2;
                   % Reshape the computed Legendre polynomials to account
                   % for all unique latitude and all longitudes combination
                   Pnm{n+1} = Pnm_u(:,in_uniq_ocean);                                 
                   clear Pnm_u
                end
            end
            % find current degree (in the loaded file)
            for n = 0:cmax_deg
                % find current degree (== all orders/m)
                r = find(degree == n);
                % compute for all orders
                for m = 1:length(r)                                             
                    cur_row = r(m);                                             
                    cur_m = order(cur_row);
                    delta_wc = delta_wc + Pnm{n+1}(cur_m+1,:).*((2*n+1)./(1+k(cur_row))).*(Cnm(cur_row).*cos(cur_m*lon_vec_ocean) + ...
                                                                                           Snm(cur_row).*sin(cur_m*lon_vec_ocean));
                end
                % clear r m cur_row cur_m
            end
        else % not enough memory case
            % compute up to max. degree
            for n = 0:cmax_deg
                % compute Legendre polynomials for given degree and unique co-latitude
                Pnm_u = legendre(n,cosinus,'norm');
                % Apply full geodetic norm
                Pnm_u(1,:) = Pnm_u(1,:)*sqrt(2);                    
                Pnm_u(2:end,:) = Pnm_u(2:end,:)*2;
                % Reshape the computed Legendre polynomials for all lon/lat
                % combination over the ocean
                Pnm = Pnm_u(:,in_uniq_ocean);                                     
                clear Pnm_u
                % find current degree (== all orders/m)
                r = find(degree == n);      
                % compute for all orders
                for m = 1:length(r) 
                    % current row (for current degree and order in the input file)
                    cur_row = r(m);                                      
                    cur_m = order(cur_row);
                    % Sum over all degrees and orders (no filtering = no
                    % additional (W) multiplication)
                    delta_wc = delta_wc + Pnm(cur_m+1,:).*((2*n+1)./(1+k(cur_row))).*(Cnm(cur_row).*cos(cur_m*lon_vec_ocean) + ...
                                                                                    Snm(cur_row).*sin(cur_m*lon_vec_ocean));
                end
            end
        end
        % Reshape vector to matrix: restore for all lon/lat, i.e. also for
        % continents
        delta_wc_all = lon_vec.*NaN;
        delta_wc_all(ocean_index) = delta_wc;
        % Reshape to lon*lat matrix
        omct.obp = reshape(delta_wc_all,size(omct.lat,1),size(omct.lat,2));
        % Convert to mm
        omct.obp = term0.*omct.obp*1000;
        % Store other important info
        omct.spherHarmDegree = cmax_deg;
        omct.input_file = file_name;
        omct.time = time(i,7);
        
        %% Write
		if strcmp(v(end),')')
			save(sprintf('%s%04d%02d%02d_%02d.mat',fullfile(ghc_path,output_prefix),cyear,cmonth,cday,chour),'omct'); % save created grid
        else
			save(sprintf('%s%04d%02d%02d_%02d.mat',fullfile(ghc_path,output_prefix),cyear,cmonth,cday,chour),'omct','-mat7-binary');
		end
		omct.obp = [];  % clear for next time epoch
        if size(time,1) > 2
            out_message = sprintf('Models: converting OMCT model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
        else
            out_message = sprintf('Models: converting OMCT model ...'); % create status message
        end
        set(findobj('Tag','text_status'),'String',out_message); drawnow % write status message
        % Clear used variables used for current time epoch
        clear file_name degree order Cnm Snm cyear cmonth cday chour out_message delta_wc clear r m cur_row cur_m
    catch
        set(findobj('Tag','text_status'),'String',['Models: Could not load or convert: ',file_name]); drawnow
        fprintf('Models: Could not load or convert: %s\n',file_name);
        clear file_name cyear cmonth cday chour
    end
end
end


