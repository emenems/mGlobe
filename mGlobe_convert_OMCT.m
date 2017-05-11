function mGlobe_convert_OMCT(start_calc,end_calc,time_resol,file,ghc_path,input_path,model_version)
%MGLOBE_CONVERT_OMCT Read and covert OMCT spherical harmonics 
% Extract ocean bottom pressure data and save them to matlab format.  
% 
% ASSUMTPION:     
%   The AOD1B data has been downloaded (and unpacked) from:
%           http://isdc.gfz-potsdam.de    OR
%           ftp://podaac-ftp.jpl.nasa.gov/allData/grace/L1B/GFZ/AOD1B/RL05
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
%   model_version  ... Product switch: 5 = oba, 6 = ocn, 7 = atm
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
%                                                                06.01.2015
%   
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
%% Constants
resolution = 1.0;                                                           % spatial resolution of the new grid (deg)
max_deg = 100;                                                              % maximal degree of developement
a = 0.63781364600000E+07;                                                   % major axis (m)
ro_ave = 5517;                                                              % mean density of the Earth
ro_water = 1000;                                                            % mean water density
[omct.lon,omct.lat] = meshgrid(0+resolution/2:resolution:360-resolution/2,... % longitude (deg)
                           -90+resolution/2:resolution:90-resolution/2);    % latitude (deg)
omct.units = 'mm';
term0 = ((a*ro_ave)/(3*ro_water));

switch model_version
    case 5
        aod1b_type = 'oba';
        output_prefix = 'OMCT_oba_6H_';
    case 6
        aod1b_type = 'ocn';
        output_prefix = 'OMCT_ocn_6H_';
    case 7
        aod1b_type = 'atm';                                                 % AOD1B product type
        output_prefix = 'OMCT_atm_6H_';
end
file_prefix = 'AOD1B_';                                                     % file prefix 
file_suffix = file(end-8:end);                                              % file suffix
output_path = ghc_path;                                                     % output path (save to)

file_love = 'mGlobe_DATA_Load_degree_k.txt';                                % load love numbers
file_oceans = 'mGlobe_DATA_OceanGrid.mat';                                  % auxiliary grid for Continents/Oceans

%% Prepare variables
lat_uniq = 90 - unique(omct.lat(:));                                        % unique co-latitude
cosinus = cos(lat_uniq*pi/180);                                                   % cosinus
lon_vec = omct.lon(:)'*pi/180;                                                     % all longitude grid cells
in_uniq(1:length(omct.lat(:))) = 0;                                         % prepare index
for l = 1:length(omct.lat(:))
    in_uniq(l) = find(90-omct.lat(l)==lat_uniq);                            % index used for unique co-latitude. This index vector will be used to re-shape the results computed for unique co-latitude values
end
clear l lat_uniq                                                            % remove used variables

%% Load required files
if ~isempty(file_love)                                                      % load love numbers
    love.data = load(file_love);
    love.degree = love.data(:,1);                                           % degree in first column
    love.k = love.data(:,2);                                                % k in second column !!!!
end

if ~isempty(file_oceans)
    oceans = importdata(file_oceans);                                       % load auxiliary grid used for continents/oceans
    if size(oceans.lon,1) == 1 || size(oceans.lon,2) == 1
       [oceans.lon,oceans.lat] = meshgrid(oceans.lon,oceans.lat); 
    end
    new = omct;
    new.lon(new.lon>=180) = new.lon(new.lon>=180) - 360;                    % transform coordinates/longitude to (-180,180) system  
    omct_id = interp2(oceans.lon,oceans.lat,oceans.id,new.lon,omct.lat);clear new;
    omct_id(omct_id>0) = 1;
    % Additional polygon (inland lakes)
    polyx = [30.75 27.25 132.3 100.8 38.75 42.75 40.25 30.75];              % Euro-Asia
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
    polyx = [-155.8 -131.8 -123.8 -119.8 -111.8 -95.75 -95.45 -79.75 -75.75 -67.75 -61.5 -62.75 -71.75 -74.75 -80.75 -110.8 -119.8 -119.8 -134.8 -155.8]+360; % North America
    polyy = [ 69.75  65.75  67.75  67.75  65.75  65.75  57.75  49.75  52.0   57.75  53.9   50.75 50.75  41.75  30.75  29.75  38.75  50.75  59.75  69.75];
    ida = inpolygon(omct.lon,omct.lat,polyx,polyy);
    omct_id(ida) = 0;clear polyx polyy ida
end

%% for loop
[mem1,mem2] = memory;                                               % futher computation depends on available memory
for i = 1:length(time(:,7))
    [cyear,cmonth,cday,chour] = datevec(time(i,7));                         % get current year, month and day
    file_name = fullfile(input_path,sprintf('%s%04d-%02d-%02d%s',file_prefix,cyear,cmonth,cday,file_suffix)); % create file name that will be loaded
    try
        [degree,order,Cnm,Snm] = mGlobe_readAOD1B(file_name,aod1b_type,time(i,7));  % load AOD1B data
        
        if max(degree) < max_deg
           cmax_deg = max(degree);
           fprintf('Oceans: Conversion OMCT model: %s max degree only %3.0f\n',file_name,cmax_deg);
        else
           cmax_deg = max_deg;
        end
        
    %     Cnm(degree<=1) = 0;                                               % set degree 0 - 1 coefficients to zero
    %     Snm(degree<=1) = 0;                                               % set degree 0 - 1 coefficients to zero
    
        k = interp1(love.degree,love.k,degree);                             % get k values for loaded degree
        delta_wc(1,1:size(omct.lat(:))) = 0;                                % prepare variable
        
        if mem2.PhysicalMemory.Available/1e+9 > 3
            if ~exist('Pnm','var')
                for n = 0:cmax_deg
                   Pnm_u = legendre(n,cosinus,'norm');                          % compute legendre polynomials for given degree and unique co-latitude
                   Pnm_u(1,:) = Pnm_u(1,:)*sqrt(2);                             % full geodetic norm
                   Pnm_u(2:end,:) = Pnm_u(2:end,:)*2;                           % full geodetic norm
                   Pnm{n+1} = Pnm_u(:,in_uniq);                                 % reshape the computed legendre polynomials
                   clear Pnm_u
                end
            end
            for n = 0:cmax_deg
                r = find(degree == n);                                          % find current degree (in the loaded file)
                for m = 1:length(r)                                             % compute for all orders
                    cur_row = r(m);                                             % current row (for current degree and order)
                    cur_m = order(cur_row);
                    delta_wc = delta_wc + term0.*(Pnm{n+1}(cur_m+1,:).*((2*n+1)./(1+k(cur_row))).*(Cnm(cur_row).*cos(cur_m*lon_vec) + ...
                                                                                                   Snm(cur_row).*sin(cur_m*lon_vec)));
                
                end
                % clear r m cur_row cur_m
            end
            
        else                                                                % not enough memory 
            for n = 0:cmax_deg                                              % compute up to max. degree
                Pnm_u = legendre(n,cosinus,'norm');                         % compute legendre polynomials for given degree and unique co-latitude
                Pnm_u(1,:) = Pnm_u(1,:)*sqrt(2);                            % full geodetic norm
                Pnm_u(2:end,:) = Pnm_u(2:end,:)*2;                          % full geodetic norm
                Pnm = Pnm_u(:,in_uniq);                                     % reshape the computed legendre polynomials
                clear Pnm_u
                r = find(degree == n);                                      % find current degree (in the loaded file)
                for m = 1:length(r)                                         % compute for all orders
                    cur_row = r(m);                                         % current row (for current degree and order)
                    cur_m = order(cur_row);
                    delta_wc = delta_wc + term0.*(Pnm(cur_m+1,:).*((2*n+1)./(1+k(cur_row))).*(Cnm(cur_row).*cos(cur_m*lon_vec) + ...
                                                                                              Snm(cur_row).*sin(cur_m*lon_vec)));
                end
                % clear r m Pnm cur_row cur_m
            end
        end
        
        omct.obp = reshape(delta_wc,size(omct.lat,1),size(omct.lat,2));     % reshape vector to matrix
        omct.obp = omct.obp*1000;                                           % convert to mm
        omct.spherHarmDegree = cmax_deg;
        omct.input_file = file_name;
        omct.time = time(i,7);
        
        %% Land/Ocean ID
        if ~isempty(file_oceans)
            omct.obp(omct_id==0) = NaN;                                     % set grid cells over continents to NaN
        end
        %% Write
        save(sprintf('%s%04d%02d%02d_%02d.mat',fullfile(output_path,output_prefix),cyear,cmonth,cday,chour),'omct'); % save created grid
        omct.obp = [];                                                      % remove values for next time epoch
        if size(time,1) > 2
            out_message = sprintf('Models: converting OMCT model ... (%3.0f%%)',100*((i-1)/size(time,1))); % create status message
        else
            out_message = sprintf('Models: converting OMCT model ...'); % create status message
        end
        set(findobj('Tag','text_status'),'String',out_message); drawnow % write status message
        clear file_name degree order Cnm Snm cyear cmonth cday chour out_message delta_wc clear r m cur_row cur_m
        
    catch
        set(findobj('Tag','text_status'),'String',['Models: Could not load or convert: ',file_name]); drawnow
        fprintf('Models: Could not load or convert: %s\n',file_name);
        clear file_name cyear cmonth cday chour
    end
end
end


