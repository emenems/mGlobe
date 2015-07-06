function mGlobe_convert_OTHER(start_calc,end_calc,time_resol,input_name,input_path,output_path,file_type,header_lines)
%MGLOBE_CONVERT_OTHER Read and convert hydrological/ocean models
% Extract required water column data and save it to supported format.
% 
% INPUT:
%   start_calc     ... starting time in matlab format (days)
%                      Example: datenum([2012,1,1,12,0,0]);
%   end_calc       ... finish time in matlab format (days)
%                      Example: datenum([2013,1,1,12,0,0]);
%   time_resol     ... time resolution switcher (not in time units): 1 == 3 hours, 2 == 6 hours,
%                      3 == 12 hours, 4 == 24 hours, 5 == 48 hours, 6 == month
%                      Example: 4
%   input_name     ... one of the input file names (string)
%                      Example: 'OTHERmodel_20120128_12.txt'
%   input_path     ... input file path (string)
%                      Example: fullfile('EXAMPLES','Models');
%   output_path    ... output file path
%                      Example: fullfile('GHM','OTHER');
%   file_type      ... input file type (1-4):
%                      1 = txt: Lat,Lon,h (m)
%                      2 = txt: Lat,Lon,h (mm)
%                      3 = txt: Lon,Lat,h (m)
%                      4 = txt: Lon,Lat,h (mm)
%                      Example: 3
%   header_lines   ... number of header lines in input file (scalar)
%                      Example: 3
% 
% OUTPUT (automatically saved):
%   out_mat        ... structure array (several matrices) containing:
%   out_mat.lon    ... longitude (degrees)
%   out_mat.lat    ... latitude  (degrees)
%   out_mat.time   ... matlab format time (days)
%   out_mat.total  ... water column (mm)
%   out_mat.input_file  ... input file name
%   out_mat.units  ... out_mat.total units
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

%% Load data
set(findobj('Tag','text_status'),'String','Models: converting...');drawnow  % write status message
perc_message = linspace(1,99,size(time,1));                                 % status percentage
for i = 1:size(time,1);
    check_out = 0;                                                          % control variable
    new_input = sprintf('%s_%04d%02d%02d_%02d.%s',...                       % create new file name for loading
        input_name(1:10),time(i,1),time(i,2),time(i,3),time(i,4),input_name(24:end));
    try
        if file_type <= 2                                                   % Lat,Lon,value structure
            data_original = dlmread(fullfile(input_path,new_input),'',header_lines,0); % read the file with given header lines
            out_mat.lat = linspace(min(min(data_original(:,1))),max(max(data_original(:,1))),length(unique(data_original(:,1)))); % equally spaced vector of latitude
            out_mat.lon = linspace(min(min(data_original(:,2))),max(max(data_original(:,2))),length(unique(data_original(:,2)))); % equally spaced vector of longitude
            [out_mat.lon,out_mat.lat] = meshgrid(out_mat.lon,out_mat.lat);  % create new grid
            out_mat.total = griddata(data_original(:,2),data_original(:,1),data_original(:,3),out_mat.lon,out_mat.lat); % interpolate the input values to new grid
			out_mat.input_file = new_input;                                 % store used input file name
			out_mat.units = 'mm';                                           % store used units
        else                                                                % Lon,Lat,value structure
            data_original = dlmread(fullfile(input_path,new_input),'',header_lines,0);
            out_mat.lon = linspace(min(min(data_original(:,1))),max(max(data_original(:,1))),length(unique(data_original(:,1))));
            out_mat.lat = linspace(min(min(data_original(:,2))),max(max(data_original(:,2))),length(unique(data_original(:,2))));
            [out_mat.lon,out_mat.lat] = meshgrid(out_mat.lon,out_mat.lat);
            out_mat.total = griddata(data_original(:,1),data_original(:,2),data_original(:,3),out_mat.lon,out_mat.lat);
			out_mat.input_file = new_input;
			out_mat.units = 'mm';
        end
        if file_type == 1 || file_type == 3
            out_mat.total = out_mat.total*1000;                             % transform meters to mm
        end
        out_mat.time = time(i,7);                                           % store model time
           
    catch exception
        % set(findobj('Tag','text_status'),'String','Models: Could not load input file (check format)'); drawnow % warn user
        % fprintf('Could not load input file %s\n',new_input);
        check_out = 1;
    end

    if check_out == 0;                                                      % write file if the loading was successful
        if time_resol == 6                                                  % Monthly data
            write_file =  sprintf('%s_M_%04d%02d.mat',...
            input_name(1:10),time(i,1),time(i,2));
        else                                                                % Daily/hourly data
            write_file =  sprintf('%s_D_%04d%02d%02d_%02d.mat',...
            input_name(1:10),time(i,1),time(i,2),time(i,3),time(i,4));
        end
        out_file_save = fullfile(output_path,write_file);
        save(out_file_save,'out_mat');										% save the created file
    else
        message = sprintf('Conversion: Could not convert: %s',new_input);
        set(findobj('Tag','text_status'),'String',message);
        fprintf('%s\n',message);
    end 
    clear out_file_save out_mat check_out data_original write_file 
    message = sprintf('Models: converting ...(%02d%%)',round(perc_message(i)));
    set(findobj('Tag','text_status'),'String',message);drawnow;
    clear message

end

