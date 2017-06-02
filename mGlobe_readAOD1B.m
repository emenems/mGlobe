function [degree,order,Cnm,Snm] = mGlobe_readAOD1B(input_file,read_type,time_in)
%MGLOBE_READAOD1B read AOD1B RL05 and RL06 data
% 
% Input:
%   input_file  ... full file name with AOD1B data
%   read_type   ... reading coefficients 'atm'|'ocn'|'glo'|'oba'
%   time_in     ... matlab time of AOD1B coefficents
% 
% Output:
%   degree      ... spherical harmonics degree
%   order       ... spherical harmonics order
%   Cnm         ... Cnm or delta Cnm ceofficients
%   Snm         ... Snm or delta Snm ceofficients
% 
% Example:
% [degree,order,Cnm,Snm] = readAOD1B('AOD1B_2013-01-01_X_05.asc','atm',datenum(2013,1,1,18,0,0));
% 
%                                                   M. Mikolaj
%                                                   mikolaj@gfz-potsdam.de

% Prepare output variables
degree = [];
order = [];
Cnm = [];
Snm = [];

%% Read header info
count_header = 1;
% open file and try to find header info
fid = fopen(input_file,'r');
find_header = fgetl(fid);
% find the end of header
while ~strcmp(find_header,'END OF HEADER')  
    % Find total number of datasets
    if strcmp(find_header(1:19),'NUMBER OF DATA SETS')
        temp = strsplit(find_header,':');
        data_per_file = str2double(temp{2})/4;
    end
    % Find maximum degree
    if strcmp(find_header(1:14),'MAXIMUM DEGREE')  
        temp = strsplit(find_header,':');
        max_deg = str2double(temp{2});
    end
    % Finc total number of records
    if strcmp(find_header(1:22),'NUMBER OF DATA RECORDS')  
        temp = strsplit(find_header,':');
        max_data = str2double(temp{2});
    end
    find_header = fgetl(fid);
    count_header = count_header + 1;
end
% Close file (will be opened later for reading coefficients)
fclose(fid);

%% First try to read only the requested dataset assuming fixed file length
% Get input time stamp and estimate the position of this time stamp inside
% the file
[year,month,day,hour] = datevec(time_in);
% The sorting of type coefficients depend on Release (RL06: max_deg = 100)
% Prepare the position for {'atm','ocn','glo','oba'};
if max_deg == 180
    type_pos = [1 2 3 4];
else
    type_pos = [1 4 2 3];
end
% Estimate the time resolution
time_resol = 24/data_per_file;
% Compute offset between TYPE/hour coef. = get number of coefficients
coef_offset = (max_data/(data_per_file*4)+1); 
% Compute the position for certain hour and all types
time_pos = (type_pos-1)*coef_offset + (count_header + 1) + coef_offset*4*(hour/time_resol);
% Use only user input type
if strcmp(read_type,'atm')
    time_pos = time_pos(1);
elseif strcmp(read_type,'ocn')
    time_pos = time_pos(2);
elseif strcmp(read_type,'glo')
    time_pos = time_pos(3);
elseif strcmp(read_type,'oba')
    time_pos = time_pos(4);
end
% Read using give position
fid = fopen(input_file,'r');
% Check if correct line was loaded by reading TYPE info only
data_check = textscan(fid,'%s',11,'HeaderLines',time_pos-1);
fclose(fid);
if strcmp(data_check{1,1}(end),read_type) && ...
   strcmp(data_check{1,1}(8),sprintf('%02d:00:00',hour)) && ...
   strcmp(data_check{1,1}(7),sprintf('%04d-%02d-%02d',year,month,day))
    % Read only required set of coefficients
    fid = fopen(input_file,'r');
    data = textscan(fid,'%f %f %f %f',coef_offset,'HeaderLines',time_pos);
    degree = data{1,1};
    order = data{1,2};
    Cnm = data{1,3};
    Snm = data{1,4};
    fclose(fid);
else
    % Read using old/slow approach
    fid = fopen(input_file,'r');
    curr_row = fgetl(fid);
    % Find desired product
    data_count = 1;
    while (data_count <= data_per_file) && sum(curr_row) ~= -1
        curr_row = fgetl(fid);  
        curr_type = curr_row(end-2:end);
        if ~strcmp(curr_type,read_type)
            while ~strcmp(curr_type,read_type)
                curr_row = fgetl(fid);
                curr_type = curr_row(end-2:end);
            end
        end
        % Get time
        curr_time = datenum(str2double(curr_row(37:41)),str2double(curr_row(43:44)),str2double(curr_row(46:47)),...
                            str2double(curr_row(49:50)),str2double(curr_row(52:53)),str2double(curr_row(55:56)));
        % read coefficients
        if curr_time == time_in
            data = fgetl(fid);
            while ~strcmp(data(1),'D')
                if data ~= -1
                    degree = vertcat(degree,str2double(data(1:3)));
                    order = vertcat(order,str2double(data(5:7)));
                    Cnm = vertcat(Cnm,str2double(data(9:24)));
                    Snm = vertcat(Snm,str2double(data(25:end)));
                    data = fgetl(fid);
                else
                    data = 'DATA';
                end
            end
            % stop the loop
            data_count = 9e+37;
        else
            % continue to next dataset
            data_count = data_count+1;
        end
    end
    fclose(fid);
end
end % function
