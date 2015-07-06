function [degree,order,Cnm,Snm] = mGlobe_readAOD1B(input_file,read_type,time_in)
%MGLOBE_READAOD1B read AOD1B data
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
%                                                   M. Mikolaj, 26.11.2014
%                                                   mikolaj@gfz-potsdam.de

fid = fopen(input_file);                                                    % open file
data_per_file = 4;                                                          % max. number of datasets per file (default 4)
%% Read the header
find_header = fgetl(fid);                                                   % start to read the header
while ~strcmp(find_header,'END OF HEADER')                                  % find the end of header
    find_header = fgetl(fid);
end


%% Find desired product
curr_row = find_header;
degree = [];                                                                % prepare output variables
order = [];
Cnm = [];
Snm = [];
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
                    
    if curr_time == time_in                                                 % read coefficients
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
        data_count = 9e+37;                                                 % stop the loop
    else
        data_count = data_count+1;                                          % continue to next dataset
    end
end

fclose(fid);

end
