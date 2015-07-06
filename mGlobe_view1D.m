function check_out = mGlobe_view1D(plot_file,plot_option,file_type,time_column,data_column,corr_option,save_option,export_option,print_file,export_file,demean_type,trend_type)
%MGLOBE_VIEW1D Function for the 1D visualization
%   Function is used for the results visualization of 1D data, i.e. time
%   dependent variations.
% 
% Input:
%   plot_file      ...  Full input file names
%   plot_option    ...  switch (0,1) for on/off plotting
%   file_type      ...  switch (1,2,3,4) for txt/mat/xls/tsf input file
%                       types
%   time_column    ...  identifies data column with MATLAB time (days)
%   data_column    ...  identifies data column (channel for tsf format)
%   corr_option    ...  switch (1 to 4) for add/subtract/plot+/plot-
%   save_option    ...  switch (1 to 4) for no/fig/eps/tiff print file
%   export_option  ...  switch (1 to 5) for no/txt/mat/xls/tsf output file
%   print_file     ...  Full output file name
%   demean_type    ...  switch (0,1) for on/off mean value subtraction
%   trend_type     ...  switch (0,X) for on/off trend subtraction
%                       (polynomial Xdegree)
% 
% Output:
%   check_out      ...  check number (1 - OK, 0 - not loaded)
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0


set(findobj('Tag','text_status'),'String','Plot: Visualizing...');drawnow   % write status message
for i = 1:length(file_type);
    if (~isempty(plot_file{i}) && plot_option(i) == 1);
        switch file_type(i)                                                 % switch between supported file formats
            case 1                                                          % txt file format                                
                try
                    data(i).original = load(char(plot_file(i)));
                    data(i).original(isnan(data(i).original(:,data_column(i))),:) = [];
                    time(i).pick = data(i).original(:,time_column(i));
                    data(i).pick = data(i).original(:,data_column(i));
                catch
                    check_out = 0;
                    return
                end
            case 3                                                          % xls file format
                try
                    data(i).original = xlsread(char(plot_file(i)));
                    data(i).original(isnan(data(i).original(:,data_column(i))),:) = [];
                    time(i).pick = data(i).original(:,time_column(i));
                    data(i).pick = data(i).original(:,data_column(i));
                catch 
                    check_out = 0;
                    return
                end
            case 2                                                          % mat file format 
                try
                    data(i).original = importdata(char(plot_file(i)));
                    data(i).original(isnan(data(i).original(:,data_column(i))),:) = [];
                    time(i).pick = data(i).original(:,time_column(i));
                    data(i).pick = data(i).original(:,data_column(i));
                catch 
                    check_out = 0;
                    return
                end
            case 4                                                          % tsoft file format
                try
                    [time(i).pick,tdata] = mGlobe_loadtsf(char(plot_file(i)));
                    data(i).pick = tdata(:,data_column(i));
                catch 
                    check_out = 0;
                    return
                end
                    
        end
        if abs(trend_type(i)) >= 1                                          % remove trend if required
            p_temp = polyfit(time(i).pick(~isnan(data(i).pick)),data(i).pick(~isnan(data(i).pick)),abs(trend_type(i))); % polyfit n-th degree
            data(i).pick = data(i).pick - polyval(p_temp,time(i).pick);
        end
        if demean_type(i) == 1                                              % subtract the average value
            data(i).pick = data(i).pick - mean(data(i).pick(~isnan(data(i).pick)));
        end
        if i == 1 && plot_option(1) == 1                                    % set reference time = gravity time series
            time_gravity = time(i).pick;
            gravity = data(i).pick;
            corr_all(1:length(time_gravity),1:3) = 0;                       % prepare variable for other time series
%         else
%             time_gravity = [];
        end
    else
        data(i).original = 0;
        time(i).pick = 0;
        data(i).pick = 0;
    end
end
clear i
for i = 2:4                                                                 % load/interpolate other time series
   if (exist('time_gravity','var') && ~isempty(plot_file{i})) && sum(plot_option(2:4))~=0
       try
       switch corr_option(i-1)
           case 1
               corr_all(:,i-1) = interp1(time(i).pick,data(i).pick,time_gravity);
           case 2
               corr_all(:,i-1) = -interp1(time(i).pick,data(i).pick,time_gravity);
       end
       catch
           set(findobj('Tag','text_status'),'String','Plot: Could not interpolated given values (series 1 to 3) ...');drawnow
       end
           
   end
end 
if (exist('gravity','var') && sum(corr_option(corr_option<=2)) ~=0) && sum(plot_option(2:4))~=0
    gravity_corrected = gravity + sum(corr_all,2);
else
    gravity_corrected = [];
end
%% Plot
legend_names = {'gravity','series 1','series 2','series 3','hydro','rain','snow','gravity corrected'}; % plot legend
legend_id(1:8) = 0;
figid = figure;
for i = 1:length(file_type)
    if i == 5 && ~isempty(gravity_corrected)
        plot(time_gravity,gravity_corrected,'r')
        legend(legend_names{legend_id==1})
    elseif i == 5
        legend(legend_names{legend_id==1})
    end
    if ~isempty(plot_file(i)) && plot_option(i) == 1;
        switch i
            case 1
                sub_plot1 = subplot(10,1,1:4);
                plot(time_gravity,gravity,'k');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                legend_id(1) = 1;
                time_limit = get(gca,'XLim');
            case 2
                switch corr_option(1)
                    case 3
                        if exist('sub_plot1','var') ~= 1
                            sub_plot1 = subplot(10,1,1:4);
                        end
                        plot(time(i).pick,data(i).pick,'b');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                        legend_id(2) = 1;
                        if exist('time_limit','var') == 1
                            set(gca,'Xlim',time_limit);
                        end
                    case 4  
                        if exist('sub_plot1','var') ~= 1
                            sub_plot1 = subplot(10,1,1:4);
                        end
                        plot(time(i).pick,-data(i).pick,'b');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                        legend_id(2) = 1;
                        if exist('time_limit','var') == 1
                            set(gca,'Xlim',time_limit);
                        end
                    otherwise
                        legend_id(8) = 1;
                end
            case 3
                switch corr_option(2)
                    case 3
                        if exist('sub_plot1','var') ~= 1
                            sub_plot1 = subplot(10,1,1:4);
                        end
                        plot(time(i).pick,data(i).pick,'m');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                        legend_id(3) = 1;
                        if exist('time_limit','var') == 1
                            set(gca,'Xlim',time_limit);
                        end
                    case 4  
                        if exist('sub_plot1','var') ~= 1
                            sub_plot1 = subplot(10,1,1:4);
                        end
                        plot(time(i).pick,-data(i).pick-0.75,'m');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                        legend_id(3) = 1;
                        if exist('time_limit','var') == 1
                            set(gca,'Xlim',time_limit);
                        end
                    otherwise
                        legend_id(8) = 1;
                end
            case 4
                switch corr_option(3)
                    case 3
                        if exist('sub_plot1','var') ~= 1
                            sub_plot1 = subplot(10,1,1:4);
                        end
                        plot(time(i).pick,data(i).pick,'g');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                        legend_id(4) = 1;
                        if exist('time_limit','var') == 1
                            set(gca,'Xlim',time_limit);
                        end
                    case 4  
                        if exist('sub_plot1','var') ~= 1
                            sub_plot1 = subplot(10,1,1:4);
                        end
                        plot(time(i).pick,-data(i).pick-0.75,'g');hold on;grid on;datetick('x','dd/mm/yy','keepticks');
                        legend_id(4) = 1;
                        if exist('time_limit','var') == 1
                            set(gca,'Xlim',time_limit);
                        end
                    otherwise
                        legend_id(8) = 1;
                end
            case 5
                sub_plot2 = subplot(10,1,6:7);
                plot(time(i).pick,data(i).pick,'b');grid on;datetick('x','dd/mm/yy','keepticks');
                legend('hydro');
                if exist('sub_plot1','var') == 1
                   time_limit = get(sub_plot1,'XLim');
                   xlim(time_limit);
                else
                   time_limit = get(gca,'XLim');
                end
            otherwise
                if exist('sub_plot3','var') ~= 1
                    sub_plot3 = subplot(10,1,9:10);
                end
                if (~isempty(plot_file{5}) && ~isempty(plot_file{6})) && exist('plot3','var') ~= 1
                    [plot3 a1 a2] = plotyy(time(6).pick,data(6).pick,time(7).pick,data(7).pick);grid on;
                    datetick(plot3(1),'x','dd/mm/yy','keepticks');
                    set(plot3(2),'XTick',[]);
                    legend('rain','snow');
                    if exist('sub_plot1','var') == 1
                       time_limit = get(sub_plot1,'XLim');
                       set(plot3(1),'XLim',time_limit);
                       set(plot3(2),'XLim',time_limit);
                    else
                       time_limit = get(plot3(1),'XLim');
                    end
                    break
                else
                    plot(time(i).pick,data(i).pick,'b');
                    legend(legend_names{i});
                    if exist('sub_plot1') == 1
                       time_limit = get(sub_plot1,'XLim');
                       xlim(time_limit);
                    else
                       time_limit = get(gca,'XLim');
                    end
                    
                end
                    
        end
    end
end
%% Printing
set(figid, 'PaperPositionMode', 'auto');         
switch save_option
    case 2                                                                  % fig format
        saveas(figid,print_file);
    case 3                                                                  % eps format
        print(figid,'-depsc','-r300',print_file);                           % default resolution = 300 DPI
    case 4                                                                  % tiff format
        print(figid,'-dtiffn','-r300',print_file);
end
check_out = 1;

%% Exporting
set(findobj('Tag','text_status'),'String','Plot: Writing output file...');drawnow
if exist('time_gravity','var')
    out_mat = [time_gravity,datevec(time_gravity),gravity];
    for i = 2:7
        if ~isempty(plot_file{i}) && plot_option(i) ~= 0
            out_mat(:,7+i) = interp1(time(i).pick,data(i).pick,time_gravity);
        else
            out_mat(:,7+i) = zeros(length(time_gravity),1);
        end
    end
    if ~isempty(gravity_corrected)
        out_mat(:,end+1) = gravity_corrected;
    else
        out_mat(:,end+1) = gravity;
    end
    switch export_option
        case 2
            fid = fopen(export_file,'w');
            fprintf(fid,'%% time_matlab   \tDate       \tTime\t\tgravity\t\tseries1\t\tseries2\t\tseries3\t\thydro\t\train\t\tsnow\t\tgravity corrected\n');
            for i = 1:length(time_gravity)
                fprintf(fid,'%12.6f   \t%4d%02d%02d   \t%02d%02d%02d\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\n',...
                    out_mat(i,1),out_mat(i,2),out_mat(i,3),out_mat(i,4),out_mat(i,5),out_mat(i,6),out_mat(i,7),out_mat(i,8),out_mat(i,9),out_mat(i,10),...
                    out_mat(i,11),out_mat(i,12),out_mat(i,13),out_mat(i,14),out_mat(i,15));
            end
            fclose('all');
        case 3
            save(export_file,'out_mat');
        case 4
            table_out = {'time_matlab','year','month','day','hour','minute','second','gravity','series1','series2','series3','hydro','rain','snow','gravity corrected'};
            if size(out_mat,1) >=65536                                      % write xls only if length of the output matrix is below the allowed Excel length
                table_out(2,1) = {'Data to long for excel file-for results,see created txt file'};
                fid = fopen([export_file(1:end-3),'txt'],'w');
                fprintf(fid,'%% time_matlab   \tDate       \tTime\t\tgravity\t\tseries1\t\tseries2\t\tseries3\t\thydro\t\train\t\tsnow\t\tgravity corrected\n');
                for i = 1:length(time_gravity)
                    fprintf(fid,'%12.6f   \t%4d%02d%02d   \t%02d%02d%02d\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\t\t%7.2f\n',...
                    out_mat(i,1),out_mat(i,2),out_mat(i,3),out_mat(i,4),out_mat(i,5),out_mat(i,6),out_mat(i,7),out_mat(i,8),out_mat(i,9),out_mat(i,10),...
                    out_mat(i,11),out_mat(i,12),out_mat(i,13),out_mat(i,14),out_mat(i,15));
                end
                fclose('all');
            else
                table_out(2:size(out_mat,1)+1,1:size(out_mat,2)) = num2cell(out_mat);
            end
            xlswrite(export_file,table_out);
        case 5
            fid = fopen(export_file,'w');
            fprintf(fid,'[TSF-file] v01.0\n\n');
            fprintf(fid,'[UNDETVAL] 1234567.89\n\n');
            out_mat(isnan(out_mat)) = 1234567.89;
            fprintf(fid,'[TIMEFORMAT] DATETIME\n\n');
            fprintf(fid,'[INCREMENT] %8.3f\n\n',(out_mat(2)-out_mat(1))*24*60*60);
            fprintf(fid,'[CHANNELS]\n');
            fprintf(fid,'  Location:GLOBEAT:gravity\n'); 
            fprintf(fid,'  Location:GLOBEAT:series1\n'); 
            fprintf(fid,'  Location:GLOBEAT:series2\n'); 
            fprintf(fid,'  Location:GLOBEAT:series3\n'); 
            fprintf(fid,'  Location:GLOBEAT:hydro\n');  
            fprintf(fid,'  Location:GLOBEAT:rain\n');  
            fprintf(fid,'  Location:GLOBEAT:snow\n');
            fprintf(fid,'  Location:GLOBEAT:gravity_corrected\n\n');
            fprintf(fid,'[UNITS]\n  nm/s^2\n  nm/s^2\n  nm/s^2\n  nm/s^2\n  ?\n  ?\n  ?\n  nm/s^2\n\n');
            fprintf(fid,'[COMMENT]\n\n');
            fprintf(fid,'[COUNTINFO] %8.0f\n\n',size(out_mat,1));
            fprintf(fid,'[DATA]\n');
            for i = 1:size(out_mat,1);
            fprintf(fid,'%04d %02d %02d  %02d %02d %02d   %17.3f %17.3f %17.3f %17.3f %17.3f %17.3f %17.3f %17.3f\n',...
                out_mat(i,2),out_mat(i,3),out_mat(i,4),out_mat(i,5),out_mat(i,6),round(out_mat(i,7)),...
                out_mat(i,8),out_mat(i,9),out_mat(i,10),out_mat(i,11),out_mat(i,12),out_mat(i,13),out_mat(i,14),out_mat(i,15));
            end
            fclose('all');      
    end
                
end
end
function [time,data,channels,units,undetval,increment,countinfo] = mGlobe_loadtsf(input_tsf)
%MGLOBE_LOADTSF load tsf files to matlab
% Input:
%  input_tsf    ...     full file name 
%
% Output:
%   time        ...     matlab time (vector)
%   data        ...     data matrix (for all channels)
%   channels    ...     channels info (cell array)
%   units       ...     units for each channel (cell array)
%   undetval    ...     undefined values (scalar)
%   increment   ...     time step in seconds (scalar)
%   countinfo   ...     total number of observations (scalar)
% 
% Example:
%   [time,data] = mGlobe_loadtsf('VI_gravity.tsf');
% 
%                                                   M.Mikolaj, 27.1.2015
%                                                   mikolaj@gfz-potsdam.de

if isempty(input_tsf)
   [name,path] = uigetfile({'*.tsf'},'Select a TSoft file');
   input_tsf = fullfile(path,name);
end
%% Initialize
time = [];
data = [];
channels = [];
units = [];
undetval = [];
increment = [];
countinfo = [];

%% Get undetval
cr = 0;
num_chan = 1;
fid = fopen(input_tsf);
row = fgetl(fid);cr = cr+1;
while ischar(row)
    if length(row) >= 6
        switch row(1:6)
            case '[UNDET'
                undetval = str2double(row(11:end));
            case '[INCRE'
                increment = str2double(row(12:end));
            case '[CHANN'
                row = fgetl(fid);cr = cr+1;
                while ~isempty(row)
                    channels{num_chan,1} = row;num_chan = num_chan+1;
                    row = fgetl(fid);cr = cr+1;
                end
            case '[UNITS'
                num_mm = 1;
                row = fgetl(fid);cr = cr+1;
                while ~isempty(row)
                    units{num_mm,1} = row;num_mm = num_mm+1;
                    row = fgetl(fid);cr = cr+1;
                end
            case '[COUNT'
                countinfo = str2double(row(12:end));
            case '[DATA]'
                data_start = cr;
%                 while ~isempty(row) && ~ischar(row)
%                     row = fgetl(fid);cr = cr+1;
%                     di = di + 1;
%                     time(di,1) = datenum(str2double(row(1:4)),str2double(row(5:7)),str2double(row(8:10)),...
%                                          str2double(row(11:13)),str2double(row(14:16)),str2double(row(17:19)));
%                 end
                break
        end
    end
    row = fgetl(fid);cr = cr+1;
end
fclose(fid);
% create format specification
formatSpec = '%d%d%d%d%d%d';
for i = 1:length(channels);
   formatSpec = [formatSpec,'%f'];
end
try 
    % Get Data
    
    try                                                                     % assumed, file contains COUNTINFO 
        fid = fopen(input_tsf,'r');
        for i = 1:data_start
            row = fgetl(fid);
        end
        if isempty(countinfo)
            count = 0;
            row = fgetl(fid);
            if isempty(row)
                while isempty(row)
                    row = fgetl(fid);
                    data_start = data_start + 1;
                end
            end
            while ischar(row)
                row = fgetl(fid);
                count = count + 1;
            end
            fclose(fid);
            fid = fopen(input_tsf,'r');
            countinfo = count;
            for i = 1:data_start
                row = fgetl(fid);
            end
        end
        dataArray = textscan(fid, formatSpec, countinfo);
        time = datenum(double(dataArray{1,1}),double(dataArray{1,2}),double(dataArray{1,3}),double(dataArray{1,4}),double(dataArray{1,5}),double(dataArray{1,6}));
        data = cell2mat(dataArray(7:end));
        if ~isempty(undetval)
            data(data == undetval) = NaN;
        end
    catch   
        dataArray = dlmread(input_tsf,'',data_start,0); % warning no footer info are allowed
        time = datenum(dataArray(:,1:6));
        data = dataArray(:,7:end);
        if ~isempty(undetval)
            data(data == undetval) = NaN;
        end
    end
        
%     % Get footer info
%     row = fgetl(fid);
%     while ischar(row)
%         if length(row) >= 5
%         end
%         row = fgetl(fid);
%     end
catch
    %fprintf('Could not load the required file. Checkt the format (file must contain: COUNTINFO, CHANNEL, UNITS, UNDETVAL)\n');
end
end   
