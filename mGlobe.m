function mGlobe(in_switch)
%MGLOBE MAIN FUNCTION GENERATING GUI FOR THE mGlobE TOOLBOX
% mGlobe toolbox allows user to calculate the global hydrological, 
% atmospheric and non-tidal ocean loading effects. All required inputs
% can be converted via the provided GUI. Obtained results can
% be visualized using the provided 1D and 2D plotting functions. 
% Please read the mGlobe_USER_MANUAL.pdf file before using mGlobe.
% This mGlobe version is designed for reviewing purposes (submitted to 
% Computers & Geosciences Journal). Always check for the latest 
% version: github.com/emenems/mGlobe.
% 
% This function requires/calls following sub-functions stored in the 
% current folder:
%   mGlobe_calc_Atmo_ERA
%   mGlobe_calc_Atmo_MERRA
%   mGlobe_calc_atmo_loading
%   mGlobe_calc_atmo_newton
%   mGlobe_calc_Hydro
%   mGlobe_calc_Ocean
%   mGlobe_convert_DEM
%   mGlobe_convert_ECCO
%   mGlobe_convert_ERA
%   mGlobe_convert_GRACE_tellus
%   mGlobe_convert_OMCT
%   mGlobe_convert_OTHER
%   mGlobe_convert_GLDAS
%   mGlobe_elip2sphere
%   mGlobe_elip2xyz
%   mGlobe_Global
%   mGlobe_interpolation
%   mGlobe_Local
%   mGlobe_readAOD1B
%   mGlobe_tesseroid
%   mGlobe_view1D
%   mGlobe_view2D
% 
% Additionally, following data files are required:
%   mGlobe_DATA_dgE_Atmo.txt
%   mGlobe_DATA_dgE_Hydro.txt
%   mGlobe_DATA_OceanGrid.mat
%   mGlobe_DATA_Load_degree_k.txt
%   mGlobe_PATH_Settings.txt
% 
% See mGlobe_USER_MANUAL.pdf for details
% 
% System requirements:
%   MATLAB R2012a or later
%   Mapping toolbox
%   Statistics toolbox
%   4 GB of RAM (8 GB recommended)
% 
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de

    if nargin == 0
%% GUI
        check_open_window = get(findobj('Tag','radio_switch_hydr'),'Value'); % check if no mGlobe window is open (works only with one window)
        if numel(check_open_window)>0
            fprintf('Please use only one app window in Matlab\n')           % send message to command window
        else
        % WINDOW
        F1 = figure('Position',[520 300 561 402],...                        % create main window
                    'Tag','main_manu','Menubar','none','Resize','off',...
                    'NumberTitle','off','Color',[0.941 0.941 0.941],...
                    'Name','mGlobe: global mass effects in hydrology, atmosphere and oceans on surface gravity');
        % SWITCHER                                                          % switch between main panels
        uicontrol(F1,'Position',[16,365,72,25],'Style','Radiobutton',...
            'Value',1,'FontSize',10,'FontWeight','bold','String','Hydro',...
            'Tag','radio_switch_hydr','CallBack','mGlobe switch_hydro');
        uicontrol(F1,'Position',[233,365,85,25],'Style','Radiobutton',...
            'Value',0,'FontSize',10,'FontWeight','bold','String','Atmo',...
            'Tag','radio_switch_atmo','CallBack','mGlobe switch_atmo');
        uicontrol(F1,'Position',[463,365,85,25],'Style','Radiobutton',...
            'Value',0,'FontSize',10,'FontWeight','bold','String','Plot',...
            'Tag','radio_switch_view','CallBack','mGlobe switch_view');
        uicontrol(F1,'Position',[348,365,85,25],'Style','Radiobutton',...
            'Value',0,'FontSize',10,'FontWeight','bold','String','Models',...
            'Tag','radio_switch_down','CallBack','mGlobe switch_down');
        uicontrol(F1,'Position',[118,365,85,25],'Style','Radiobutton',...
            'Value',0,'FontSize',10,'FontWeight','bold','String','Ocean',...
            'Tag','radio_switch_ocean','CallBack','mGlobe switch_ocean');
        
        % DESCRIPTION                                                       % lower information panel
        uicontrol(F1,'Position',[9,1,492,26],'Style','Text',...
            'FontSize',8,'FontAngle','italic','String','Set the global hydrological effect',...
            'Tag','text_status');
        uicontrol(F1,'Position',[499,6,55,19],'Style','Pushbutton',...      % Help button: no callbak
            'FontSize',8,'String','mGlobe','ForegroundColor','blue',...
            'Tag','push_name','FontAngle','italic','FontWeight','bold');
        
        % PANEL_MAIN                                                        % sub-panels associated with SWITCHER, i.e. PANEL_HYDRO, PANEL_ATMO ...
        p1 = uipanel('Title','Continental Water Storage Effect','Units','characters',...
                    'Position',[1.4,1.923,109.4,26],...
                    'Tag','uipanel_main_hydr','Visible','on','FontSize',9);
        p2 = uipanel('Title','Visualize 1D and 2D results','Units','characters',...
                    'Position',[1.4,1.923,109.4,26],...
                    'Tag','uipanel_main_view','Visible','off','FontSize',9);
        p3 = uipanel('Title','Conversion of Global Model','Units','characters',...
                    'Position',[1.4,1.923,109.4,26],...
                    'Tag','uipanel_main_down','Visible','off','FontSize',9);
        p4 = uipanel('Title','Global Atmospheric Effect','Units','characters',...
                    'Position',[1.4,1.923,109.4,26],...
                    'Tag','uipanel_main_atmo','Visible','off','FontSize',9);
        p5 = uipanel('Title','Non-tidal Ccean Loading Effect','Units','characters',...
                    'Position',[1.4,1.923,109.4,26],...
                    'Tag','uipanel_main_ocean','Visible','off','FontSize',9);
                
        % PANEL_HYDRO       
        p1_1 = uipanel(p1,'Title','Position','Units','characters',...
            'Position',[2.2,15.308,51,9.077],...
            'Tag','uipanel_calc_position');       
        p1_2 = uipanel(p1,'Title','Time','Units','characters',...
            'Position',[56 15.308 51 9.077],...
            'Tag','uipanel_calc_time');       
        p1_3 = uipanel(p1,'Title','Topography','Units','characters',...
            'Position',[2.2 5.408 104.8 3.9],...
            'Tag','uipanel_calc_topography');       
        p1_4 = uipanel(p1,'Title','Output','Units','characters',...
            'Position',[2.2 0.462 62.4 4.93],...
            'Tag','uipanel_calc_output');       
        p1_5 = uipanel(p1,'Title','Hydrological Model','Units','characters',...
            'Position',[2.2 9.308 104.8 5.9],...
            'Tag','uipanel_calc_model');
        uicontrol(p1,'Style','PushButton','Units','characters',...
            'position',[92.6 1.385 13.6 2.077],...
            'Tag','push_calc_calc','CallBack','mGlobe hydro_calc',...
            'String','Calculate','FontWeight','bold');
        uicontrol(p1,'Style','Text','Units','characters',...
            'position',[64.4 1.846 16.6 1.075],...
            'String','Threshold (deg):');
        uicontrol(p1,'Style','Edit','Units','characters',...
            'position',[81.4 1.538 10.2 1.692],'BackgroundColor','white',...
            'String','0.1','Tag','edit_hydro_treshold');
        
        % PANEL_VIEW
        p2_1 = uipanel(p2,'Title','1D data','Units','characters',...
            'Position',[1.8 6.462 105.2 18],...
            'Tag','uipanel_view_1D');  
        p2_2 = uipanel(p2,'Title','2D data','Units','characters',...
            'Position',[1.8 1.077 105.2 5],...
            'Tag','uipanel_view_2D'); 
        
        % PANEL_MODEL
        p3_1 = uipanel(p3,'Title','GLDAS/MERRA convert','Units','characters',...
            'Position',[2.2 12.308 104.8 3.923],...
            'Tag','uipanel_down_gldas'); 
        p3_2 = uipanel(p3,'Title','ERA/NCEP convert (netCDF)','Units','characters',...
            'Position',[2.2 8.308 104.8 3.923],...
            'Tag','uipanel_conv_era'); 
        p3_3 = uipanel(p3,'Title','Other Models (e.g. WGHM)','Units','characters',...
            'Position',[2.2 4.308 104.8 3.923],...
            'Tag','uipanel_conv_other');
        p3_4 = uipanel(p3,'Title','Time','Units','characters',...
            'Position',[2.2 16.308 51 8.308],...
            'Tag','uipanel_conv_other');
        p3_5 = uipanel(p3,'Title','DEM conversion','Units','characters',...
            'Position',[2.2 0.308 104.8 3.923],...
            'Tag','uipanel_conv_dem');
        p3_6 = uipanel(p3,'Title','GRACE-TELLUS (grid)','Units','characters',...
            'Position',[54.2 16.308 52.8 8.308],...
            'Tag','uipanel_conv_grace');
        
        % PANEL_OCEAN       
        p5_1 = uipanel(p5,'Title','Position','Units','characters',...
            'Position',[2.2 15.308 51 9.077],...
            'Tag','uipanel_calc_position_ocean');       
        p5_2 = uipanel(p5,'Title','Time','Units','characters',...
            'Position',[56.0 15.308 51 9.077],...
            'Tag','uipanel_calc_time_ocean');       
        p5_3 = uipanel(p5,'Title','Model setup and conversion','Units','characters',...
             'Position',[2.2 5.408 104.8 9.8],...
             'Tag','uipanel_calc_topography_ocean');       
        p5_4 = uipanel(p5,'Title','Output','Units','characters',...
            'Position',[2.2 0.462 62.4 4.93],...
            'Tag','uipanel_calc_output_ocean');     
        uicontrol(p5,'Style','PushButton','Units','characters',...
            'position',[92.6 1.385 13.6 2.077],...
            'Tag','push_calc_calc_ocean','CallBack','mGlobe ocean_calc',...
            'String','Calculate','FontWeight','bold');
        uicontrol(p5,'Style','Text','Units','characters',...
            'position',[64.4 1.846 16.6 1.075],...
            'String','Threshold (deg):');
        uicontrol(p5,'Style','Edit','Units','characters',...
            'position',[81.4 1.538 10.2 1.692],'BackgroundColor','white',...
            'String','0.1','Tag','edit_ocean_treshold');
        
        % PANEL_ATMO       
        p4_1 = uipanel(p4,'Title','Position','Units','characters',...
            'Position',[2.2,15.308,51,9.077],...
            'Tag','uipanel_calc_position_atmo');       
        p4_2 = uipanel(p4,'Title','Time','Units','characters',...
            'Position',[56 15.308 51 9.077],...
            'Tag','uipanel_calc_time_atmo');             
        p4_4 = uipanel(p4,'Title','Output','Units','characters',...
            'Position',[2.2 0.462 62.4 4.93],...
            'Tag','uipanel_calc_output_atmo');       
        p4_3 = uipanel(p4,'Title','Atmospheric model data','Units','characters',...
            'Position',[2.2 5.692 104.8 9.538],...
            'Tag','uipanel_calc_model_atmo');
        uicontrol(p4,'Style','PushButton','Units','characters',...
            'position',[65.0 1.385 20.5 2.077],...
            'Tag','push_calc_atmo_era','CallBack','mGlobe atmo_calc_era',...
            'String','Calculate ERA','FontWeight','bold');
        uicontrol(p4,'Style','PushButton','Units','characters',...
            'position',[86.1 1.385 20.5 2.077],...
            'Tag','push_calc_atmo_merra','CallBack','mGlobe atmo_calc_merra',...
            'String','Calculate MERRA','FontWeight','bold','FontAngle','italic');
        
        % IN HYDRO
        % position
        uicontrol(p1_1,'units','characters','Position',[0.9,5.615,10.0,1.077],...
                    'Style','Text','String','Latitude');
        uicontrol(p1_1,'units','characters','Position',[1.7,3.385,10.0,1.077],...
                    'Style','Text','String','Longitude');
        uicontrol(p1_1,'units','characters','Position',[0.02,1.154,10.0,1.077],...
                    'Style','Text','String','Height');
        uicontrol(p1_1,'units','characters','Position',[36.8,5.692,5.6,1],...
                    'Style','Text','String','deg');
        uicontrol(p1_1,'units','characters','Position',[36.8,1.385,5.6,1],...
                    'Style','Text','String','m');
        uicontrol(p1_1,'units','characters','Position',[36.8,3.538,5.6,1],...
                    'Style','Text','String','deg');
        uicontrol(p1_1,'units','characters','Position',[13.6,5.308,23.2,1.692],...
                    'Style','Edit','Tag','edit_hydro_pos_lat',...
                    'String','48.24885','BackgroundColor','white');
        uicontrol(p1_1,'units','characters','Position',[13.6,3.077,23.2,1.692],...
                    'Style','Edit','Tag','edit_hydro_pos_lon',...
                    'String','16.35650','BackgroundColor','white')
        uicontrol(p1_1,'units','characters','Position',[13.6,0.846,23.2,1.692],...
                    'Style','Edit','Tag','edit_hydro_pos_hei',...
                    'String','192.70','BackgroundColor','white')
        % time
        uicontrol(p1_2,'units','characters','Position',[1.4,5.385,7,1.077],...
                    'Style','Text','units','characters',...
                    'String','Start');
        uicontrol(p1_2,'units','characters','Position',[1,3.154,7,1.077],...
                    'Style','Text','String','End');
        uicontrol(p1_2,'units','characters','Position',[1.4,0.923,7,1.077],...
                    'Style','Text','String','Step');
        uicontrol(p1_2,'units','characters','Position',[10.8 7.054 5.6,1],...
                    'Style','Text','String','year');
        uicontrol(p1_2,'units','characters','Position',[18.8 7.054 5.9,1],...
                    'Style','Text','String','month');
        uicontrol(p1_2,'units','characters','Position',[28.4 7.054 5.6,1],...
                    'Style','Text','String','day');
        uicontrol(p1_2,'units','characters','Position',[37 7.054 5.6,1],...
                    'Style','Text','String','hour');
        uicontrol(p1_2,'units','characters','Position',[10.4 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_hydro_time_start_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p1_2,'units','characters','Position',[19.2 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_hydro_time_start_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p1_2,'units','characters','Position',[28 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_hydro_time_start_day',...
                    'String','01','BackgroundColor','white')
        uicontrol(p1_2,'units','characters','Position',[36.8 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_hydro_time_start_hour',...
                    'String','12','BackgroundColor','white')
        uicontrol(p1_2,'units','characters','Position',[10.4,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_hydro_time_end_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p1_2,'units','characters','Position',[19.2,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_hydro_time_end_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p1_2,'units','characters','Position',[28,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_hydro_time_end_day',...
                    'String','03','BackgroundColor','white');
        uicontrol(p1_2,'units','characters','Position',[36.8,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_hydro_time_end_hour',...
                    'String','12','BackgroundColor','white');
        uicontrol(p1_2,'units','characters','Position',[10.4,0.615,18,1.692],...
                    'Style','Popupmenu','Tag','popup_hydro_time_step',...
                    'String','3 hours|6 hours|12 hours|Day|2 days|Month','Value',6,'BackgroundColor','white');
        % Topography
        uicontrol(p1_3,'units','characters','Position',[1.4,0.692,12,1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','DEM','UserData',[],...
                    'Tag','push_hydro_topo_dem','CallBack',...
                    'mGlobe load_hydro_dem');
        uicontrol(p1_3,'units','characters','Position',[13.6,0.538,77.6,1.769],...
                    'Style','Text','Tag','text_hydro_topo_dem',...
                    'String','If required, load DEM up to 1.05 deg from gravimeter');
        uicontrol(p1_3,'units','characters','Position',[91.2 0.692 12 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Show','UserData',[],...
                    'Tag','push_hydro_topo_show','CallBack',...
                    'mGlobe load_hydro_show');
        % Output
        uicontrol(p1_4,'units','characters','Position',[1.6 1.88 12 1.689],...
                    'Style','Pushbutton','units','characters',...
                    'String','File','UserData','output.txt',...
                    'Tag','push_hydro_out','CallBack',...
                    'mGlobe load_hydro_out');
        uicontrol(p1_4,'units','characters','Position',[14 1.83 46.6 1.538],...
                    'Style','Text','Tag','text_hydro_out',...
                    'String','Chose your output file (output.txt default)');
        uicontrol(p1_4,'units','characters','Position',[1.8,0.215,8,1.692],...
                    'Style','Checkbox','Tag','check_hydro_xls',...
                    'String','xls','Value',0);
        uicontrol(p1_4,'units','characters','Position',[10.6,0.215,7.6,1.692],...
                    'Style','Checkbox','Tag','check_hydro_txt',...
                    'String','txt','Value',1);
        uicontrol(p1_4,'units','characters','Position',[19,0.215,7.6,1.692],...
                    'Style','Checkbox','Tag','check_hydro_tsf',...
                    'String','tsf','Value',0);
        uicontrol(p1_4,'units','characters','Position',[28 0.231 22.8 1.692],...
                    'Style','Checkbox','Tag','check_hydro_average',...
                    'String','subtract average','Value',0);
        uicontrol(p1_4,'units','characters','Position',[46 0.231 13.5 1.692],...
                    'Style','Pushbutton','units','characters',...
                    'String','SHP (optional)','UserData',0,...
                    'Tag','push_hydro_load_shp','CallBack','mGlobe push_hydro_shp');
        % Hydrological model
        uicontrol(p1_5,'units','characters','Position',[0.6 3.0 7 1.077],...
                    'Style','Text','String','Model');
        uicontrol(p1_5,'units','characters','Position',[13.6 2.769 27.6 1.692],...
                    'Style','Popupmenu','Tag','popup_hydro_model',...
                    'String','GLDAS/CLM|GLDAS/MOS|GLDAS/NOAH (0.25�)|GLDAS/NOAH (1�)|GLDAS/VIC|ERA Interim|MERRA (assimilation)|Other|GRACE|NCEP Reanalysis-2|MERRA2|NCEP Reanalysis-1|GLDASv2.1/NOAH (0.25�)',...
                    'Value',3,'BackgroundColor','white','CallBack','mGlobe select_hydro_model');
        uicontrol(p1_5,'units','characters','Position',[67.2 3 8 1.077],...
                    'Style','Text','Tag','push_hydro_model_path',...
                    'String','Path:','UserData',fullfile('GHM','NOAH025'));
        uicontrol(p1_5,'units','characters','Position',[77 3 18 1.08],...
                    'Style','Text','Tag','text_hydro_model_path',...
                    'String','/GHM/NOAH025/');
        uicontrol(p1_5,'units','characters','Position',[42 3.0 8 1.08],...
                    'Style','Text','Tag','text_hydro_model_layer',...
                    'String','Layer:');
        uicontrol(p1_5,'units','characters','Position',[51 2.769 14 1.692],...
                    'Style','Popupmenu','Tag','popup_hydro_model_layer',...
                    'String','total|soilm1|soilm2|soilm3|soilm4|swe',...
                    'Value',1,'BackgroundColor','white','CallBack','mGlobe select_hydro_layer');
        uicontrol(p1_5,'units','characters','Position',[0.8 0.85 9.0 1.08],...
                    'Style','Text','String','Exclude:');
        uicontrol(p1_5,'units','characters','Position',[57.5 0.85 18.5 1.08],...
                    'Style','Text','String','Mass coserv.:');
        uicontrol(p1_5,'units','characters','Position',[10.5 0.55 15 1.769],...
                    'Style','Checkbox','String','Greenland','Value',1,...
                    'Tag','check_hydro_model_green');
        uicontrol(p1_5,'units','characters','Position',[27 0.55 15 1.769],...
                    'Style','Checkbox','String','Antarctica','Value',0,...
                    'Tag','check_hydro_model_perma');
        uicontrol(p1_5,'units','characters','Position',[75 0.538 27.4 1.769],...
                    'Style','Popupmenu','String','off|Ocean layer (from mass excess)|Continent + Ocean model','Value',2,...
                    'Tag','check_hydro_model_ocean','BackgroundColor','white');
        uicontrol(p1_5,'units','characters','Position',[43.0 0.6 14.7 1.569],...
                    'Style','Pushbutton','units','characters',...
                    'String','Include only','UserData',[],...
                    'Tag','push_hydro_include','CallBack',...
                    'mGlobe load_hydro_include');
        
        % IN VIEW
        % 1D data
        uicontrol(p2_1,'units','characters','Position',[31.6 15.883 10.4 1.154],...
                    'Style','Text','String','Load file');
        uicontrol(p2_1,'units','characters','Position',[43.6 15.883 9.8 1.154],...
                    'Style','Text','String','- trend');
        uicontrol(p2_1,'units','characters','Position',[53 15.883 9.6 1.154],...
                    'Style','Text','String','- mean');
        uicontrol(p2_1,'units','characters','Position',[62.8 15.883 13.6 1.154],...
                    'Style','Text','String','Time column');
        uicontrol(p2_1,'units','characters','Position',[76.6 15.883 13.65 1.154],...
                    'Style','Text','String','Data column');
        uicontrol(p2_1,'units','characters','Position',[89.2 15.883 13.6 1.154],...
                    'Style','Text','String','File type');
        uicontrol(p2_1,'units','characters','Position',[1.6 14.308 22 1.077],...
                    'Style','Checkbox','String','Gravity variation',...
                    'Tag','check_view_gravity','Value',0);
        uicontrol(p2_1,'units','characters','Position',[1.6,12.077 15,1.077],...
                    'Style','Checkbox','String','Series 1',...
                    'Tag','check_view_ghc','Value',0);
        uicontrol(p2_1,'units','characters','Position',[1.6,9.846,15,1.077],...
                    'Style','Checkbox','String','Series 2',...
                    'Tag','check_view_ntol','Value',0);
        uicontrol(p2_1,'units','characters','Position',[1.6 7.615 15 1.077],...
                    'Style','Checkbox','String','Series 3',...
                    'Tag','check_view_atme','Value',0);
        uicontrol(p2_1,'units','characters','Position',[1.6 5.385 22 1.077],...
                    'Style','Checkbox','String','Hydro data',...
                    'Tag','check_view_hydro','Value',0);
        uicontrol(p2_1,'units','characters','Position',[1.6 3.154 20 1.077],...
                    'Style','Checkbox','String','Rain rate/snow',...
                    'Tag','check_view_rain','Value',0);
        uicontrol(p2_1,'units','characters','Position',[14.6 11.692 13 1.692],...
                    'Style','Popupmenu','String','Add|Subtract|Plot|Plot-',...
                    'Tag','popup_view_add_ghc','Value',3,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[14.6 9.538 13 1.692],...
                    'Style','Popupmenu','String','Add|Subtract|Plot|Plot-',...
                    'Tag','popup_view_add_ntol','Value',3,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[14.6 7.385 13 1.692],...
                    'Style','Popupmenu','String','Add|Subtract|Plot|Plot-',...
                    'Tag','popup_view_add_atme','Value',3,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[91 14.154 11.2 1.692],...
                    'Style','Popupmenu','String','txt|mat|xls|tsf',...
                    'Tag','popup_view_type_gravity','Value',1,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[91 11.923 11.2 1.692],...
                    'Style','Popupmenu','String','txt|mat|xls|tsf',...
                    'Tag','popup_view_type_ghc','Value',1,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[91 9.692 11.2 1.692],...
                    'Style','Popupmenu','String','txt|mat|xls|tsf',...
                    'Tag','popup_view_type_ntol','Value',1,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[91 7.462 11.2 1.692],...
                    'Style','Popupmenu','String','txt|mat|xls|tsf',...
                    'Tag','popup_view_type_atme','Value',1,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[91 5.231 11.2 1.692],...
                    'Style','Popupmenu','String','txt|mat|xls|tsf',...
                    'Tag','popup_view_type_hydro','Value',1,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[91 3 11.2 1.692],...
                    'Style','Popupmenu','String','txt|mat|xls|tsf',...
                    'Tag','popup_view_type_rain','Value',1,'BackgroundColor','white');
        uicontrol(p2_1,'units','characters','Position',[71.6 0.615 17 1.692],...
                    'Style','Popupmenu','String','Print as|fig|eps|tiff',...
                    'Tag','popup_view_printas_1D','Value',1,'BackgroundColor','white',...
                    'UserData',[]);
        uicontrol(p2_1,'units','characters','Position',[53 0.615 17 1.692],...
                    'Style','Popupmenu','String','Export as|txt|mat|xls|tsf',...
                    'Tag','popup_view_exportas_1D','Value',1,'BackgroundColor','white',...
                    'UserData',[]);
        uicontrol(p2_1,'units','characters','Position',[28.8 14 15 1.769],...
                    'Style','Pushbutton','String','gravity',...
                    'Tag','push_view_file_gravity','UserData',[],...
                    'CallBack','mGlobe load_view_gravity');
        uicontrol(p2_1,'units','characters','Position',[28.8 11.769 15 1.769],...
                    'Style','Pushbutton','String','e.g. GHE',...
                    'Tag','push_view_file_ghc','UserData',[],...
                    'CallBack','mGlobe load_view_ghc');
        uicontrol(p2_1,'units','characters','Position',[28.8 9.538 15 1.769],...
                    'Style','Pushbutton','String','e.g. NTOL',...
                    'Tag','push_view_file_ntol','UserData',[],...
                    'CallBack','mGlobe load_view_ntol');
        uicontrol(p2_1,'units','characters','Position',[28.8 7.308 15 1.769],...
                    'Style','Pushbutton','String','e.g. ATMO',...
                    'Tag','push_view_file_atme','UserData',[],...
                    'CallBack','mGlobe load_view_atme');
        uicontrol(p2_1,'units','characters','Position',[28.8 5.077 15 1.769],...
                    'Style','Pushbutton','String','Hydro',...
                    'Tag','push_view_file_hydro','UserData',[],...
                    'CallBack','mGlobe load_view_hydro');
        uicontrol(p2_1,'units','characters','Position',[28.8 2.846 15 1.769],...
                    'Style','Pushbutton','String','Rain / Snow',...
                    'Tag','push_view_file_rain','UserData',[],...
                    'CallBack','mGlobe load_view_rain');
        uicontrol(p2_1,'units','characters','Position',[90.2 0.538 12 1.769],...
                    'Style','Pushbutton','String','1D plot','FontWeight','bold',...
                    'Tag','push_view_1Dplot','UserData',[],...
                    'CallBack','mGlobe view_1Dplot');
        uicontrol(p2_1,'units','characters','Position',[64.6 14.154 11.2 1.615],...
                    'Style','Edit','String','1','BackgroundColor','white',...
                    'Tag','edit_view_time_gravity');
        uicontrol(p2_1,'units','characters','Position',[64.6 11.923 11.2 1.615],...
                    'Style','Edit','String','1','BackgroundColor','white',...
                    'Tag','edit_view_time_ghc');
        uicontrol(p2_1,'units','characters','Position',[64.6 9.692 11.2 1.615],...
                    'Style','Edit','String','1','BackgroundColor','white',...
                    'Tag','edit_view_time_ntol');
        uicontrol(p2_1,'units','characters','Position',[64.6 7.462 11.2 1.615],...
                    'Style','Edit','String','1','BackgroundColor','white',...
                    'Tag','edit_view_time_atme');
        uicontrol(p2_1,'units','characters','Position',[64.6 5.231 11.2 1.615],...
                    'Style','Edit','String','1','BackgroundColor','white',...
                    'Tag','edit_view_time_hydro');
        uicontrol(p2_1,'units','characters','Position',[64.6 3 11.2 1.615],...
                    'Style','Edit','String','1','BackgroundColor','white',...
                    'Tag','edit_view_time_rain');
        uicontrol(p2_1,'units','characters','Position',[77.8 14.154 11.2 1.615],...
                    'Style','Edit','String','4','BackgroundColor','white',...
                    'Tag','edit_view_data_gravity');
        uicontrol(p2_1,'units','characters','Position',[77.8 11.923 11.2 1.615],...
                    'Style','Edit','String','4','BackgroundColor','white',...
                    'Tag','edit_view_data_ghc');
        uicontrol(p2_1,'units','characters','Position',[77.8 9.692 11.2 1.615],...
                    'Style','Edit','String','4','BackgroundColor','white',...
                    'Tag','edit_view_data_ntol');
        uicontrol(p2_1,'units','characters','Position',[77.8 7.461 11.2 1.615],...
                    'Style','Edit','String','4','BackgroundColor','white',...
                    'Tag','edit_view_data_atme');
        uicontrol(p2_1,'units','characters','Position',[77.8 5.230 11.2 1.615],...
                    'Style','Edit','String','2','BackgroundColor','white',...
                    'Tag','edit_view_data_hydro');
        uicontrol(p2_1,'units','characters','Position',[77.8 3 5 1.615],...
                    'Style','Edit','String','2','BackgroundColor','white',...
                    'Tag','edit_view_data_rain');
        uicontrol(p2_1,'units','characters','Position',[84 3 5 1.615],...
                    'Style','Edit','String','3','BackgroundColor','white',...
                    'Tag','edit_view_data_snow');
        uicontrol(p2_1,'units','characters','Position',[45.4 14.154 7 1.615],...
                    'Style','Edit','String','0','BackgroundColor','white',...
                    'Tag','edit_view_trend_gravity');
        uicontrol(p2_1,'units','characters','Position',[45.4 11.923 7 1.615],...
                    'Style','Edit','String','0','BackgroundColor','white',...
                    'Tag','edit_view_trend_ghc');
        uicontrol(p2_1,'units','characters','Position',[45.4 9.692 7 1.615],...
                    'Style','Edit','String','0','BackgroundColor','white',...
                    'Tag','edit_view_trend_ntol');
        uicontrol(p2_1,'units','characters','Position',[45.4 7.461 7 1.615],...
                    'Style','Edit','String','0','BackgroundColor','white',...
                    'Tag','edit_view_trend_atme');
        uicontrol(p2_1,'units','characters','Position',[45.4 5.230 7 1.615],...
                    'Style','Edit','String','0','BackgroundColor','white',...
                    'Tag','edit_view_trend_hydro');
                
        uicontrol(p2_1,'units','characters','Position',[56.4 14.385 5 1.077],...
                    'Style','CheckBox','Value',1,...
                    'Tag','check_view_demean_gravity');
        uicontrol(p2_1,'units','characters','Position',[56.4 12.154 5 1.077],...
                    'Style','CheckBox','Value',1,...
                    'Tag','check_view_demean_ghc');
        uicontrol(p2_1,'units','characters','Position',[56.4 9.923 5 1.077],...
                    'Style','CheckBox','Value',1,...
                    'Tag','check_view_demean_ntol');
        uicontrol(p2_1,'units','characters','Position',[56.4 7.692 5 1.077],...
                    'Style','CheckBox','Value',1,...
                    'Tag','check_view_demean_atme');
        uicontrol(p2_1,'units','characters','Position',[56.4 5.461 5 1.077],...
                    'Style','CheckBox','Value',1,...
                    'Tag','check_view_demean_hydro');
        % 2D data
        uicontrol(p2_2,'units','characters','Position',[3.6 2.769 10.4 1.154],...
                    'Style','Text','String','Load file');
        uicontrol(p2_2,'units','characters','Position',[16.8 2.769 15.6 1.154],...
                    'Style','Text','String','Exclude/Include');
        uicontrol(p2_2,'units','characters','Position',[47 2.769,10.4,1.154],...
                    'Style','Text','String','Max val.');
        uicontrol(p2_2,'units','characters','Position',[34.4 2.769,10.4,1.154],...
                    'Style','Text','String','File type');
        uicontrol(p2_2,'units','characters','Position',[57.4,2.769,13.4,1.154],...
                    'Style','Text','String','Map extend');
        uicontrol(p2_2,'units','characters','Position',[1.2 0.923,15,1.769],...
                    'Style','Pushbutton','String','grid',...
                    'Tag','push_view_file_spatial','UserData',[],...
                    'CallBack','mGlobe load_view_spatial');
        uicontrol(p2_2,'units','characters','Position',[90.2 1 12 1.769],...
                    'Style','Pushbutton','String','2D plot','FontWeight','bold',...
                    'Tag','push_view_2Dplot','UserData',[],...
                    'CallBack','mGlobe view_2Dplot');
        uicontrol(p2_2,'units','characters','Position',[34.6 1 11.2 1.692],...
                    'Style','Popupmenu','String','txt (Lon,Lat,data)|txt (Lat,Lon,data)|mat|arc ASCII|grd 6 text|netCDF',...
                    'BackgroundColor','white',...
                    'Tag','popup_view_type_spatial','Value',3);
        uicontrol(p2_2,'units','characters','Position',[58.8,1,11.2,1.692],...
                    'Style','Popupmenu','String','World|Europe|North America|South America|Africa|Asia|min-max|other',...
                    'BackgroundColor','white',...
                    'Tag','popup_view_mapextend','Value',1);
        uicontrol(p2_2,'units','characters','Position',[16.8,1,16.5,1.692],...
                    'Style','Popupmenu','String','All included|Exclude Greenland|Exclude Antarctica|Exclude Greenland+Antarctica|Load inlclusion polygon (*.txt)',...
                    'BackgroundColor','white',...
                    'Tag','popup_view_exclude','Value',1);
        uicontrol(p2_2,'units','characters','Position',[47.4,1,9.8,1.692],...
                    'Style','Edit','String','1000','BackgroundColor','white',...
                    'Tag','edit_view_data_max');
        uicontrol(p2_2,'units','characters','Position',[71.1 1 17 1.692],...
                    'Style','Popupmenu','String','Save As|fig|eps|tiff','BackgroundColor','white',...
                    'Tag','popup_view_printas_2D','Value',1,'UserData',[]);   

        % IN CONVERT
        % Gldas/Merra
        uicontrol(p3_1,'units','characters','Position',[1.8 1.11 14 1.08],...
                    'Style','Text','String','Model version');
        uicontrol(p3_1,'units','characters','Position',[62.2 1.108 18 1.08],...
                    'Style','Text','String','/GHM/NOAH025/',...
                    'Tag','text_down_gldas_path');
        uicontrol(p3_1,'units','characters','Position',[18.8 0.8 25 1.692],...
                    'Style','Popupmenu','String','GLDAS/CLM|GLDAS/MOS|GLDAS/NOAH025|GLDAS/NOAH10|GLDAS/VIC|MERRA Land|MERRA2 Land|GLDASv2.1/NOAH025',...
                    'Tag','popup_down_gldas_model','Value',3,'BackgroundColor','white',...
                    'CallBack','mGlobe select_down_model');
        uicontrol(p3_1,'units','characters','Position',[49.6 1.11 12.8 1.08],...
                    'Style','Text','String','Output path','UserData',fullfile('GHM','NOAH025'),...
                    'Tag','push_down_gldas_path');
        uicontrol(p3_1,'units','characters','Position',[81 0.88 21.4 1.7],...
                    'Style','Pushbutton','String','Convert','UserData',[],...
                    'Tag','push_down_gldas_run','CallBack','mGlobe calc_down_gldas');
        % Era/Ncel
        uicontrol(p3_2,'units','characters','Position',[1.2 0.6 8.4 1.692],...
                    'Style','Pushbutton','String','Input','UserData',[],...
                    'CallBack','mGlobe load_down_era_input','Tag',...
                    'push_down_era_input');
        uicontrol(p3_2,'units','characters','Position',[62.2 1.11 18 1],...
                    'Style','Text','String','/GHM/ERA|NCEP/',...
                    'Tag','text_down_era_path');
        uicontrol(p3_2,'units','characters','Position',[9.8,1.11,35.6,1.08],...
                    'Style','Text','String','input data *.nc',...
                    'Tag','text_down_era_input');
        uicontrol(p3_2,'units','characters','Position',[49.6 1.11 12.8 1],...
                    'Style','Text','String','Output path','UserData',fullfile('GHM','ERA'),...
                    'Tag','push_down_era_path');
        uicontrol(p3_2,'units','characters','Position',[81 0.88 10 1.7],...
                    'Style','Pushbutton','String','ERA',...
                    'CallBack','mGlobe calc_down_era');
        uicontrol(p3_2,'units','characters','Position',[92.4 0.88 10 1.7],...
                    'Style','Pushbutton','String','NCEP',...
                    'CallBack','mGlobe calc_down_ncep');
        % Other
        uicontrol(p3_3,'units','characters','Position',[1.2 0.6 8.4 1.692],...
                    'Style','Pushbutton','String','Input','UserData',[],...
                    'CallBack','mGlobe load_down_other_input','Tag',...
                    'push_down_other_input');
        uicontrol(p3_3,'units','characters','Position',[10.4 1.11 8 1.077],...
                    'Style','Text','String','Type:');
        uicontrol(p3_3,'units','characters','Position',[33 1.11 8.2 1.077],...
                    'Style','Text','String','Header:');
        uicontrol(p3_3,'units','characters','Position',[49.6 1.11 12.8 1.08],...
                    'Style','Text','String','Output path','UserData',fullfile('GHM','OTHER'),...
                    'Tag','push_down_other_path');
        uicontrol(p3_3,'units','characters','Position',[62.2 1.108 18 1.08],...
                    'Style','Text','String','/GHM/OTHER/',...
                    'Tag','text_down_other_path');
        uicontrol(p3_3,'units','characters','Position',[81 0.88 21.4 1.7],...
                    'Style','Pushbutton','String','Convert model',...
                    'CallBack','mGlobe calc_down_other');
        uicontrol(p3_3,'units','characters','Position',[18.2 0.8 13.6 1.692],...
                    'Style','Popupmenu','Value',1,'String',...
                    'Lat,Lon,Water (m)|Lat,Lon,Water (mm)|Lon,Lat,Water (m)|Lon,Lat,Water (mm)',...
                    'tag','popup_down_other_type','BackgroundColor','white');
        uicontrol(p3_3,'units','characters','Position',[41.6 0.8 7 1.692],...
                    'Style','Edit','String','0','BackgroundColor','white',...
                    'tag','edit_down_other_header');
        % time
        uicontrol(p3_4,'units','characters','Position',[1.2 4.846 7 1.077],...
                    'Style','Text','units','characters','String','Start');
        uicontrol(p3_4,'units','characters','Position',[0.8 2.923 7 1.077],...
                    'Style','Text','String','End');
        uicontrol(p3_4,'units','characters','Position',[1 0.923 7 1.077],...
                    'Style','Text','String','Step');
        uicontrol(p3_4,'units','characters','Position',[9.6,6.308,7,1.077],...
                    'Style','Text','String','year');
        uicontrol(p3_4,'units','characters','Position',[18.5,6.308,7,1.077],...
                    'Style','Text','String','month');
        uicontrol(p3_4,'units','characters','Position',[27.4,6.308,7,1.077],...
                    'Style','Text','String','day');
        uicontrol(p3_4,'units','characters','Position',[36.2,6.308,7,1.077],...
                    'Style','Text','String','hour');
        uicontrol(p3_4,'units','characters','Position',[10,4.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_start_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p3_4,'units','characters','Position',[18.8,4.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_start_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p3_4,'units','characters','Position',[27.6,4.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_start_day',...
                    'String','01','BackgroundColor','white')
        uicontrol(p3_4,'units','characters','Position',[36.4,4.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_start_hour',...
                    'String','12','BackgroundColor','white')
        uicontrol(p3_4,'units','characters','Position',[10,2.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_end_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p3_4,'units','characters','Position',[18.8,2.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_end_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p3_4,'units','characters','Position',[27.6,2.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_end_day',...
                    'String','03','BackgroundColor','white');
        uicontrol(p3_4,'units','characters','Position',[36.4,2.538,7,1.692],...
                    'Style','Edit','Tag','edit_down_time_end_hour',...
                    'String','12','BackgroundColor','white');
        uicontrol(p3_4,'units','characters','Position',[10,0.615,17.8,1.692],...
                    'Style','Popupmenu','Tag','popup_down_time_step',...
                    'String','3 hours|6 hours|12 hours|Day|2 days|Month','Value',6,'BackgroundColor','white');
        % GRACE
        uicontrol(p3_6,'units','characters','Position',[11.8 5.154 38.8 1.077],...
                    'Style','Text','units','characters','Tag','text_down_grace_input',...
                    'String','input (TELLUS) data *.nc','UserData',[]);
        uicontrol(p3_6,'units','characters','Position',[0.4 3.231 28.8 1],...
                    'Style','Text','Tag','text_down_grace_land',...
                    'String','Path+prefix:/GRACE/LAND/');
        uicontrol(p3_6,'units','characters','Position',[1.4 4.769 8.4 1.692],...
                    'Style','Pushbutton','String','Input','UserData',[],...
                    'CallBack','mGlobe load_down_grace_input','Tag',...
                    'push_down_grace_input')
        uicontrol(p3_6,'units','characters','Position',[29.4 2.923 21 1.692],...
                    'Style','Edit','Tag','edit_down_grace_name',...
                    'String','GRC_GFZ_RL05_CONv1409s','BackgroundColor','white');
        uicontrol(p3_6,'units','characters','Position',[1.4 0.846 15 1.692],...
                    'Style','Popupmenu','Tag','popup_down_grace_land',...
                    'String','Land Grid|Ocean Grid','Value',1,'BackgroundColor','white','CallBack','mGlobe select_grace_model');
        uicontrol(p3_6,'units','characters','Position',[17.4 0.924 11 1.692],...
                    'Style','Checkbox','Tag','check_down_grace_land',...
                    'String','Scale','Value',1);
        uicontrol(p3_6,'units','characters','Position',[29 0.846 21.4 1.692],...
                    'Style','Pushbutton','String','Convert GRACE',...
                    'CallBack','mGlobe calc_down_grace');
        % Dem convert
        uicontrol(p3_5,'units','characters','Position',[1.2 0.6 8.4 1.692],...
                    'Style','Pushbutton','String','Input',...
                    'Tag','push_down_dem_input','CallBack','mGlobe load_down_dem_input');
        uicontrol(p3_5,'units','characters','Position',[9.8 1.11 23.6 1.08],...
                    'Style','Text','String','input DEM file *.',...
                    'Tag','text_down_dem_input');
        uicontrol(p3_5,'units','characters','Position',[48 0.8 8.4 1.692],...
                    'Style','Pushbutton','String','Output',...
                    'Tag','push_down_dem_output','CallBack','mGlobe load_down_dem_output');
        uicontrol(p3_5,'units','characters','Position',[56.6 1.077 23.6 1.077],...
                    'Style','Text','String','output file *.mat',...
                    'Tag','text_down_dem_output');
        uicontrol(p3_5,'units','characters','Position',[81 0.88 21.4 1.7],...
                    'Style','Pushbutton','String','Convert DEM',...
                    'Tag','push_down_dem_run','CallBack','mGlobe calc_down_dem');
        uicontrol(p3_5,'units','characters','Position',[33.6 0.8 13 1.69],...
                    'Style','Popupmenu','String','txt (Lon,Lat,H)|txt (Lat,Lon,H)|arc ascii|grd 6 text|netCDF',...
                    'Tag','popup_down_dem_type','BackgroundColor','white');      
        % IN ATMO 
        % position
        uicontrol(p4_1,'units','characters','Position',[0.9,5.615,10,1.077],...
                    'Style','Text','String','Latitude');
        uicontrol(p4_1,'units','characters','Position',[1.7,3.385,10,1.077],...
                    'Style','Text','String','Longitude');
        uicontrol(p4_1,'units','characters','Position',[0.02,1.154,10,1.077],...
                    'Style','Text','String','Height');
        uicontrol(p4_1,'units','characters','Position',[36.8,5.692,5.6,1],...
                    'Style','Text','String','deg');
        uicontrol(p4_1,'units','characters','Position',[36.8,1.385,5.6,1],...
                    'Style','Text','String','m');
        uicontrol(p4_1,'units','characters','Position',[36.8,3.538,5.6,1],...
                    'Style','Text','String','deg');
        uicontrol(p4_1,'units','characters','Position',[13.6,5.308,23.2,1.692],...
                    'Style','Edit','Tag','edit_atmo_pos_lat',...
                    'String','48.24885','BackgroundColor','white');
        uicontrol(p4_1,'units','characters','Position',[13.6,3.077,23.2,1.692],...
                    'Style','Edit','Tag','edit_atmo_pos_lon',...
                    'String','16.35650','BackgroundColor','white')
        uicontrol(p4_1,'units','characters','Position',[13.6,0.846,23.2,1.692],...
                    'Style','Edit','Tag','edit_atmo_pos_hei',...
                    'String','192.70','BackgroundColor','white')
        % time
        uicontrol(p4_2,'units','characters','Position',[1.4,5.385,7,1.077],...
                    'Style','Text','units','characters','String','Start');
        uicontrol(p4_2,'units','characters','Position',[1,3.154,7,1.077],...
                    'Style','Text','String','End');
        uicontrol(p4_2,'units','characters','Position',[1.4,0.923,7,1.077],...
                    'Style','Text','String','Step');
        uicontrol(p4_2,'units','characters','Position',[10.8 7.054 5.6,1],...
                    'Style','Text','String','year');
        uicontrol(p4_2,'units','characters','Position',[18.8 7.054 5.9,1],...
                    'Style','Text','String','month');
        uicontrol(p4_2,'units','characters','Position',[28.4 7.054 5.6,1],...
                    'Style','Text','String','day');
        uicontrol(p4_2,'units','characters','Position',[37 7.054 5.6,1],...
                    'Style','Text','String','hour');
        uicontrol(p4_2,'units','characters','Position',[10.4 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_atmo_time_start_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p4_2,'units','characters','Position',[19.2 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_atmo_time_start_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p4_2,'units','characters','Position',[28 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_atmo_time_start_day',...
                    'String','01','BackgroundColor','white')
        uicontrol(p4_2,'units','characters','Position',[36.8 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_atmo_time_start_hour',...
                    'String','12','BackgroundColor','white')
        uicontrol(p4_2,'units','characters','Position',[10.4,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_atmo_time_end_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p4_2,'units','characters','Position',[19.2,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_atmo_time_end_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p4_2,'units','characters','Position',[28,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_atmo_time_end_day',...
                    'String','03','BackgroundColor','white');
        uicontrol(p4_2,'units','characters','Position',[36.8,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_atmo_time_end_hour',...
                    'String','12','BackgroundColor','white');
        uicontrol(p4_2,'units','characters','Position',[10.4,0.615,18,1.692],...
                    'Style','Popupmenu','Tag','popup_atmo_time_step',...
                    'String','6 hours|12 hours|Day|2 days','Value',2,'BackgroundColor','white');
        % Output
        uicontrol(p4_4,'units','characters','Position',[1.6 1.88 12 1.689],...
                    'Style','Pushbutton','units','characters',...
                    'String','File','UserData','output.txt',...
                    'Tag','push_atmo_out','CallBack','mGlobe load_atmo_out');
        uicontrol(p4_4,'units','characters','Position',[14 1.83 46.6 1.538],...
                    'Style','Text','Tag','text_atmo_out',...
                    'String','Chose your output file (output.txt default)');
        uicontrol(p4_4,'units','characters','Position',[1.8,0.215,8,1.692],...
                    'Style','Checkbox','Tag','check_atmo_xls',...
                    'String','xls','Value',0);
        uicontrol(p4_4,'units','characters','Position',[10.6,0.215,7.6,1.692],...
                    'Style','Checkbox','Tag','check_atmo_txt',...
                    'String','txt','Value',1);
        uicontrol(p4_4,'units','characters','Position',[19,0.215,7.6,1.692],...
                    'Style','Checkbox','Tag','check_atmo_tsf',...
                    'String','tsf','Value',0);
        uicontrol(p4_4,'units','characters','Position',[28 0.231 22.8 1.692],...
                    'Style','Checkbox','Tag','check_atmo_average',...
                    'String','subtract average','Value',0);
%         uicontrol(p4_4,'units','characters','Position',[46 0.231 13.5 1.692],...
%                     'Style','Pushbutton','units','characters',...
%                     'String','SHP (optional)','UserData',0,...
%                     'Tag','push_atmo_load_shp','CallBack','mGlobe push_atmo_shp');
        % ECMWF data
        uicontrol(p4_3,'units','characters','Position',[1.8 6.923 41.4 1.077],...
                    'Style','Text','String','Pressure levels (37 for ERA, MERRA 42)');
        uicontrol(p4_3,'units','characters','Position',[52.2 6.923 20.8 1.077],...
                    'Style','Text','String','Surface data');
        uicontrol(p4_3,'units','characters','Position',[21.2 5.231 30.8 1.077],...
                    'Style','Text','Tag','text_atmo_load_geopotential',...
                    'String','select geopotential height data');
        uicontrol(p4_3,'units','characters','Position',[21.2 3.077 30.8 1.077],...
                    'Style','Text','Tag','text_atmo_load_humidity',...
                    'String','select spec.humidity data');
        uicontrol(p4_3,'units','characters','Position',[21.2 0.923 30.8 1.077],...
                    'Style','Text','Tag','text_atmo_load_temperature',...
                    'String','select temperature data');
        uicontrol(p4_3,'units','characters','Position',[72.2 5.231 30.8 1.077],...
                    'Style','Text','Tag','text_atmo_load_surface',...
                    'String','select surface data file');
        uicontrol(p4_3,'units','characters','Position',[72.2 3.077 30.8 1.077],...
                    'Style','Text','Tag','text_atmo_load_orography',...
                    'String','select orography file');
        uicontrol(p4_3,'units','characters','Position',[72.2 0.692 30.8 1.077],...
                    'Style','Text','Tag','text_atmo_load_merra',...
                    'String','select surface pressure file','Visible','on','FontAngle','italic'); 
        uicontrol(p4_3,'units','characters','Position',[2 4.846 19 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Geopotential','UserData',[],...
                    'Tag','push_atmo_geopotential','CallBack',...
                    'mGlobe load_atmo_geopotential');
        uicontrol(p4_3,'units','characters','Position',[2 2.692 19 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Specific humidity','UserData',[],...
                    'Tag','push_atmo_humidity','CallBack',...
                    'mGlobe load_atmo_humidity');
        uicontrol(p4_3,'units','characters','Position',[2 0.538 19 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Temperature','UserData',[],...
                    'Tag','push_atmo_temperature','CallBack',...
                    'mGlobe load_atmo_temperature');
        uicontrol(p4_3,'units','characters','Position',[52.8 4.846 19 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Surface data','UserData',[],...
                    'Tag','push_atmo_surface','CallBack',...
                    'mGlobe load_atmo_surface');
        uicontrol(p4_3,'units','characters','Position',[52.8 2.692 19 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Orography','UserData',[],...
                    'Tag','push_atmo_orography','CallBack',...
                    'mGlobe load_atmo_orography');
        uicontrol(p4_3,'units','characters','Position',[52.8 0.538 19 1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','MERRA surf.p','UserData',[],...
                    'Tag','push_atmo_merra','FontAngle','italic','CallBack',...
                    'mGlobe load_atmo_merra','Visible','on'); 
                
        % IN OCEAN
        % position
        uicontrol(p5_1,'units','characters','Position',[0.9,5.615,10,1.077],...
                    'Style','Text','String','Latitude');
        uicontrol(p5_1,'units','characters','Position',[1.7,3.385,10,1.077],...
                    'Style','Text','String','Longitude');
        uicontrol(p5_1,'units','characters','Position',[0.02,1.154,10,1.077],...
                    'Style','Text','String','Height');
        uicontrol(p5_1,'units','characters','Position',[36.8,5.692,5.6,1],...
                    'Style','Text','String','deg');
        uicontrol(p5_1,'units','characters','Position',[36.8,1.385,5.6,1],...
                    'Style','Text','String','m');
        uicontrol(p5_1,'units','characters','Position',[36.8,3.538,5.6,1],...
                    'Style','Text','String','deg');
        uicontrol(p5_1,'units','characters','Position',[13.6,5.308,23.2,1.692],...
                    'Style','Edit','Tag','edit_ocean_pos_lat',...
                    'String','48.24885','BackgroundColor','white');
        uicontrol(p5_1,'units','characters','Position',[13.6,3.077,23.2,1.692],...
                    'Style','Edit','Tag','edit_ocean_pos_lon',...
                    'String','16.35650','BackgroundColor','white')
        uicontrol(p5_1,'units','characters','Position',[13.6,0.846,23.2,1.692],...
                    'Style','Edit','Tag','edit_ocean_pos_hei',...
                    'String','192.70','BackgroundColor','white')
        % time
        uicontrol(p5_2,'units','characters','Position',[1.4,5.385,7,1.077],...
                    'Style','Text','units','characters','String','Start');
        uicontrol(p5_2,'units','characters','Position',[1,3.154,7,1.077],...
                    'Style','Text','String','End');
        uicontrol(p5_2,'units','characters','Position',[1.4,0.923,7,1.077],...
                    'Style','Text','String','Step');
        uicontrol(p5_2,'units','characters','Position',[10.8 7.054 5.6,1],...
                    'Style','Text','String','year');
        uicontrol(p5_2,'units','characters','Position',[18.8 7.054 5.9,1],...
                    'Style','Text','String','month');
        uicontrol(p5_2,'units','characters','Position',[28.4 7.054 5.6,1],...
                    'Style','Text','String','day');
        uicontrol(p5_2,'units','characters','Position',[37 7.054 5.6,1],...
                    'Style','Text','String','hour');
        uicontrol(p5_2,'units','characters','Position',[10.4 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_ocean_time_start_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p5_2,'units','characters','Position',[19.2 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_ocean_time_start_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p5_2,'units','characters','Position',[28 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_ocean_time_start_day',...
                    'String','01','BackgroundColor','white')
        uicontrol(p5_2,'units','characters','Position',[36.8 5.077 7 1.692],...
                    'Style','Edit','Tag','edit_ocean_time_start_hour',...
                    'String','12','BackgroundColor','white')
        uicontrol(p5_2,'units','characters','Position',[10.4,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_ocean_time_end_year',...
                    'String','2002','BackgroundColor','white');
        uicontrol(p5_2,'units','characters','Position',[19.2,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_ocean_time_end_month',...
                    'String','11','BackgroundColor','white');
        uicontrol(p5_2,'units','characters','Position',[28,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_ocean_time_end_day',...
                    'String','03','BackgroundColor','white');
        uicontrol(p5_2,'units','characters','Position',[36.8,2.846,7,1.692],...
                    'Style','Edit','Tag','edit_ocean_time_end_hour',...
                    'String','12','BackgroundColor','white');
        uicontrol(p5_2,'units','characters','Position',[10.4,0.615,18,1.692],...
                    'Style','Popupmenu','Tag','popup_ocean_time_step',...
                    'String','3 hours|6 hours|12 hours|Day|2 days|Month','Value',6,'BackgroundColor','white');
        % Model Conversion
        uicontrol(p5_3,'units','characters','Position',[1.6,2.623+4,12,1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Input','UserData',[],...
                    'Tag','push_ocean_load_convert_input','CallBack',...
                    'mGlobe load_ocean_convert_input');
        uicontrol(p5_3,'units','characters','Position',[13.6,2.515+4,60,1.6],...
                    'Style','Text','Tag','text_ocean_convert_input',...
                    'String','Input file for conversion only (ECCO and OMCT)','UserData',[]);
        uicontrol(p5_3,'units','characters','Position',[60,0.238+3.8,8,1.769],...
                    'Style','Text','units','characters',...
                    'String','Paht:','UserData',fullfile('OBPM'),...
                    'Tag','push_ocean_convert_output');
        uicontrol(p5_3,'units','characters','Position',[69,0.238+3.8,28,1.769],...
                    'Style','Text','Tag','text_ocean_convert_output',...
                    'String','/OBPM/ECCO1/','UserData',fullfile('OBPM','ECCO1'));
        uicontrol(p5_3,'units','characters','Position',[1.45,1.9,15.25,1.769],...
                    'Style','Text',...
                    'String','Mass conserv.:');
        uicontrol(p5_3,'units','characters','Position',[78.9,2.623+4,24.1,1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Convert (ECCO/OMCT)','UserData',[],...
                    'Tag','push_ocean_convert','CallBack',...
                    'mGlobe load_ocean_convert_calc');
        uicontrol(p5_3,'units','characters','Position',[1.6,0.238+3.8,14.5,1.769],...
                    'Style','Text','Tag','text_ocean_convert_model',...
                    'String','Model version:');
        uicontrol(p5_3,'units','characters','Position',[18,0.238+4.17,40,1.769],...
                    'Style','Popupmenu','units','characters','BackgroundColor','white',...
                    'String','ECCO-JPL (conversion: *.txt for monthly and *.cdf for 12 hours)|Other (calculate only)|GRACE (calculate only)|ECCO2_beta (conversion: *.nc)|OMCT_oba (conversion: *.txt AOD1B spherical harmonics)|OMCT_ocn (conversion: *.txt AOD1B spherical harmonics)|OMCT_atm (conversion: *.txt AOD1B spherical harmonics)|OMCT6_oba (conversion: *.txt AOD1B spherical harmonics)|OMCT6_ocn (conversion: *.txt AOD1B spherical harmonics)|OMCT6_atm (conversion: *.txt AOD1B spherical harmonics)','UserData',[],...
                    'Tag','popup_ocean_load_convert_model','CallBack','mGlobe popup_ocean_model');
        uicontrol(p5_3,'units','characters','Position',[18,2.2,40,1.769],...
                    'Style','Popupmenu','units','characters','BackgroundColor','white',...
                    'String','Off|Subtract area average|Use mean pressure time series','UserData',[],'Value',2,...
                    'Tag','popup_ocean_load_conserv','CallBack','mGlobe popup_ocean_conservation');
        uicontrol(p5_3,'units','characters','Position',[1.6,0.32,12,1.769],...
                    'Style','Pushbutton','units','characters',...
                    'String','Select','UserData',[],...
                    'Tag','push_ocean_load_time_series','CallBack',...
                    'mGlobe load_ocean_time_series','Visible','off');
        uicontrol(p5_3,'units','characters','Position',[25,0.01,30,1.769],...
                    'Style','Text','units','characters',...
                    'String','Select input file (*.txt)','UserData',[],...
                    'Tag','text_ocean_load_time_series','Visible','off');
        uicontrol(p5_3,'units','characters','Position',[59,0.01,14,1.769],...
                    'Style','Text','units','characters',...
                    'String','Time column:','UserData',[],...
                    'Tag','text_ocean_load_time_series_time','Visible','off');
        uicontrol(p5_3,'units','characters','Position',[79,0.01,18,1.769],...
                    'Style','Text','units','characters',...
                    'String','Pressure column:','UserData',[],...
                    'Tag','text_ocean_load_time_series_press','Visible','off');
        uicontrol(p5_3,'units','characters','Position',[74,0.32,4,1.769],...
                    'Style','Edit','units','characters',...
                    'String','1','BackgroundColor','w',...
                    'Tag','edit_ocean_load_time_series_time','Visible','off');
        uicontrol(p5_3,'units','characters','Position',[79+18+1,0.32,4,1.769],...
                    'Style','Edit','units','characters',...
                    'String','9','BackgroundColor','w',...
                    'Tag','edit_ocean_load_time_series_press','Visible','off');
        % Output
        uicontrol(p5_4,'units','characters','Position',[1.6 1.88 12 1.689],...
                    'Style','Pushbutton','units','characters',...
                    'String','File','UserData','output.xls',...
                    'Tag','push_ocean_out','CallBack','mGlobe load_ocean_out');
        uicontrol(p5_4,'units','characters','Position',[14,1.83,46.6,1.538],...
                    'Style','Text','Tag','text_ocean_out',...
                    'String','Chose your output file (output.txt default)');
        uicontrol(p5_4,'units','characters','Position',[1.8,0.215,8,1.692],...
                    'Style','Checkbox','Tag','check_ocean_xls',...
                    'String','xls','Value',0);
        uicontrol(p5_4,'units','characters','Position',[10.6,0.215,7.6,1.692],...
                    'Style','Checkbox','Tag','check_ocean_txt',...
                    'String','txt','Value',1);
        uicontrol(p5_4,'units','characters','Position',[19,0.215,7.6,1.692],...
                    'Style','Checkbox','Tag','check_ocean_tsf',...
                    'String','tsf','Value',0);
        uicontrol(p5_4,'units','characters','Position',[28 0.231 22.8 1.692],...
                    'Style','Checkbox','Tag','check_ocean_average',...
                    'String','subtract average','Value',0);
        uicontrol(p5_4,'units','characters','Position',[46 0.231 13.5 1.692],...
                    'Style','Pushbutton','units','characters',...
                    'String','SHP (optional)','UserData',0,...
                    'Tag','push_ocean_load_shp','CallBack','mGlobe push_ocean_shp');
        % Set model paths
        mGlobe('select_hydro_model');
        mGlobe('select_down_model');
        mGlobe('popup_ocean_model');
        % END GUI
        end
%% Button response
        
    else
        switch in_switch                                                    % switch between callbacks  
            %% Panel Visibility
            case 'switch_hydro'                                             
                set(findobj('Tag','uipanel_main_hydr'),'Visible','on');     % show only required panel
                set(findobj('Tag','uipanel_main_atmo'),'Visible','off');
                set(findobj('Tag','uipanel_main_view'),'Visible','off');
                set(findobj('Tag','uipanel_main_down'),'Visible','off');
                set(findobj('Tag','uipanel_main_ocean'),'Visible','off');
                set(findobj('Tag','radio_switch_hydr'),'Value',1);          % switch radio button on/off
                set(findobj('Tag','radio_switch_atmo'),'Value',0);
                set(findobj('Tag','radio_switch_view'),'Value',0);
                set(findobj('Tag','radio_switch_down'),'Value',0);
                set(findobj('Tag','radio_switch_ocean'),'Value',0);
                set(findobj('Tag','text_status'),'String','Set your global hydrological effect');
            case 'switch_atmo'
                set(findobj('Tag','uipanel_main_hydr'),'Visible','off');
                set(findobj('Tag','uipanel_main_atmo'),'Visible','on');
                set(findobj('Tag','uipanel_main_view'),'Visible','off');
                set(findobj('Tag','uipanel_main_down'),'Visible','off');
                set(findobj('Tag','uipanel_main_ocean'),'Visible','off');
                set(findobj('Tag','radio_switch_hydr'),'Value',0);
                set(findobj('Tag','radio_switch_atmo'),'Value',1);
                set(findobj('Tag','radio_switch_view'),'Value',0);
                set(findobj('Tag','radio_switch_down'),'Value',0);
                set(findobj('Tag','radio_switch_ocean'),'Value',0);
                set(findobj('Tag','text_status'),'String','Set your global atmospheric effect');
            case 'switch_view'
                set(findobj('Tag','uipanel_main_hydr'),'Visible','off');
                set(findobj('Tag','uipanel_main_atmo'),'Visible','off');
                set(findobj('Tag','uipanel_main_view'),'Visible','on');
                set(findobj('Tag','uipanel_main_down'),'Visible','off');
                set(findobj('Tag','uipanel_main_ocean'),'Visible','off');
                set(findobj('Tag','radio_switch_hydr'),'Value',0);
                set(findobj('Tag','radio_switch_atmo'),'Value',0);
                set(findobj('Tag','radio_switch_view'),'Value',1);
                set(findobj('Tag','radio_switch_down'),'Value',0);
                set(findobj('Tag','radio_switch_ocean'),'Value',0);
                set(findobj('Tag','text_status'),'String','Set your visualisation options');
            case 'switch_down'
                set(findobj('Tag','uipanel_main_hydr'),'Visible','off');
                set(findobj('Tag','uipanel_main_atmo'),'Visible','off');
                set(findobj('Tag','uipanel_main_view'),'Visible','off');
                set(findobj('Tag','uipanel_main_down'),'Visible','on');
                set(findobj('Tag','uipanel_main_ocean'),'Visible','off');
                set(findobj('Tag','radio_switch_hydr'),'Value',0);
                set(findobj('Tag','radio_switch_atmo'),'Value',0);
                set(findobj('Tag','radio_switch_view'),'Value',0);
                set(findobj('Tag','radio_switch_down'),'Value',1);
                set(findobj('Tag','radio_switch_ocean'),'Value',0);
                set(findobj('Tag','text_status'),'String','Set your conversion options');
            case 'switch_ocean'
                set(findobj('Tag','uipanel_main_hydr'),'Visible','off');
                set(findobj('Tag','uipanel_main_atmo'),'Visible','off');
                set(findobj('Tag','uipanel_main_view'),'Visible','off');
                set(findobj('Tag','uipanel_main_down'),'Visible','off');
                set(findobj('Tag','uipanel_main_ocean'),'Visible','on');
                set(findobj('Tag','radio_switch_hydr'),'Value',0);
                set(findobj('Tag','radio_switch_atmo'),'Value',0);
                set(findobj('Tag','radio_switch_view'),'Value',0);
                set(findobj('Tag','radio_switch_down'),'Value',0);
                set(findobj('Tag','radio_switch_ocean'),'Value',1);
                set(findobj('Tag','text_status'),'String','Set your non-tidal ocean effect');
            %% HYROLOGICAL EFFECT
            case 'select_hydro_model'                                       % change the shown path according to selected model
                % Read path set file
                [ghm_main,~,grace_main] = mGlobe_getModelPath;
                % %% I have created GHM/ERA5 folder. this will just concatenate base path + ERA5 suffix
                ghm_path = {fullfile(ghm_main,'CLM'),fullfile(ghm_main,'MOS'),...
                            fullfile(ghm_main,'NOAH025'),fullfile(ghm_main,'NOAH10'),...
                            fullfile(ghm_main,'VIC'),fullfile(ghm_main,'ERA'),fullfile(ghm_main,'MERRA'),...
                            fullfile(ghm_main,'OTHER'),fullfile(grace_main,'LAND'),fullfile(ghm_main,'NCEP'),...
                            fullfile(ghm_main,'MERRA2'),fullfile(ghm_main,'NCEP'),...
                            fullfile(ghm_main,'NOAH025v21'),fullfile(ghm_main,'ERA5')};
                val = get(findobj('Tag','popup_hydro_model'),'Value');
                set(findobj('Tag','push_hydro_model_path'),'UserData',ghm_path{val}); % each button stores data about the path to model data
                set(findobj('Tag','text_hydro_model_path'),'String',ghm_path{val}); % show path for current model
                set(findobj('Tag','popup_hydro_model_layer'),'Value',1); % always set the layer to 'total' (by default)
                set(findobj('Tag','push_down_other_path'),'UserData',fullfile(ghm_main,'OTHER'));
                set(findobj('Tag','text_down_other_path'),'String',fullfile(ghm_main,'OTHER'));
                % %% Will need to adjust for ERA5 (need to have Matlab to see what it does)
                if val == 6
                    set(findobj('Tag','push_down_era_path'),'UserData',fullfile(ghm_main,'ERA'));
                    set(findobj('Tag','text_down_era_path'),'String',fullfile(ghm_main,'ERA'));
                else
                    set(findobj('Tag','push_down_era_path'),'UserData',fullfile(ghm_main,'NCEP'));
                    set(findobj('Tag','text_down_era_path'),'String',fullfile(ghm_main,'NCEP'));
                end
                switch val
                    case 1
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilm1|soilm2|soilm3|soilm4|soilm5|soilm6|soilm7|soilm8|soilm9|soilm10|swe'); % show different layers for different models
                    case 2
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilm1|soilm2|soilm3|swe');
                    case 3
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilm1|soilm2|soilm3|soilm4|swe');
                    case 4
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilm1|soilm2|soilm3|soilm4|swe');
                    case 5
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilm1|soilm2|soilm3|swe');
                    case 6
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|swvl1|swvl2|swvl3|swvl4|sd');
                    case 7
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total (twland)');
                    case 8
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total');
                    case 9
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total');
                    case 10
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilw1|soilw2|weasd');
                    case 11
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total (twland)');
                    case 12
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilw1|soilw2|weasd');
                    case 13
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|soilm1|soilm2|soilm3|soilm4|swe');
                    case 14
                        % %% This should be now OK = 4 SM layers + snow: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation
                        set(findobj('Tag','popup_hydro_model_layer'),'String','total|swvl1|swvl2|swvl3|swvl4|sd');
                end
            case 'load_hydro_dem'                                          % Load DEM for Hydro effect (Continental water storage)
                [name,path] = uigetfile('*.*','Load DEM up to 1 deg from point of observation');
                if name == 0                                               % If cancelled-> no DEM used
                    set(findobj('Tag','push_hydro_topo_dem'),'UserData',[]); % empty if no data selected
                    set(findobj('Tag','text_hydro_topo_dem'),'String','No DEM used');
                else
                    set(findobj('Tag','push_hydro_topo_dem'),'UserData',[path,name]); % store the full file name
                    set(findobj('Tag','text_hydro_topo_dem'),'String',name); % display the file name (without path)
                end
            case 'push_hydro_shp'                                          % Load Coastline shapefile
                [name,path] = uigetfile('*.*','Load SHP for points up to 1 deg from point of observation');
                if name == 0                                               % If cancelled-> no Shapefile
                    set(findobj('Tag','push_hydro_load_shp'),'UserData',0); % empty if no data selected
                else
                    set(findobj('Tag','push_hydro_load_shp'),'UserData',[path,name]); % store the full file name
                end
            case 'push_ocean_shp'                                          % Load Coastline shapefile
                [name,path] = uigetfile('*.*','Load SHP for points up to 1 deg from point of observation');
                if name == 0                                               % If cancelled-> no Shapefile
                    set(findobj('Tag','push_ocean_load_shp'),'UserData',0); % empty if no data selected
                else
                    set(findobj('Tag','push_ocean_load_shp'),'UserData',[path,name]); % store the full file name
                end
            case 'push_atmo_shp'                                          % Load Coastline shapefile
                [name,path] = uigetfile('*.*','Load SHP for points up to 1 deg from point of observation');
                if name == 0                                               % If cancelled-> no Shapefile
                    set(findobj('Tag','push_atmo_load_shp'),'UserData',0); % empty if no data selected
                else
                    set(findobj('Tag','push_atmo_load_shp'),'UserData',[path,name]); % store the full file name
                end
            case 'load_hydro_show'                                          % Show loaded DEM + point of observation
                dem_path = get(findobj('Tag','push_hydro_topo_dem'),'UserData'); % first, get full file name
                if isempty(dem_path)                                        % check if a DEM has been selected
                    set(findobj('Tag','text_status'),'String','Load your DEM file first');
                else
                    try                                                     % try to display the DEM
                        set(findobj('Tag','text_status'),'String','Loading DEM ...');drawnow % send message to GUI status line
                        dem = importdata(dem_path);                             % load DEM
                        Input = [str2double(get(findobj('Tag','edit_hydro_pos_lat'),'String')),... % Get input coordinates
                        str2double(get(findobj('Tag','edit_hydro_pos_lon'),'String')),...
                        str2double(get(findobj('Tag','edit_hydro_pos_hei'),'String'))];
                        figure;                                                 % open new figure
                        worldmap([min(min(dem.lat)),max(max(dem.lat))],...      % create optimized (min-max) map projection
                             [min(min(dem.lon)),max(max(dem.lon))]);
                        surfm(dem.lat,dem.lon,dem.height,dem.height);           % display the DEM
                        land = shaperead('landareas', 'UseGeoCoords', true);    % load continents
                        plot3m([land.Lat],[land.Lon],max(max(dem.height))+100,'Color','black'); % plot the continents boarders
                        plot3m(Input(1),Input(2),Input(3)+100,'k.','MarkerSize',18); % add the point of observation/calculation
                        color_id = colorbar('location','southoutside');         % show colorbar
                        set(get(color_id,'xlabel'),'string','Altitude (m)');clear color_id % add colorbar name
                    catch
                        set(findobj('Tag','text_status'),'String','Load valid DEM file first (including dem.lon,dem.lat,dem.height)');drawnow % send warning if the DEM has not been loaded/displayed
                    end
                    set(findobj('Tag','text_status'),'String','Set your global hydrological effect');drawnow % confirmation
                end
            case 'load_hydro_out'                                           % Set output file = Hydro panel
                [name,path] = uiputfile('*.*','Output file: excel sheet or txt');
                if name == 0                                                % If cancelled-> default file (output.txt)
                    set(findobj('Tag','push_hydro_out'),'UserData','output.txt'); % store default file name if no file has been selected
                    set(findobj('Tag','text_hydro_out'),'String','output.txt');
                else
                    set(findobj('Tag','push_hydro_out'),'UserData',[path,name]); % store the selected full file name
                    set(findobj('Tag','text_hydro_out'),'String',name);
                end
            case 'load_hydro_include'                                       % load inclusion polygon (grid cells out of this polygon will be set to zero) = Hydro panel
                [name,path] = uigetfile('*.txt','Select your inclusion polygon *.txt [Lon (deg) Lat (deg)]');
                if name == 0                                                % If cancelled-> no inclusion polygon
                    set(findobj('Tag','push_hydro_include'),'UserData',[]);
                else
                    set(findobj('Tag','push_hydro_include'),'UserData',[path,name]);
                end
            case 'hydro_calc'                                              % Calculate the Continental Water Storage Effect
                Input = [str2double(get(findobj('Tag','edit_hydro_pos_lat'),'String')),... % Get input coordinates
                    str2double(get(findobj('Tag','edit_hydro_pos_lon'),'String')),...
                    str2double(get(findobj('Tag','edit_hydro_pos_hei'),'String'))];
                output_file = get(findobj('Tag','push_hydro_out'),'UserData'); % Get output file
                output_file_type = [get(findobj('Tag','check_hydro_xls'),'Value'),...
                                    get(findobj('Tag','check_hydro_txt'),'Value'),...
                                    get(findobj('Tag','check_hydro_tsf'),'Value')];
                subtract_average = get(findobj('Tag','check_hydro_average'),'Value');
                
                step_calc = get(findobj('Tag','popup_hydro_time_step'),'Value'); % Get time step
                if output_file_type(3) == 1  && step_calc == 6              % uncheck the TSF output if monthly data are selected (TSF requires constant sampling)
                    set(findobj('Tag','check_hydro_tsf'),'Value',0)
                    output_file_type(3) =  0;                               % turn of the TSF output
                end
                
                if sum(output_file_type) == 0                               % set txt output file if no output file type is selected by user 
                    set(findobj('Tag','check_hydro_txt'),'Value',1);
                    output_file_type(2) = 1;
                    drawnow                                                 % update the GUI
                end
                
                DEM_file = get(findobj('Tag','push_hydro_topo_dem'),'UserData'); % Get DEM file (full file name, empty if not selected)
                INCLUDE_file = get(findobj('Tag','push_hydro_include'),'UserData'); % Get inclusion polygon (full file name, empty if not selected)
                
                
                start_calc = datenum(str2double(get(findobj('Tag','edit_hydro_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_hydro_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_hydro_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_hydro_time_start_hour'),'String')),0,0);
                
                end_calc = datenum(str2double(get(findobj('Tag','edit_hydro_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_hydro_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_hydro_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_hydro_time_end_hour'),'String')),0,0);
                
                
                exclude_calc = [get(findobj('Tag','check_hydro_model_green'),'Value'),... % Get the exclusion areas
                    get(findobj('Tag','check_hydro_model_perma'),'Value')];
                
                model_calc = get(findobj('Tag','popup_hydro_model'),'Value'); % Get model name
                
                model_layer = get(findobj('Tag','popup_hydro_model_layer'),'Value');  % Get model layer
                
                mass_conserv = get(findobj('Tag','check_hydro_model_ocean'),'Value'); % Get mass conservation condition
                
                ghc_treshold = str2double(get(findobj('Tag','edit_hydro_treshold'),'String')); % Get calculation (distance) threshold
                ghc_path = get(findobj('Tag','push_hydro_model_path'),'UserData');
                name = 1;
                shp_file = get(findobj('Tag','push_hydro_load_shp'),'UserData');
                
                if model_calc == 8                                          % prompt user to pick OTHER model (fixed prefix)
                    set(findobj('Tag','text_status'),'String','You have chosen OTHER model => pick file with *.mat grid, first TEN letters will be used as a PREFIX (e.g. WGHM05_All)');
                    [name,~] = uigetfile(ghc_path,'Pick file with *.mat grid, first TEN letters will be used as a PREFIX123 for data loading (e.g. WGHM05_All)');
                    if name ~=0                                             % some file name is expected
                        ghc_path = fullfile(ghc_path,name(1:10));           % create new input file prefix = path + file prefix
                        set(findobj('Tag','text_status'),'String',['File path + prefix = ',ghc_path]); % send to status line
                    else
                        name = 0;                                           % the 'name' variable will be use to check if everything was selected correctly
                    end
                end
                
                if model_calc == 9                                          % prompt user to pick GRACE model (prefix will be determined)
                    set(findobj('Tag','text_status'),'String','You have chosen GRACE model => pick file with *.mat grid, first 22 letters will be used as a PREFIX (e.g. GRC_GFZ_RL05_FILT0_CON))');
                    [name,~] = uigetfile(fullfile(ghc_path),'Pick file with *.mat GRACE grid, first 22 letters will be used as a PREFIX (e.g. GRC_GFZ_RL05_FILT0_CON)');
                    if name ~=0                                             % some file name is expected
                        ghc_path = fullfile(ghc_path,name(1:22));           % create new input file prefix = path + file prefix
                        set(findobj('Tag','text_status'),'String',['File path + prefix = ',ghc_path]);
                    else
                        name = 0;
                    end
                end
                
                if mass_conserv == 3 && sum(exclude_calc) ~= 0              % change settings/exclusion area if coupled model is selected
                    exclude_calc(:,:) = 0;
                    set(findobj('Tag','text_status'),'String',...
                    'Warning: if coupled model (continent + ocean) is choosen, no areas will be excluded'); 
                    set(findobj('Tag','check_hydro_model_green'),'Value',0);
                    set(findobj('Tag','check_hydro_model_perma'),'Value',0);
                    drawnow 
                    pause(5);                                               % wait... 
                end
                if (ghc_treshold > 1 || ghc_treshold < 0.05) || (start_calc > end_calc) || (sum(double(name)) == 0) % warn user that the threshold must be within 0.05-1 degree and the starting date must be < than the end date
                    set(findobj('Tag','text_status'),'String','Threshold must be within <0.05,1.00> degree, start time <= end time and model prefix set correctly');
                elseif exist('mGlobe_DATA_dgE_Hydro.txt','file')==2 && exist('mGlobe_DATA_OceanGrid.mat','file')==2
                    set(findobj('Tag','text_status'),'String','Hydro: starting the computation...');drawnow
                    mGlobe_calc_Hydro(Input,output_file,output_file_type,DEM_file,...
                                        start_calc,end_calc,step_calc,exclude_calc,...
                                        model_calc,model_layer,mass_conserv,...
                                        ghc_treshold,ghc_path,subtract_average,...
                                        INCLUDE_file,shp_file); % calculate the hydrological effect
                    pause(5)                                               % wait 5 sec, than write message
                    set(findobj('Tag','text_status'),'String','Set the global hydrological effect');
                else
                    set(findobj('Tag','text_status'),'String',...           % warn user that these two files must be in the current folder
                    'Please ensure that both mGlobe_DATA_dgE_Hydro.txt and mGlobe_DATA_OceanGrid.mat file are on your MATLAB current folder'); drawnow 
                end
            %% VIEW RESULTS
            case 'load_view_gravity'                                        % Load gravity variation for visualisation
                [name,path] = uigetfile('*.*','Load observed gravity variation');
                if name == 0                                                % If cancelled-> no data
                    set(findobj('Tag','push_view_file_gravity'),'UserData',[]);
                    set(findobj('Tag','check_view_gravity'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_gravity'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_gravity'),'Value',1);
                end
            case 'load_view_ghc'                                            % Load gravity time series 1
                [name,path] = uigetfile('*.*','Load time series 1');
                if name == 0                                                % If cancelled-> no data
                    set(findobj('Tag','push_view_file_ghc'),'UserData',[]);
                    set(findobj('Tag','check_view_ghc'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_ghc'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_ghc'),'Value',1);
                end
            case 'load_view_ntol'                                           % Load gravity time series 2
                [name,path] = uigetfile('*.*','Load time series 2');
                if name == 0                                                % If cancelled-> no data
                    set(findobj('Tag','push_view_file_ntol'),'UserData',[]);
                    set(findobj('Tag','check_view_ntol'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_ntol'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_ntol'),'Value',1);
                end
            case 'load_view_atme'                                           % Load gravity time series 3
                [name,path] = uigetfile('*.*','Load time series 3');
                if name == 0                                                % If cancelled -> no data
                    set(findobj('Tag','push_view_file_atme'),'UserData',[]);
                    set(findobj('Tag','check_view_atme'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_atme'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_atme'),'Value',1);
                end
            case 'load_view_hydro'                                          % Load hydro. observations for visualisation
                [name,path] = uigetfile('*.*','Load observed hydrological parameter');
                if name == 0                                                % If cancelled-> no data
                    set(findobj('Tag','push_view_file_hydro'),'UserData',[]);
                    set(findobj('Tag','check_view_hydro'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_hydro'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_hydro'),'Value',1);
                end
            case 'load_view_rain'                                           % Load rain data for visualisation
                [name,path] = uigetfile('*.*','Load rain/snow data');
                if name == 0                                                % If cancelled-> no data
                    set(findobj('Tag','push_view_file_rain'),'UserData',[]);
                    set(findobj('Tag','check_view_rain'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_rain'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_rain'),'Value',1);
                end
            case 'load_view_spatial'                                       % Load 2D data
                [name,path] = uigetfile('*.*','Load 2D data (grid)');
                if name == 0                                               % If cancelled-> no data
                    set(findobj('Tag','push_view_file_spatial'),'UserData',[]);
                    set(findobj('Tag','check_view_spatial'),'Value',0);
                else
                    set(findobj('Tag','push_view_file_spatial'),'UserData',[path,name]);
                    set(findobj('Tag','check_view_spatial'),'Value',1);
                end
            case 'view_1Dplot'                                              % 1D PLOT
                plot_file{1} = get(findobj('Tag','push_view_file_gravity'),'UserData'); % get all plotting options
                plot_file{2} = get(findobj('Tag','push_view_file_ghc'),'UserData');
                plot_file{3} = get(findobj('Tag','push_view_file_ntol'),'UserData');
                plot_file{4} = get(findobj('Tag','push_view_file_atme'),'UserData');
                plot_file{5} = get(findobj('Tag','push_view_file_hydro'),'UserData');
                plot_file{6} = get(findobj('Tag','push_view_file_rain'),'UserData');
                plot_file{7} = get(findobj('Tag','push_view_file_rain'),'UserData');
                plot_option(1) = get(findobj('Tag','check_view_gravity'),'Value');
                plot_option(2) = get(findobj('Tag','check_view_ghc'),'Value');
                plot_option(3) = get(findobj('Tag','check_view_ntol'),'Value');
                plot_option(4) = get(findobj('Tag','check_view_atme'),'Value');
                plot_option(5) = get(findobj('Tag','check_view_hydro'),'Value');
                plot_option(6) = get(findobj('Tag','check_view_rain'),'Value');
                plot_option(7) = plot_option(6);
                file_type(1) = get(findobj('Tag','popup_view_type_gravity'),'Value');
                file_type(2) = get(findobj('Tag','popup_view_type_ghc'),'Value');
                file_type(3) = get(findobj('Tag','popup_view_type_ntol'),'Value');
                file_type(4) = get(findobj('Tag','popup_view_type_atme'),'Value');
                file_type(5) = get(findobj('Tag','popup_view_type_hydro'),'Value');
                file_type(6) = get(findobj('Tag','popup_view_type_rain'),'Value');
                file_type(7) = file_type(6);
                trend_type(1) = str2double(get(findobj('Tag','edit_view_trend_gravity'),'String'));
                trend_type(2) = str2double(get(findobj('Tag','edit_view_trend_ghc'),'String'));
                trend_type(3) = str2double(get(findobj('Tag','edit_view_trend_ntol'),'String'));
                trend_type(4) = str2double(get(findobj('Tag','edit_view_trend_atme'),'String'));
                trend_type(5) = str2double(get(findobj('Tag','edit_view_trend_hydro'),'String'));
                trend_type(6:7) = 0;
                demean_type(1) = get(findobj('Tag','check_view_demean_gravity'),'Value');
                demean_type(2) = get(findobj('Tag','check_view_demean_ghc'),'Value');
                demean_type(3) = get(findobj('Tag','check_view_demean_ntol'),'Value');
                demean_type(4) = get(findobj('Tag','check_view_demean_atme'),'Value');
                demean_type(5) = get(findobj('Tag','check_view_demean_hydro'),'Value');
                demean_type(6:7) = 0;
                time_column(1) = str2double(get(findobj('Tag','edit_view_time_gravity'),'String'));
                time_column(2) = str2double(get(findobj('Tag','edit_view_time_ghc'),'String'));
                time_column(3) = str2double(get(findobj('Tag','edit_view_time_ntol'),'String'));
                time_column(4) = str2double(get(findobj('Tag','edit_view_time_atme'),'String'));
                time_column(5) = str2double(get(findobj('Tag','edit_view_time_hydro'),'String'));
                time_column(6) = str2double(get(findobj('Tag','edit_view_time_rain'),'String'));
                time_column(7) = time_column(6);
                data_column(1) = str2double(get(findobj('Tag','edit_view_data_gravity'),'String'));
                data_column(2) = str2double(get(findobj('Tag','edit_view_data_ghc'),'String'));
                data_column(3) = str2double(get(findobj('Tag','edit_view_data_ntol'),'String'));
                data_column(4) = str2double(get(findobj('Tag','edit_view_data_atme'),'String'));
                data_column(5) = str2double(get(findobj('Tag','edit_view_data_hydro'),'String'));
                data_column(6) = str2double(get(findobj('Tag','edit_view_data_rain'),'String'));
                data_column(7) = str2double(get(findobj('Tag','edit_view_data_snow'),'String'));
                corr_option(1) = get(findobj('Tag','popup_view_add_ghc'),'Value');
                corr_option(2) = get(findobj('Tag','popup_view_add_ntol'),'Value');
                corr_option(3) = get(findobj('Tag','popup_view_add_atme'),'Value');
                save_option = get(findobj('Tag','popup_view_printas_1D'),'Value');
                export_option = get(findobj('Tag','popup_view_exportas_1D'),'Value');
                if save_option ~= 1                                         % prompt user to select an output file for printing
                    set(findobj('Tag','text_status'),'String',...
                    'Please select output file for printing, otherwise result will be printed in default file: output_print.*');
                    [name,path] = uiputfile('*.*','Choose printing output file (default output_print.*)');
                    if name == 0                                            % If cancelled-> print to default file
                        switch save_option
                            case 2
                                print_file = 'output_print.fig';
                            case 3
                                print_file = 'output_print.eps';
                            case 4
                                print_file = 'output_print.tiff';
                        end
                    else
                        print_file = [path,name];
                    end
                else
                    print_file = [];
                end
                if export_option ~= 1                                       % prompt user to select an output file plotting results
                    set(findobj('Tag','text_status'),'String',...
                    'Please select output file for results export, otherwise result will be written in default file: output_export.*');
                    [name,path] = uiputfile('*.*','Choose writing output file (default output_export.*)');
                    if name == 0                                            % If cancelled-> export to default file
                        switch export_option
                            case 2
                                export_file = 'output_export.txt';
                            case 3
                                export_file = 'output_print.mat';
                            case 4
                                export_file = 'output_print.xls';
                            case 5
                                export_file = 'output_print.tsf';
                        end
                    else
                        export_file = [path,name];
                    end
                else
                    export_file = [];
                end
                plot_file_control(1:7) = 0;                                 % check selected and loaded files                   
                id_row = find(plot_option == 1);
                if ~isempty(id_row)
                    for i = 1:length(id_row)
                        if isempty(plot_file{id_row(i)})
                            plot_file_control(id_row(i)) = 0;
                        else
                            plot_file_control(id_row(i)) = 1;
                        end
                    end
                end
                if sum(plot_option) == 0 || (sum(plot_file_control) ~= length(id_row)) % warn user that at leas one file must be selected
                    set(findobj('Tag','text_status'),'String',...
                    'Please select at least one dataset for plotting and load file for each checked dataset');
                else
                    check_out = mGlobe_view1D(plot_file,plot_option,file_type,time_column,data_column,corr_option,save_option,export_option,print_file,export_file,demean_type,trend_type); % perform 1D plot
                    if check_out == 1
                        set(findobj('Tag','text_status'),'String','Set your visualisation options');
                    else
                        set(findobj('Tag','text_status'),'String',...       % warn user that all files must be in valid format
                        'Please ensure that all input files are in correct format (file type and columns)');
                    end
                end  
            case 'view_2Dplot'                                              % 2D plot (surfm)
                plot_file = get(findobj('Tag','push_view_file_spatial'),'UserData'); % get all required data
                plot_exclude = get(findobj('Tag','popup_view_exclude'),'Value');
                file_type = get(findobj('Tag','popup_view_type_spatial'),'Value');
                max_value = str2double(get(findobj('Tag','edit_view_data_max'),'String'));
                map_extend = get(findobj('Tag','popup_view_mapextend'),'Value');
                save_option = get(findobj('Tag','popup_view_printas_2D'),'Value');
                if save_option ~= 1                                         % prompt user to select an output file for printing
                    set(findobj('Tag','text_status'),'String',...
                    'Please select output file for printing, otherwise result will be printed in default file: output_print.*');
                    [name,path] = uiputfile('*.*','Choose printing output file (default output_print.*)');
                    if name == 0                                            % If cancelled-> result will be printed to default file
                        switch save_option
                            case 2
                                print_file = 'output_print.fig';
                            case 3
                                print_file = 'output_print.eps';
                            case 4
                                print_file = 'output_print.tiff';
                        end
                    else
                        print_file = [path,name];
                    end
                else
                    print_file = [];
                end
                if isempty(plot_file)                                       % warn user to load input grid
                    set(findobj('Tag','text_status'),'String',...
                    'Please load input file');
                else
                    set(findobj('Tag','text_status'),'String','Plot: Visualizing...');drawnow
                    check_out = mGlobe_view2D(plot_file,file_type,max_value,map_extend,save_option,print_file,plot_exclude); % perform 2D visualisation
                    if check_out == 1                                       % print output message
                        set(findobj('Tag','text_status'),'String',...      
                        'Plot: Visualizing...');pause(3);
                        set(findobj('Tag','text_status'),'String','Set your visualisation options');
                    else
                        set(findobj('Tag','text_status'),'String',...
                        'Please ensure that all input files are in correct format (file type)');
                    end
                end
                
            %%  CONVERT
            % GLDAS
            case 'select_down_model'
                [ghm_main,~,~] = mGlobe_getModelPath;
                val = get(findobj('Tag','popup_down_gldas_model'),'Value'); % get GLDAS model version
                switch val                                                  % switch output path according to the selected model version
                    case 1
                    set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'CLM'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'CLM'));
                    case 2
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'MOS'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'MOS'));
                    case 3
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'NOAH025'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'NOAH025'));
                    case 4
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'NOAH10'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'NOAH10'));
                    case 5
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'VIC'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'VIC'));
                    case 6
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'MERRA'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'MERRA'));
                    case 7
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'MERRA2'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'MERRA2'));
                    case 8
                        set(findobj('Tag','push_down_gldas_path'),'UserData',fullfile(ghm_main,'NOAH025v21'));
                        set(findobj('Tag','text_down_gldas_path'),'String',fullfile(ghm_main,'NOAH025v21'));
                end
            case 'calc_down_gldas'
                [name,path] = uigetfile({'*.nc;*.nc4','NetCDF (*.nc,*.nc4)'},...
                            'Choose ONE input GLDAS or MERRA netCDF file'); % select MERRA/GLDAS interim netcdf (*.nc) input file
                if name == 0                                                
                    set(findobj('Tag','text_status'),'String','You must select one netcdf file');
                else
                    input_path = path;
                    input_file = name;
                    model_version = get(findobj('Tag','popup_down_gldas_model'),'Value'); % get all required settings
                    output_path = get(findobj('Tag','push_down_gldas_path'),'UserData');
                    start_calc = datenum(str2double(get(findobj('Tag','edit_down_time_start_year'),'String')),... % Get date of start
                        str2double(get(findobj('Tag','edit_down_time_start_month'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_start_day'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_start_hour'),'String')),0,0);
                    end_calc = datenum(str2double(get(findobj('Tag','edit_down_time_end_year'),'String')),... % Get date of end
                        str2double(get(findobj('Tag','edit_down_time_end_month'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_end_day'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_end_hour'),'String')),0,0);
                    step_calc = get(findobj('Tag','popup_down_time_step'),'Value'); % Get time step
                    if start_calc > end_calc                                    % check valid settings and perform conversion
                        set(findobj('Tag','text_status'),'String','Start time must be <= End time');
                    else
                        set(findobj('Tag','text_status'),'String','Models: Starting the conversion...');drawnow
                        if model_version <=8 
                            mGlobe_convert_GLDAS(start_calc,end_calc,model_version,step_calc,output_path,input_path,input_file);
                        end
                        set(findobj('Tag','text_status'),'String','Conversion completed');
                        pause(5);
                        set(findobj('Tag','text_status'),'String','Set your conversion options');
                    end
                end
            % ERA    
            case 'load_down_era_input'                                     
                [name,path] = uigetfile('*.nc','Choose your input ERA/NCEP netCDF file'); % select ERA interim netcdf (*.nc) input file
                if name == 0                                                
                    set(findobj('Tag','push_down_era_input'),'UserData',[]);
                    set(findobj('Tag','text_down_era_input'),'String',[]);
                else
                    set(findobj('Tag','push_down_era_input'),'UserData',[path,name]);
                    set(findobj('Tag','text_down_era_input'),'String',name);
                end 
            case 'calc_down_era'
                [ghm_main,~,~] = mGlobe_getModelPath;
                output_path = fullfile(ghm_main,'ERA');
                start_calc = datenum(str2double(get(findobj('Tag','edit_down_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_down_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_start_hour'),'String')),0,0);
                end_calc = datenum(str2double(get(findobj('Tag','edit_down_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_down_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_end_hour'),'String')),0,0);
                step_calc = get(findobj('Tag','popup_down_time_step'),'Value'); % Get time step
                input_file = get(findobj('Tag','push_down_era_input'),'UserData');
                if start_calc > end_calc 									% check if starting time > end time
                    set(findobj('Tag','text_status'),'String','Start time must be <= End time');
                elseif ~isempty(input_file)
                    set(findobj('Tag','text_status'),'String','Models: conversion starting...');drawnow
                    mGlobe_convert_ERA(start_calc,end_calc,step_calc,input_file,output_path) % covert ERA model
                    pause(5);
                    set(findobj('Tag','text_status'),'String','Set your conversion options');
                else
                    set(findobj('Tag','text_status'),'String','Please choose your ERA netCDF input file');
                end
            case 'calc_down_ncep'
                model_ver = menu('Choose NCEP version','Reanalysis 1 (beta)','Reanalysis 2');
                if model_ver ~= 0
                    [ghm_main,~,~] = mGlobe_getModelPath;
                    output_path = fullfile(ghm_main,'NCEP');
                    start_calc = datenum(str2double(get(findobj('Tag','edit_down_time_start_year'),'String')),... % Get date of start
                        str2double(get(findobj('Tag','edit_down_time_start_month'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_start_day'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_start_hour'),'String')),0,0);
                    end_calc = datenum(str2double(get(findobj('Tag','edit_down_time_end_year'),'String')),... % Get date of end
                        str2double(get(findobj('Tag','edit_down_time_end_month'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_end_day'),'String')),...
                        str2double(get(findobj('Tag','edit_down_time_end_hour'),'String')),0,0);
                    step_calc = get(findobj('Tag','popup_down_time_step'),'Value'); % Get time step
                    input_file = get(findobj('Tag','push_down_era_input'),'UserData');
                    if start_calc > end_calc 									% check if starting time > end time
                        set(findobj('Tag','text_status'),'String','Start time must be <= End time');
                    elseif ~isempty(input_file)
                        input_path = fileparts(input_file);
                        set(findobj('Tag','text_status'),'String','Models: conversion starting...');drawnow
                        mGlobe_convert_NCEP(start_calc,end_calc,step_calc,input_path,output_path,model_ver) % covert NCEP model
                        pause(5);
                        set(findobj('Tag','text_status'),'String','Set your conversion options');
                    else
                        set(findobj('Tag','text_status'),'String','Please choose your NCEP netCDF input file.');
                    end
                else
                    set(findobj('Tag','text_status'),'String','Please select NCEP Reanalysis version');
                end
            % DEM
            case 'load_down_dem_input'                                      % load input dem file                          
                [name,path] = uigetfile('*.*','Choose your input DEM file - Latitude,Longitude,Height');
                if name == 0
                    set(findobj('Tag','push_down_dem_input'),'UserData',[]);
                    set(findobj('Tag','text_down_dem_input'),'String','no DEM loaded');
                else
                    set(findobj('Tag','push_down_dem_input'),'UserData',[path,name]);
                    set(findobj('Tag','text_down_dem_input'),'String',name);
                end   
            case 'load_down_dem_output'                                     % choose DEM output file     
                [name,path] = uiputfile('*.*','Choose your output DEM file - matlab file *.mat');
                if name == 0
                    set(findobj('Tag','push_down_dem_output'),'UserData',[]);
                    set(findobj('Tag','text_down_dem_output'),'String','no output file');
                else
                    set(findobj('Tag','push_down_dem_output'),'UserData',[path,name]);
                    set(findobj('Tag','text_down_dem_output'),'String',name);
                end 
            case 'calc_down_dem'                                            % transform loaded DEM to supported file format (structure area)
                DEM_input = get(findobj('Tag','push_down_dem_input'),'UserData');
                DEM_output = get(findobj('Tag','push_down_dem_output'),'UserData');
                DEM_type = get(findobj('Tag','popup_down_dem_type'),'Value');
                if ~isempty(DEM_input) && ~isempty(DEM_output)
                    set(findobj('Tag','text_status'),'String','Models: Starting conversion...');drawnow
                    mGlobe_convert_DEM(DEM_input,DEM_output,DEM_type); 		% perform transformation
                    pause(5);
                    set(findobj('Tag','text_status'),'String','Set your conversion options');
                else
                    set(findobj('Tag','text_status'),'String','Choose your input and output DEM file');
                end
            % OTHER                                         
            case 'load_down_other_input'                                    % load an input file for the conversion (fixed prefix is required)                           
                [name,path] = uigetfile('*.*','Choose your input file: PREFIX1234_YYYYMMDD_HH.txt|dat|...');
                if name == 0
                    set(findobj('Tag','push_down_other_input'),'UserData',[]);
                    set(findobj('Tag','text_status'),'String','No input file');
                else
                    set(findobj('Tag','push_down_other_input'),'UserData',path);
                    set(findobj('Tag','popup_down_other_type'),'UserData',name);
                    if length(name)<24
                        set(findobj('Tag','text_status'),'String','Wrong input file name: too short');
                    else
                        message = sprintf('File prefix: %s, date: %s, time: %s, extension: %s',name(1:10),name(12:19),name(21:22),name(24:end));
                        set(findobj('Tag','text_status'),'String',message);clear message
                    end
                end 
            case 'calc_down_other'                                          % transform loaded input OTHER model to supported file format (structure area)
                output_path = get(findobj('Tag','push_down_other_path'),'UserData');
                start_calc = datenum(str2double(get(findobj('Tag','edit_down_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_down_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_start_hour'),'String')),0,0);
                end_calc = datenum(str2double(get(findobj('Tag','edit_down_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_down_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_end_hour'),'String')),0,0);
                step_calc = get(findobj('Tag','popup_down_time_step'),'Value'); % Get time step
                input_path = get(findobj('Tag','push_down_other_input'),'UserData');
                input_name = get(findobj('Tag','popup_down_other_type'),'UserData');
                file_type = get(findobj('Tag','popup_down_other_type'),'Value');
                header_lines = str2double(get(findobj('Tag','edit_down_other_header'),'String'));
                
                if start_calc > end_calc
                    set(findobj('Tag','text_status'),'String','Start time must be <= End time');
                elseif ~isempty(input_name) && length(input_name)>23
                    set(findobj('Tag','text_status'),'String','Models: starting conversion...');drawnow
                    mGlobe_convert_OTHER(start_calc,end_calc,step_calc,input_name,input_path,output_path,file_type,header_lines);
                    set(findobj('Tag','text_status'),'String','Models: Conversion completed');
                    pause(5);
                    set(findobj('Tag','text_status'),'String','Set your conversion options');
                else
                    set(findobj('Tag','text_status'),'String','Please choose your (correct) input file');
                end
            % GRACE
            case 'load_down_grace_input'                                    % load GRACE mass grid (one file for all time epochs)                        
                [name,path] = uigetfile('*.nc','Choose your input TELLUS GRACE netCDF file');
                if name == 0
                    set(findobj('Tag','push_down_grace_input'),'UserData',[]);
                    set(findobj('Tag','text_down_grace_input'),'String',[]);
                else
                    set(findobj('Tag','push_down_grace_input'),'UserData',[path,name]);
                    set(findobj('Tag','text_down_grace_input'),'String',name);
                end 
            case 'select_grace_model'                                       % switch between the Land and Ocean mass grid
                ocean_land = get(findobj('Tag','popup_down_grace_land'),'Value');
                if ocean_land == 1
                    set(findobj('Tag','text_down_grace_land'),'String','Path+prefix:/GRACE/LAND/');
                    set(findobj('Tag','check_down_grace_land'),'Value',1);
                else
                    set(findobj('Tag','text_down_grace_land'),'String','Path+prefix:/GRACE/OCEAN/');
                    set(findobj('Tag','check_down_grace_land'),'Value',0);
                end
            case 'calc_down_grace'                                          % convert GRACE mass grids to supported file format
                file_path = get(findobj('Tag','push_down_grace_input'),'UserData');
                start_calc = datenum(str2double(get(findobj('Tag','edit_down_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_down_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_start_hour'),'String')),0,0);
                end_calc = datenum(str2double(get(findobj('Tag','edit_down_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_down_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_down_time_end_hour'),'String')),0,0);
                output_name = get(findobj('Tag','edit_down_grace_name'),'String'); % Get time step
                input_file = get(findobj('Tag','push_down_grace_input'),'UserData');
                ocean_land = get(findobj('Tag','popup_down_grace_land'),'Value');
                scale_mat = get(findobj('Tag','check_down_grace_land'),'Value');
                if scale_mat == 1 && ~isempty(input_file)                   % load scaling matrix if required (only for land mass grids)
                    [name,path] = uigetfile('*.*','Choose your input SCALE MATRIX (TELLUS GRACE) netCDF file');
                    set(findobj('Tag','text_status'),'String','Choose your input SCALE MATRIX (TELLUS GRACE) netCDF file');
                    if name == 0
                        scale_file = [];
                    else
                        scale_file = [path,name];
                    end 
                else
                    scale_file = [];
                end
                if (start_calc > end_calc) || length(output_name) ~= 22     % fixed file format is required, i.e. file prefix
                    set(findobj('Tag','text_status'),'String','Start time must be <= End time and Output PREFIX must have 22 letters exactly');
                elseif ~isempty(input_file)
                    set(findobj('Tag','text_status'),'String','Models: conversion starting...');drawnow
                    [~,~,grace_main] = mGlobe_getModelPath;
                    mGlobe_convert_GRACE_tellus(start_calc,end_calc,file_path,output_name,grace_main,ocean_land,scale_file); % convert GRACE model
                    set(findobj('Tag','text_status'),'String','Conversion completed');
                    pause(5);
                    set(findobj('Tag','text_status'),'String','Set your conversion options');
                else
                    set(findobj('Tag','text_status'),'String','Please choose your (valid) TELLUS GRACE netCDF input file');
                end
                
                %% OCEAN
            case 'load_ocean_convert_input'                                 % load files required for ECCO1/2 model conversion
                model_version = get(findobj('Tag','popup_ocean_load_convert_model'),'Value');
                step_calc = get(findobj('Tag','popup_ocean_time_step'),'Value'); % Get time step
                switch model_version
                    case 1
                        if step_calc == 6
                            [name,path] = uigetfile('*.*','Choose your input file: ECCO_kf080_YYYYDOY_YYYYDOY_AveRmvd_OBP.txt');
                        else
                            [name,path] = uigetfile('*.*','Choose your input file: e.g. OBPano_08_08.02160_04320_012.cdf');
                        end
                    case 4
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. PHIBOT.1440x720.20000107.nc');
                    case 5
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. AOD1B_2009-01-01_X_05.asc');
                    case 6
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. AOD1B_2009-01-01_X_05.asc');
                    case 7
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. AOD1B_2009-01-01_X_05.asc');
                    case 8
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. AOD1B_2009-01-01_X_06.asc');
                    case 9
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. AOD1B_2009-01-01_X_06.asc');
                    case 10
                        [name,path] = uigetfile('*.*','Choose your input file: e.g. AOD1B_2009-01-01_X_06.asc');
                end
                if name == 0
                    set(findobj('Tag','push_ocean_load_convert_input'),'UserData',[]);
                    set(findobj('Tag','text_status'),'String','No input file');
                    set(findobj('Tag','text_ocean_convert_input'),'String','Please choose ECCO file');
                else
                    set(findobj('Tag','push_ocean_load_convert_input'),'UserData',name);
                    set(findobj('Tag','text_ocean_convert_input'),'String',name);
                    set(findobj('Tag','text_ocean_convert_input'),'UserData',path);
                end 
            case 'popup_ocean_conservation'
                model_conserv = get(findobj('Tag','popup_ocean_load_conserv'),'Value');
                if model_conserv == 3
                    set(findobj('Tag','push_ocean_load_time_series'),'Visible','on');
                    set(findobj('Tag','text_ocean_load_time_series'),'Visible','on');
                    set(findobj('Tag','text_ocean_load_time_series_time'),'Visible','on');
                    set(findobj('Tag','text_ocean_load_time_series_press'),'Visible','on');
                    set(findobj('Tag','edit_ocean_load_time_series_time'),'Visible','on');
                    set(findobj('Tag','edit_ocean_load_time_series_press'),'Visible','on');
                else
                    set(findobj('Tag','push_ocean_load_time_series'),'UserData',[]);
                    set(findobj('Tag','push_ocean_load_time_series'),'Visible','off');
                    set(findobj('Tag','text_ocean_load_time_series'),'UserData',[]);
                    set(findobj('Tag','text_ocean_load_time_series'),'Visible','off');
                    set(findobj('Tag','text_ocean_load_time_series_time'),'Visible','off');
                    set(findobj('Tag','text_ocean_load_time_series_press'),'Visible','off');
                    set(findobj('Tag','edit_ocean_load_time_series_time'),'Visible','off');
                    set(findobj('Tag','edit_ocean_load_time_series_press'),'Visible','off');
                    set(findobj('Tag','text_ocean_load_time_series'),'String','Select input file (*.txt)');
                end
            case 'load_ocean_time_series'
                [name,path] = uigetfile('*.*','Choose your pressure time series (units: Pa, format: *.txt)');
                if name == 0
                    set(findobj('Tag','push_ocean_load_time_series'),'UserData',[]);
                    set(findobj('Tag','text_status'),'String','No input file');
                    set(findobj('Tag','text_ocean_load_time_series'),'String','Please select your pressure time series');
                else
                    set(findobj('Tag','push_ocean_load_time_series'),'UserData',[path,name]);
                    set(findobj('Tag','text_ocean_load_time_series'),'String',name);
                end
            case 'load_ocean_out'                                           % Set NTOL output file
                [name,path] = uiputfile('*.*','Output file: excel sheet or txt');
                if name == 0                                                % If cancelled-> default file
                    set(findobj('Tag','push_ocean_out'),'UserData','output.txt');
                    set(findobj('Tag','text_ocean_out'),'String','output.txt');
                else
                    set(findobj('Tag','push_ocean_out'),'UserData',[path,name]);
                    set(findobj('Tag','text_ocean_out'),'String',name);
                end
                
            case 'popup_ocean_model'                                        % switch between OBP models
                val = get(findobj('Tag','popup_ocean_load_convert_model'),'Value');
                [~,obpm_main,grace_main] = mGlobe_getModelPath;
                switch val
                    case 1
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'ECCO1'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'ECCO1'));
                    case 2
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OTHER'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OTHER'));
                    case 3
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(grace_main,'OCEAN'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(grace_main,'OCEAN'));
                    case 4
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'ECCO2'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'ECCO2'));
                    case 5
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OMCT'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OMCT'));
                    case 6
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OMCT'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OMCT'));
                    case 7
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OMCT'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OMCT'));
                    case 8
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OMCT6'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OMCT6'));
                    case 9
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OMCT6'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OMCT6'));
                    case 10
                        set(findobj('Tag','text_ocean_convert_output'),'UserData',fullfile(obpm_main,'OMCT6'));
                        set(findobj('Tag','text_ocean_convert_output'),'String',fullfile(obpm_main,'OMCT6'));
                end
            case 'load_ocean_convert_calc'                                  % Convert ECCO1/2 ocean bottom pressure models to supported file format
                file = get(findobj('Tag','push_ocean_load_convert_input'),'UserData');
                start_calc = datenum(str2double(get(findobj('Tag','edit_ocean_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_ocean_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_start_hour'),'String')),0,0);
                
                end_calc = datenum(str2double(get(findobj('Tag','edit_ocean_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_ocean_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_end_hour'),'String')),0,0);
                
                step_calc = get(findobj('Tag','popup_ocean_time_step'),'Value'); % Get time step
                model_version = get(findobj('Tag','popup_ocean_load_convert_model'),'Value');
                obp_model = get(findobj('Tag','text_ocean_convert_output'),'UserData');
                input_path = get(findobj('Tag','text_ocean_convert_input'),'UserData');
                if start_calc > end_calc
                    set(findobj('Tag','text_status'),'String','Start time must be <= End time');
                elseif (~isempty(file) && model_version == 1) || (~isempty(file) && model_version == 4)
                    set(findobj('Tag','text_status'),'String','Ocean: starting ECCO conversion...');drawnow
                    mGlobe_convert_ECCO(start_calc,end_calc,step_calc,file,obp_model,input_path,model_version); % convert ECCO model data
                    set(findobj('Tag','text_status'),'String','Ocean: Conversion completed');
                    pause(5);
                    set(findobj('Tag','text_status'),'String','Set your non-tidal ocean loading effect');
                elseif (~isempty(file) && model_version >= 5)
                    set(findobj('Tag','text_status'),'String','Ocean: starting OMCT conversion (this step requires some time)...');drawnow
                    mGlobe_convert_OMCT(start_calc,end_calc,step_calc,file,obp_model,input_path,model_version); % convert OMCT model data
                    set(findobj('Tag','text_status'),'String','Ocean: Conversion completed');
                    pause(5);
                    set(findobj('Tag','text_status'),'String','Set your non-tidal ocean loading effect');
                else
                    set(findobj('Tag','text_status'),'String','Please choose your (ECCO/OMCT) input file');
                end
                
            case 'ocean_calc'                                              % Calculate the non-tidal ocean loading effect (NTOL)
                Input = [str2double(get(findobj('Tag','edit_ocean_pos_lat'),'String')),... % Get input coordinates
                    str2double(get(findobj('Tag','edit_ocean_pos_lon'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_pos_hei'),'String'))];
                output_file = get(findobj('Tag','push_ocean_out'),'UserData'); % Get output file
                output_file_type = [get(findobj('Tag','check_ocean_xls'),'Value'),...
                                    get(findobj('Tag','check_ocean_txt'),'Value'),...
                                    get(findobj('Tag','check_ocean_tsf'),'Value')];
                
                step_calc = get(findobj('Tag','popup_ocean_time_step'),'Value'); % Get time step
                if output_file_type(3) == 1  && step_calc == 6
                    set(findobj('Tag','check_ocean_tsf'),'Value',0)
                    output_file_type(3) =  0;
                end
                if sum(output_file_type) == 0
                    set(findobj('Tag','check_ocean_txt'),'Value',1);
                    output_file_type(2) = 1;
                    drawnow
                end
                
                start_calc = datenum(str2double(get(findobj('Tag','edit_ocean_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_ocean_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_start_hour'),'String')),0,0);
                
                end_calc = datenum(str2double(get(findobj('Tag','edit_ocean_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_ocean_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_ocean_time_end_hour'),'String')),0,0);
                subtract_average = get(findobj('Tag','check_ocean_average'),'Value');
                
                ghc_path = get(findobj('Tag','text_ocean_convert_output'),'UserData');
                model_version = get(findobj('Tag','popup_ocean_load_convert_model'),'Value');                
                ghc_treshold = str2double(get(findobj('Tag','edit_ocean_treshold'),'String')); % Get calculation (distance) threshold
                mean_field = get(findobj('Tag','popup_ocean_load_conserv'),'Value');
                pressure_time_series{1} = get(findobj('Tag','push_ocean_load_time_series'),'UserData');
                pressure_time_series{2} = {str2double(get(findobj('Tag','edit_ocean_load_time_series_time'),'String'))};
                pressure_time_series{3} = {str2double(get(findobj('Tag','edit_ocean_load_time_series_press'),'String'))};
                name = 1;
                shp_file = get(findobj('Tag','push_ocean_load_shp'),'UserData');
                if model_version == 2                                       % prompt user to select OTHER model (fixed prefix is required)
                    set(findobj('Tag','text_status'),'String','You have chosen OTHER model => pick file with *.mat grid, first TEN letters will be used as a PREFIX (e.g. MODEL_1234)');
                    [name,~] = uigetfile(ghc_path,'Pick file with *.mat grid, first TEN letters will be used as a PREFIX for data loading (e.g. MODEL_1234)');
                    set(findobj('Tag','text_ocean_convert_output'),'String',ghc_path,'UserData',ghc_path);drawnow
                    if name ~=0
                        ghc_path = fullfile(ghc_path,name(1:10));
                        set(findobj('Tag','text_status'),'String',['File path + prefix = ',ghc_path]);
                    end
                end
                if model_version == 3                                       % prompt user to select GRACE model (fixed prefix is required)
                    [~,~,grace_main] = mGlobe_getModelPath;
                    set(findobj('Tag','text_status'),'String','You have chosen GRACE model => pick file with *.mat grid, first 22 letters will be used as a PREFIX (e.g. GRC_JPL_RL05_FILT0_OCE)');
                    [name,~] = uigetfile(grace_main,'Pick file with *.mat grid, first 22 letters will be used as a PREFIX for data loading (e.g. GRC_JPL_RL05_FILT0_OCE)');
                    set(findobj('Tag','text_ocean_convert_output'),'String',ghc_path,'UserData',ghc_path);drawnow
                    if name ~=0
                        ghc_path = fullfile(ghc_path,name(1:22));
                        set(findobj('Tag','text_status'),'String',['File path + prefix = ',ghc_path]);
                    end
                end
                if mean_field == 3 && isempty(pressure_time_series{1})
                    name = 0;
                elseif (mean_field == 3 && ~isempty(pressure_time_series{1})) && (cell2mat(pressure_time_series{2})+cell2mat(pressure_time_series{3})<=1)
                    name = 0;
                end
                    
                if  (ghc_treshold > 1 || ghc_treshold < 0.05) || (start_calc > end_calc) || (sum(double(name)) == 0) % check if everything is set correctly
                    set(findobj('Tag','text_status'),'String','Threshold must be within <0.05,1.00> degree, start time <= end time and model set correctly');
                elseif exist('mGlobe_DATA_dgE_Hydro.txt','file')==2 && exist('mGlobe_DATA_OceanGrid.mat','file')==2
                    set(findobj('Tag','text_status'),'String','Ocean: starting the computation...');drawnow
                    mGlobe_calc_Ocean(Input,output_file,output_file_type,start_calc,end_calc,...
                                    step_calc,ghc_treshold,ghc_path,model_version,...
                                    subtract_average,mean_field,...
                                    pressure_time_series,shp_file); % start the NTOL computation
                    pause(5)                                                % wait 5 sec, than write original message
                    set(findobj('Tag','text_status'),'String','Set your non-tidal ocean loading effect');
                else
                    set(findobj('Tag','text_status'),'String',...           % warn user that these two file must be located in the current folder
                    'Please ensure that both mGlobe_DATA_dgE_Hydro.txt and mGlobe_DATA_OceanGrid.mat file are on your MATLAB current folder'); drawnow 
                end 
                %% ATMO
                case 'load_atmo_geopotential'                               % load netCDF ERA geopotential data
                    [name,path] = uigetfile('*.nc','Geopotential height, ERA: netcdf (fixed suffix: xyz_HH_YYYY.nc), MERRA: netcdf (fixed suffix: YYYYMMDD.SUB.nc)');
                    if name == 0                                            % If cancelled-> no data
                        set(findobj('Tag','push_atmo_geopotential'),'UserData',[]);
                        set(findobj('Tag','text_atmo_load_geopotential'),'String','No data loaded');
                    else
                        set(findobj('Tag','push_atmo_geopotential'),'UserData',[path,name]);
                        set(findobj('Tag','text_atmo_load_geopotential'),'String',name);
                    end
                case 'load_atmo_humidity'                                   % load netCDF ERA humidity data
                    [name,path] = uigetfile('*.nc','Spec. humidity, ERA: netcdf (fixed suffix: xyz_HH_YYYY.nc), MERRA: netcdf (fixed suffix: YYYYMMDD.SUB.nc)');
                    if name == 0                                            % If cancelled-> no data
                        set(findobj('Tag','push_atmo_humidity'),'UserData',[]);
                        set(findobj('Tag','text_atmo_load_humidity'),'String','No data loaded');
                    else
                        set(findobj('Tag','push_atmo_humidity'),'UserData',[path,name]);
                        set(findobj('Tag','text_atmo_load_humidity'),'String',name);
                    end
                case 'load_atmo_temperature'                                % load netCDF ERA temperature data
                    [name,path] = uigetfile('*.nc','Temperature, ERA: netcdf (fixed suffix: xyz_HH_YYYY.nc), MERRA: netcdf (fixed suffix: YYYYMMDD.SUB.nc)');
                    if name == 0                                            % If cancelled-> no data
                        set(findobj('Tag','push_atmo_temperature'),'UserData',[]);
                        set(findobj('Tag','text_atmo_load_temperature'),'String','No data loaded');
                    else
                        set(findobj('Tag','push_atmo_temperature'),'UserData',[path,name]);
                        set(findobj('Tag','text_atmo_load_temperature'),'String',name);
                    end
                case 'load_atmo_surface'                                    % load netCDF ERA surface data
                    [name,path] = uigetfile('*.nc','Surface data, ERA: netcdf (2m temp.,2m dewp.,press.: suffix xyz_YYYY.nc), MERRA: netcdf (temp., humid.:suffix YYYYMMDD.SUB.nc)');
                    if name == 0                                            % If cancelled-> no data
                        set(findobj('Tag','push_atmo_surface'),'UserData',[]);
                        set(findobj('Tag','text_atmo_load_surface'),'String','No data loaded');
                    else
                        set(findobj('Tag','push_atmo_surface'),'UserData',[path,name]);
                        set(findobj('Tag','text_atmo_load_surface'),'String',name);
                    end
                case 'load_atmo_orography'                                  % load netCDF ERA orography data
                    [name,path] = uigetfile('*.*','Orography, ERA: netcdf (*.nc), MERRA: hdf or nc (*.hdf| *.nc)');
                    if name == 0                                            % If cancelled-> no data
                        set(findobj('Tag','push_atmo_orography'),'UserData',[]);
                        set(findobj('Tag','text_atmo_load_orography'),'String','No data loaded');
                    else
                        set(findobj('Tag','push_atmo_orography'),'UserData',[path,name]);
                        set(findobj('Tag','text_atmo_load_orography'),'String',name);
                    end
                case 'load_atmo_merra'                                      % load MERRA input data
                    [name,path] = uigetfile('*.nc','MERRA surface pressure file (suffix: YYYYMMDD.SUB.nc)');
                    if name == 0                                            % If cancelled-> no data 
                        set(findobj('Tag','push_atmo_merra'),'UserData',[]);
                        set(findobj('Tag','text_atmo_load_merra'),'String','No data loaded');
                    else
                        set(findobj('Tag','push_atmo_merra'),'UserData',[path,name]); % path only
                        set(findobj('Tag','text_atmo_load_merra'),'String',name);
                    end
            case 'load_atmo_out'                                            % Set output file
                [name,path] = uiputfile('*.*','Output file: excel sheet or txt');
                if name == 0                                                % If cancelled-> default output file
                    set(findobj('Tag','push_atmo_out'),'UserData','output.txt');
                    set(findobj('Tag','text_atmo_out'),'String','output.txt');
                else
                    set(findobj('Tag','push_atmo_out'),'UserData',[path,name]);
                    set(findobj('Tag','text_atmo_out'),'String',name);
                end
            case 'atmo_calc_era'                                                % Calculate atmospheric effect
                Input = [str2double(get(findobj('Tag','edit_atmo_pos_lat'),'String')),... % Get input coordinates
                    str2double(get(findobj('Tag','edit_atmo_pos_lon'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_pos_hei'),'String'))];
                output_file = get(findobj('Tag','push_atmo_out'),'UserData'); % Get output file
                output_file_type = [get(findobj('Tag','check_atmo_xls'),'Value'),...
                                    get(findobj('Tag','check_atmo_txt'),'Value'),...
                                    get(findobj('Tag','check_atmo_tsf'),'Value')];
                                
                step_calc = get(findobj('Tag','popup_atmo_time_step'),'Value')+1; % Get time step +1 (+1 because the default popupmenu is for hydrological model = 3Hour step)
                
                if sum(output_file_type) == 0                               % set txt output file if nothing is selected by user
                    set(findobj('Tag','check_atmo_txt'),'Value',1);
                    output_file_type(2) = 1;
                    drawnow
                end
                subtract_average = get(findobj('Tag','check_atmo_average'),'Value');
                file_ref = get(findobj('Tag','push_atmo_orography'),'UserData');
                file_temp = get(findobj('Tag','push_atmo_temperature'),'UserData'); 
                file_humid = get(findobj('Tag','push_atmo_humidity'),'UserData'); 
                file_height = get(findobj('Tag','push_atmo_geopotential'),'UserData'); 
                file_sp = get(findobj('Tag','push_atmo_surface'),'UserData'); 
                check_load = isempty(file_ref) + isempty(file_temp) + isempty(file_humid) + isempty(file_height) + isempty(file_sp);
                
                try                                                         % check if every file is OK
                    ncid_ref = netcdf.open(file_ref,'NC_NOWRITE');          % try to open all files (check)
                    netcdf.close(ncid_ref);
                    ncid_ref = netcdf.open(file_temp,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref= netcdf.open(file_humid,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref = netcdf.open(file_height,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref = netcdf.open(file_sp,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                catch
                    check_load = 1;
                end
                
                start_calc = datenum(str2double(get(findobj('Tag','edit_atmo_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_atmo_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_start_hour'),'String')),0,0);
                
                end_calc = datenum(str2double(get(findobj('Tag','edit_atmo_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_atmo_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_end_hour'),'String')),0,0);
                
                if start_calc > end_calc
                    check_load = 1;
                end
                if check_load ~= 0
                    set(findobj('Tag','text_status'),'String','Check time settings and load all valid input data (temperature, humidity, geopotential, surface data and orography)');
                elseif exist('mGlobe_DATA_dgE_Atmo.txt','file')==2 && exist('mGlobe_DATA_OceanGrid.mat','file')==2
                    set(findobj('Tag','text_status'),'String','Atmo: starting the computation...');drawnow
                    mGlobe_calc_Atmo_ERA(Input,output_file,output_file_type,file_ref,file_temp,file_humid,file_height,file_sp,start_calc,end_calc,step_calc,subtract_average);
                    pause(5)                                                % wait 5 sec, than write output message
                    set(findobj('Tag','text_status'),'String','Set the global atmospheric effect');
                else
                    set(findobj('Tag','text_status'),'String',...           % warn user that these two files must be located in the current folder
                    'Please ensure that both mGlobe_DATA_dgE_Atmo.txt and mGlobe_DATA_OceanGrid.mat file are on your MATLAB current folder'); drawnow 
                end
            case 'atmo_calc_merra'                                          % Calculate MERRA atmospheric effect
                Input = [str2double(get(findobj('Tag','edit_atmo_pos_lat'),'String')),... % Get input coordinates
                    str2double(get(findobj('Tag','edit_atmo_pos_lon'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_pos_hei'),'String'))];
                output_file = get(findobj('Tag','push_atmo_out'),'UserData'); % Get output file
                output_file_type = [get(findobj('Tag','check_atmo_xls'),'Value'),...
                                    get(findobj('Tag','check_atmo_txt'),'Value'),...
                                    get(findobj('Tag','check_atmo_tsf'),'Value')];
                                
                step_calc = get(findobj('Tag','popup_atmo_time_step'),'Value')+1; % Get time step +1 (+1 because the default popupmenu is for hydrological model = 3Hour step)
                
                if sum(output_file_type) == 0                               % set txt output file if nothing is selected by user
                    set(findobj('Tag','check_atmo_txt'),'Value',1);
                    output_file_type(2) = 1;
                    drawnow
                end
                subtract_average = get(findobj('Tag','check_atmo_average'),'Value');
                file_ref = get(findobj('Tag','push_atmo_orography'),'UserData');
                file_temp = get(findobj('Tag','push_atmo_temperature'),'UserData'); 
                file_humid = get(findobj('Tag','push_atmo_humidity'),'UserData'); 
                file_height = get(findobj('Tag','push_atmo_geopotential'),'UserData'); 
                file_sthd = get(findobj('Tag','push_atmo_surface'),'UserData'); 
                file_sp = get(findobj('Tag','push_atmo_merra'),'UserData'); % merra surface pressure file (not for ERA)
                check_load = isempty(file_ref) + isempty(file_temp) + isempty(file_humid) + isempty(file_height) + isempty(file_sp) + isempty(file_sthd);
                
                try                                                         % check if every file is OK
                    switch file_ref(end-2:end);
                        case 'hdf'
                            ncid_ref = hdfinfo(file_ref);                           % try to open all files (check), orography in hdf file, all others in netcdf file format
                        case '.nc'
                            ncid_ref = netcdf.open(file_ref,'NC_NOWRITE');
                            netcdf.close(ncid_ref);
                    end
                    ncid_ref = netcdf.open(file_temp,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref= netcdf.open(file_humid,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref = netcdf.open(file_height,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref = netcdf.open(file_sp,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                    ncid_ref = netcdf.open(file_sthd,'NC_NOWRITE');
                    netcdf.close(ncid_ref);
                catch
                    check_load = 1;
                end
                
                start_calc = datenum(str2double(get(findobj('Tag','edit_atmo_time_start_year'),'String')),... % Get date of start
                    str2double(get(findobj('Tag','edit_atmo_time_start_month'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_start_day'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_start_hour'),'String')),0,0);
                
                end_calc = datenum(str2double(get(findobj('Tag','edit_atmo_time_end_year'),'String')),... % Get date of end
                    str2double(get(findobj('Tag','edit_atmo_time_end_month'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_end_day'),'String')),...
                    str2double(get(findobj('Tag','edit_atmo_time_end_hour'),'String')),0,0);
                
                if start_calc > end_calc
                    check_load = 1;
                end
                if check_load ~= 0
                    set(findobj('Tag','text_status'),'String','Check time settings and load all valid input data (temperature, humidity, geopotential, surface data and orography)');
                elseif exist('mGlobe_DATA_dgE_Atmo.txt','file')==2 && exist('mGlobe_DATA_OceanGrid.mat','file')==2
                    set(findobj('Tag','text_status'),'String','Atmo: starting the computation...');drawnow
                    mGlobe_calc_Atmo_MERRA(Input,output_file,output_file_type,file_ref,file_temp,file_humid,file_height,file_sthd,file_sp,start_calc,end_calc,step_calc,subtract_average);
                    pause(5)                                                % wait 5 sec, than write output message
                    set(findobj('Tag','text_status'),'String','Set the global atmospheric effect');
                else
                    set(findobj('Tag','text_status'),'String',...           % warn user that these two files must be located in the current folder
                    'Please ensure that both mGlobe_DATA_dgE_Atmo.txt and mGlobe_DATA_OceanGrid.mat file are on your MATLAB current folder'); drawnow 
                end
        end

    end
function [ghm_main,obpm_main,grace_main] = mGlobe_getModelPath
    try
        er = 0;
        fid = fopen('mGlobe_PATH_Settings.txt','r');
        if fid > 0
            setting = textscan(fid,'%s %s','Delimiter','>','commentstyle','%');
            if length(setting{1}) ~= 3
                er = 1;
            else
                ghm_main = [];grace_main = [];
                for ii = 1:length(setting{1})
                    switch char(setting{1}(ii))
                        case 'GHM'
                            ghm_main = char(setting{2}(ii));
                        case 'GRACE'
                            grace_main = char(setting{2}(ii));
                        case 'OBPM'
                            obpm_main = char(setting{2}(ii));
                    end
                end
            end
        else
            er = 1;
        end
        fclose(fid);
    catch
        er = 1;
    end
    if er == 1 || isempty(ghm_main) || isempty(grace_main) || isempty(obpm_main)
        disp('Could not set model paths!');
        % Use default values if error occurred
        ghm_main = 'GHM';
        grace_main = 'GRACE';
        obpm_main = 'OBPM';
    end
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  KONIEC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                