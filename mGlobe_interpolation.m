function [DataN,DataID]= mGlobe_interpolation(LonI,LatI,DataI,LonN,LatN,continent)
%MGLOBE_INTERPOLATION Interpolation of spherical/ellipsoidal data
% Function is used for the interpolation of data in ellipsoidal
% coordinates, i.e. with a jump at the 180 deg meridian, e.g. for 
% longitude [..179 179.5 -179.5 -179..]
% The second part serves for the identification of continental areas ->
% oceans = -1, continent = 1. Function works only with global coverage
% data and file: mGlobe_DATA_OceanGrid.mat in the same folder!
% 
% Input:    
%   LatI      ... given (input) latitude created by meshgrid (deg)
%   LonI      ... given (input) longitude created by meshgrid (deg)
%   DataI     ... given (input) data
%   LatN      ... latitude created by meshgrid used for interpolation
%   LonI      ... longitude created by meshgrid used for interpolation
%   continent ... 1 for continent/ocean determination = 1/-1
%                 0 without continent/ocean area identification
%                 string for use of shapefile for continent/ocean ID
% Output:   
%   DataN     ... new interpolated data for LonN, LatN grid
%   DataI     ... ID for ocean (==-1) and land (==1)
% 
%                                         M.Mikolaj, mikolaj@gfz-potsdam.de
%                                                                18.06.2014
%                                                                      v1.0

%% Prepare grid = sort after longitude
% Before transforming coordinates, check if SHP file on input for
% Ocean/Land identification
if ischar(continent)
    boarders = shaperead(continent,'UseGeoCoords',true);
    DataID = LonN*0 - 0.5;
    for i = 1:length(boarders)
        DataID = DataID + double(inpolygon(LonN,LatN,[boarders(i).Lon],[boarders(i).Lat]));
    end
    DataID(DataID>0) = 1;
end
LonI(LonI<0) = LonI(LonI<0)+360;                                            % Transform given longitude to <0,360) system
LonN(LonN<0) = LonN(LonN<0)+360;                                            % Transform new longitude to <0,360) system
r = find(abs(diff(LonI(1,:)))==max(abs(diff(LonI(1,:)))));                  % find the point of discontinuity, e.g. ...,359,1,2...               
LonI = horzcat(LonI(:,r+1:end),LonI(:,1:r));                                % Connect matrices to remove discontinuity
LatI = horzcat(LatI(:,r+1:end),LatI(:,1:r));
DataI = horzcat(DataI(:,r+1:end),DataI(:,1:r));
LonI = LonI + 10;                                                           % Shift Longitude by 10 deg
LonN = LonN + 10;                                                           % Shift Longitude by 10 deg
diffLonI = abs(LonI(1,1)-LonI(1,2));                                        % determinate longitudal step
LonI = horzcat(LonI(:,1)-diffLonI,LonI,LonI(:,end)+diffLonI);               % ad first end last column => eliminate blind area, e.g. between 359 an 1
LatI = horzcat(LatI(:,1),LatI,LatI(:,end));                                 % ad first end last column
DataI = horzcat(DataI(:,end),DataI,DataI(:,1));                             % ad first end last column

%% Interpolate new data
DataN = interp2(LonI,LatI,DataI,LonN,LatN);                                 % interpolate data (linear) => no extrapolation/NaN between 359 an 1 deg
clear LonI LatI DataI diffLonI r
%% Identify areas on ocean/continent
if ~ischar(continent)
    if continent == 1
        load('mGlobe_DATA_OceanGrid.mat');                                        % load auxiliary grid of oceans
        [oceans.lon,oceans.lat] = meshgrid(oceans.lon,oceans.lat);
        oceans.lon(oceans.lon<0) = oceans.lon(oceans.lon<0) + 360;              % Transform ocean grid longitude to <0,360) system
        oceans.lon = oceans.lon + 10;                                           % Shift Longitude by 10 deg
        r = find(abs(diff(oceans.lon(1,:)))==max(abs(diff(oceans.lon(1,:)))));  % find the point of discontinuity, e.g. ...,359,1,2...               
        oceans.lon = horzcat(oceans.lon(:,r+1:end),oceans.lon(:,1:r));          % Connect matrices to remove discontinuity
        oceans.lat = horzcat(oceans.lat(:,r+1:end),oceans.lat(:,1:r));
        oceans.id = horzcat(oceans.id(:,r+1:end),oceans.id(:,1:r));
        oceans.id = oceans.id + 1;oceans.id(oceans.id~=1) = -1;                 % Switch ocean/continent = 1/0 to -1/1
        DataID = interp2(oceans.lon,oceans.lat,oceans.id,LonN,LatN);            % Interpolate ID
        DataID(DataID>0) = 1;DataID(DataID<=0) = -1;                            % Negative values for oceans, positive for continent
        DataID(LatN<=-85.95) = 1;                                               % fixed continent = antarctica
    else 
        DataID = [];
    end
end
end

