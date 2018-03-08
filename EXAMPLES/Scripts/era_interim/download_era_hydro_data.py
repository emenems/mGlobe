#!/usr/bin/env python

# Use this script to download ERA Interim data needed to compute hydrological (not atmosphe)
# Make sure the home folder contains .ecmwfapirc, see https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
# To install 'ecmwfapi': sudo pip install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

# Set year/month/day to be downloaded. This will save all data into one file for conversion.
# Thus, do not request too long interval!
start_year = '2016'
start_month = '01-01'
end_year = '2016'
stop_month = '12-31'
# Set time to be downloaded (usually 12:00:00)
# Possible values: '00:00:00/06:00:00/12:00:00/18:00:00'
time = '12:00:00';
# Set output file.
output_file = 'i:/GlobalModel/ERAinterim/Surface/OriginalForHydro/ERA_SM_and_Snow.nc'


# prepare/declare arguments for downloading (will be modified later)
args_sl = {
    "class": "ei",
    "dataset": "interim",
    "date": "2015-02-01/to/2016-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "39.128/40.128/41.128/42.128/141.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
	'format': "netcdf",
}
# Download data (update argument + download)
args_sl['target'] = output_file
args_sl['time'] = time
args_sl['date'] = start_year + '-' + start_month + '/to/' + end_year + '-' + stop_month
server.retrieve(args_sl)
	
	
			
