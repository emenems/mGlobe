#!/usr/bin/env python

# Use this script to download ERA Interim data needed to compute atmospheric effect (not hydrological)
# Make sure the home folder contains .ecmwfapirc, see https://software.ecmwf.int/wiki/display/WEBAPI/Access+ECMWF+Public+Datasets
# To install 'ecmwfapi': sudo pip install https://software.ecmwf.int/wiki/download/attachments/56664858/ecmwf-api-client-python.tgz
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

# Set years to be downloaded. 
start_year = 2003
end_year = 2003
# If you do not wish to download whole year (not recommended,see mGlobe user manual), 
# modify the start_month and end_month. This will be applied to each year!
start_month = '01-01'
stop_month = '12-31'

# Set output folder. WARNING: It is assumed that this folder contains:
#	Geopotential
#	SpecificHumidity
#	Temperature
#	Surface
# If not, modify the code below (see first for loop)
output_folder = 'd:/GlobalModel/ERAinterim/'

# Set parameter to be downloaded: 
#	geopotential = 129.128, 
#	specific humidity = 133.128, 
#	temperature = 130.128
param_id = [129.128,133.128,130.128];

# prepare/declare arguments for downloading (will be modified later)
# args_pl = pressure level
# args_sl = surface level
args_pl = {
    "class": "ei",
    "dataset": "interim",
    "date": "2016-01-01/to/2016-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levelist": "1/2/3/5/7/10/20/30/50/70/100/125/150/175/200/225/250/300/350/400/450/500/550/600/650/700/750/775/800/825/850/875/900/925/950/975/1000",
    "levtype": "pl",
    "param": "129.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00",
    "type": "an",
	'format': "netcdf",
    "target": "target_name_will_be_overwritten.nc",
}
args_sl = {
    "class": "ei",
    "dataset": "interim",
    "date": "2016-01-01/to/2016-12-31",
    "expver": "1",
    "grid": "0.75/0.75",
    "levtype": "sfc",
    "param": "134.128/167.128/168.128",
    "step": "0",
    "stream": "oper",
    "time": "00:00:00/06:00:00/12:00:00/18:00:00",
    "type": "an",
	'format': "netcdf",
    "target": "target_name_will_be_overwritten.nc",
}

# Main loop for downloading all data (the range function does not return the last value => end_year+1)
for year in range(start_year, end_year+1):
	# Download surface date (all hours and parameters in one file)
	args_sl['target'] = output_folder + '/Surface/ERA_SP_2T_2D_ALLh_06hStep_D_' + str(year) + '.nc'
	args_sl['date'] = str(year) + '-' + start_month + '/to/' + str(year) + '-' + stop_month
	server.retrieve(args_sl)
	# Run loop for all parameters in pressure level
	for para in param_id:
		args_pl['param'] = str(para)
		if para == 129.128:
			output_name_pl = output_folder + '/Geopotential/ERA_GEOPOT_24hStep_'
		elif para == 133.128:
			output_name_pl = output_folder + '/SpecificHumidity/ERA_SHUMID_24hStep_'
		else:
			output_name_pl = output_folder + '/Temperature/ERA_TEMP_24hStep_'	
		args_pl['date'] = str(year) + '-' + start_month + '/to/' + str(year) + '-' + stop_month
		
		# Run loop for all hours (separate hours in each file and for each parameter)
		for hour in range(0,24,6):
			args_pl['target'] = output_name_pl + str(hour).rjust(2,'0') + '_' + str(year) + '.nc'
			args_pl['time'] = str(hour).rjust(2,'0') + ':00:00'
			server.retrieve(args_pl)
	
	
			
