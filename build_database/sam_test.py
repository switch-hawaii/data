from sscapi import *

def setup_pv(ssc, data):
	ssc.data_set_number( data, 'system_capacity', 1 )
	ssc.data_set_number( data, 'module_type', 0 )   # 0=standard, 1=premium, 2=thin film
	ssc.data_set_number( data, 'array_type', 3 )    # single-axis, backtracked
	ssc.data_set_number( data, 'losses', 14 )   # given in % (is 14% standard?)
	ssc.data_set_number( data, 'tilt', 0 )  # 0 for tracking I think
	ssc.data_set_number( data, 'azimuth', 180 )
	ssc.data_set_number( data, 'adjust:constant', 0 )
	ssc.data_set_number( data, 'dc_ac_ratio', 1.25 )

def run_pvwattsv5( ssc, data ):
	# run PV system simulation
	mod = ssc.module_create("pvwattsv5")	
	ssc.module_exec_set_print( 0 );
	if ssc.module_exec(mod, data) == 0:
		print 'PVWatts V5 simulation error'
		idx = 1
		msg = ssc.module_log(mod, 0)
		while (msg != None):
			print '\t: ' + msg
			msg = ssc.module_log(mod, idx)
			idx = idx + 1
	else:
		ann = ssc.data_get_number(data, "ac_annual")
		print 'PVWatts V5 Simulation ok, e_net (annual kW)=', ann
	ssc.module_free(mod)

wf = '/Users/Matthias/SAM Downloaded Weather Files/lat21.55620_lon-158.03720_2007.csv';

print wf

ssc = PySSC()
dat = ssc.data_create()
setup_pv(ssc,dat);

######## Simple Setup ########
ssc.data_set_string( dat, 'solar_resource_file', wf );
run_pvwattsv5( ssc, dat );

######## Detailed Setup ########
ssc.data_clear(dat);

# read a weather file for this example program
# and extract the data from it into a bunch of Python variables
# note: this weather data could come from any source
ssc.data_set_string( dat, 'file_name', wf );
ssc.module_exec_simple_no_thread( 'wfreader', dat );	
lat = ssc.data_get_number(dat, 'lat');
lon = ssc.data_get_number(dat, 'lon');
tz = ssc.data_get_number(dat, 'tz');
elev = ssc.data_get_number(dat, 'elev');
year = ssc.data_get_array(dat, 'year');
month = ssc.data_get_array(dat, 'month')
day = ssc.data_get_array(dat, 'day');
hour = ssc.data_get_array(dat, 'hour');
minute = ssc.data_get_array(dat, 'minute');
beam = ssc.data_get_array(dat, 'beam');
diffuse = ssc.data_get_array(dat, 'diffuse');
wspd = ssc.data_get_array(dat, 'wspd');
tdry = ssc.data_get_array(dat, 'tdry');
albedo = ssc.data_get_array(dat, 'albedo');	

# create an SSC data with a bunch of fields
ssc.data_clear( dat );

wfd = ssc.data_create();
ssc.data_set_number( wfd, 'lat', lat);
ssc.data_set_number( wfd, 'lon', lon);
ssc.data_set_number( wfd, 'tz',  tz);
ssc.data_set_number( wfd, 'elev',  elev);

ssc.data_set_array( wfd, 'year',  year);
ssc.data_set_array( wfd, 'month',  month);
ssc.data_set_array( wfd, 'day',  day);
ssc.data_set_array( wfd, 'hour', hour);
ssc.data_set_array( wfd, 'minute', minute);
# note: if using an hourly TMY file with integrated/averaged
# values, do not set the minute column here. otherwise
# SSC will assume it is instantaneous data and will not adjust
# the sun position in sunrise and sunset hours appropriately
# however, if using subhourly data or instantaneous NSRDB data
# do explicitly provide the minute data column for sunpos calcs

ssc.data_set_array( wfd, 'dn', beam);
ssc.data_set_array( wfd, 'df', diffuse);
ssc.data_set_array( wfd, 'wspd', wspd);
ssc.data_set_array( wfd, 'tdry', tdry);
ssc.data_set_array( wfd, 'albedo', albedo);

# instead of setting a string weather file, simply
# set the table variable that contains the various fields
# with solar resource data
ssc.data_set_table( dat, 'solar_resource_data', wfd );	

# we can free the resource data table now, since
# the previous line copies it all into SSC
ssc.data_free(wfd);

# set up other PV parameters and run
setup_pv( ssc,dat );
run_pvwattsv5( ssc, dat );

ssc.data_free(dat);
