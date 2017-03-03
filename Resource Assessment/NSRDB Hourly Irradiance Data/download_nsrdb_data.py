"""This code can download NSRDB data via the NREL API. However, the API is limited to
50 cells per hour or 100 cells per day, per API key or IP address. So it is easier
to download larger regions via the NSRDB viewer at https://nsrdb.nrel.gov/nsrdb-viewer 
or https://maps.nrel.gov/nsrdb-viewer/. 

To use that, click on the "Download Data" tab at the top, left corner, then on "NSRDB Data Download (Box)" 
(click the box icon). Then follow the prompts. Use the PSM dataset, not Spectral TMY. Turn off "half hour intervals".
Turn on GHI and wind direction (not really needed, but could be nice to have).
Then follow the prompts from there. This will allow download of up to about 5 years of hourly data for the 
island of Oahu (232 cells). It takes 1-3 days before they send a download link. Then you can download the 
data via Globus.
"""

import urllib, os
import numpy as np

def make_nsrdb_url(lat, lon, year):
    # boilerplate from https://nsrdb.nrel.gov/api-instructions
    # Declare all variables as strings. Spaces must be replaced with '+', i.e., change 'John Smith' to 'John+Smith'.
    # Define the lat, long of the location and the year
    # You must request an NSRDB api key from the link above
    # api_key = 'qAk6u4D30wTNrQvNZ4eVI7TdioiPofQbvOOza1cq'
    # your_email = 'mfripp@hawaii.edu'

    api_key = 'wiYIiWfFDwQUwbIAbsdpXO07y6lGNFbg9Hn4rfUF'
    your_email = 'mfripp@gmail.com'

    #api_key = 'X7gJD0EDpWOZMVvYEzzZVWdNrGV9Gk0O6GXMh0nH'
    #your_email = 'dummy@dummy.com'

    # Set the attributes to extract (e.g., dhi, ghi, etc.), separated by commas.
    # attributes = 'ghi,dhi,dni,wind_speed_10m_nwp,surface_air_temperature_nwp,solar_zenith_angle'

    # Set leap year to true or false. True will return leap day data if present, false will not.
    leap_year = 'false'
    # Set time interval in minutes, i.e., '30' is half hour intervals. Valid intervals are 30 & 60.
    interval = '60'
    # Specify Coordinated Universal Time (UTC), 'true' will use UTC, 'false' will use the local time zone of the data.
    # NOTE: In order to use the NSRDB data in SAM, you must specify UTC as 'false'. SAM requires the data to be in the
    # local time zone.
    utc = 'false'
    # Your full name, use '+' instead of spaces.
    your_name = 'Matthias+Fripp'
    # Your reason for using the NSRDB.
    reason_for_use = 'power+system+planning'
    # Your affiliation
    your_affiliation = 'University+of+Hawaii'
    # Please join our mailing list so we can keep you up-to-date on new developments.
    mailing_list = 'false'

    # Declare url string
    # url = (
    #     'http://developer.nrel.gov/api/solar/nsrdb_0512_download.csv?'
    #     'wkt=POINT({lon}%20{lat})&names={year}&leap_day={leap}&interval={interval}&utc={utc}'
    #     '&full_name={name}&email={email}&affiliation={affiliation}&mailing_list={mailing_list}'
    #     '&reason={reason}&api_key={api}&attributes={attr}'.format(
    #         year=year, lat=lat, lon=lon, leap=leap_year, interval=interval,
    #         utc=utc, name=your_name, email=your_email, mailing_list=mailing_list,
    #         affiliation=your_affiliation, reason=reason_for_use, api=api_key, attr=attributes
    #     )
    # )
    url = (
        'http://developer.nrel.gov/api/solar/nsrdb_0512_download.csv?'
        'wkt=POINT({lon}%20{lat})&names={year}&leap_day={leap}&interval={interval}&utc={utc}'
        '&full_name={name}&email={email}&affiliation={affiliation}&mailing_list={mailing_list}'
        '&reason={reason}&api_key={api}'.format(
            year=year, lat=lat, lon=lon, leap=leap_year, interval=interval, 
            utc=utc, name=your_name, email=your_email, mailing_list=mailing_list, 
            affiliation=your_affiliation, reason=reason_for_use, api=api_key
        )
    )
    return url
         
for lon in np.arange(-158.3, -157.66+.04, 0.04):
    for lat in np.arange(21.25, 21.73+.04, 0.04):
        for year in [2007]:
            outfile = os.path.join('nsrdb', 'nsrdb_{:.3f}_{:.3f}_{}.csv'.format(lat, lon, year))
            if os.path.exists(outfile) and os.path.getsize(outfile) > 200:
                print "skipping {}, already downloaded".format(outfile)
            elif (lon, lat) <= (-158.14+.001, 21.53+.001):
                print "skipping {}, already downloaded on other computer".format(outfile)
            else:
                urllib.urlretrieve(make_nsrdb_url(lat, lon, year), outfile)
                print "saved {}".format(outfile)

# could look at first 2 lines to get new cell centers, but we don't
# # Return just the first 2 lines to get metadata:
# info = pd.read_csv(url, nrows=1)
# # See metadata for specified properties, e.g., timezone and elevation
# timezone, elevation = info['Local Time Zone'], info['Elevation']
