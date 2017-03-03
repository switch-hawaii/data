#!/bin/bash
set -e
set -u
set -o pipefail

echo ""
read -p "Are you sure you want to rebuild the owits hourly data file? " -n 1
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

scriptpath=$(unset CDPATH && cd "$(dirname "$0")" && echo $PWD)
datapath="/Volumes/LaCie/OWITS_DATA"
mypipe="/tmp/owits_pipe"

echo ""
echo "Creating owits.hourly table..."
mysql -e '
create database if not exists owits;
use owits;
drop table if exists hourly;
create table hourly
  (grid char(1), i smallint, j smallint,
  datetime_utc datetime,
  temp_sfc float, press_sfc float, precip float, hum_2m float,
  dswrf float, dlwrf float,
  speed10 float, dir10 float, temp10 float,
  speed50 float, dir50 float, temp50 float,
  speed80 float, dir80 float, temp80 float,
  speed100 float, dir100 float, temp100 float,
  speed200 float, dir200 float, temp200 float
  );
'
# fields in the file:
# DATE - Date
# TIME - Time (Greenwich Mean Time)
# TSFC - Surface Skin Temperature (K)
# PSFC - Surface Pressure (mb)
# PCP  - Accumulation Precipitation (mm or kg/m^2)
# Q2M  - Specific Humidity at 2M Above Ground Level (g/kg)
# DSWRF - Downward Shortwave Radiation Flux (W/m^2)
# DLWRF - Downward Longwave Radiation Flux (W/m^2)
# T10 - Temperature at 10M Above Ground Level (K)
# S10 - Wind Speed at 10M Above Ground Level (m/s)
# W10 - Wind Direction at 10M Above Ground Level (Degrees)
# T50 - Temperature at 50M Above Ground Level (K)
# S50 - Wind Speed at 50M Above Ground Level (m/s)
# W50 - Wind Direction at 50M Above Ground Level (Degrees)
# T80 - Temperature at 80M Above Ground Level (K)
# S80 - Wind Speed at 80M Above Ground Level (m/s)
# W80 - Wind Direction at 80M Above Ground Level (Degrees)
# T100 - Temperature at 100M Above Ground Level (K)
# S100 - Wind Speed at 100M Above Ground Level (m/s)
# W100 - Wind Direction at 100M Above Ground Level (Degrees)
# T200 - Temperature at 200M Above Ground Level (K)
# S200 - Wind Speed at 200M Above Ground Level (m/s)
# W200 - Wind Direction at 200M Above Ground Level (Degrees)

if [[ ! -p $mypipe ]]; then
    mkfifo -m 0666 $mypipe
fi
find "$datapath" -type f \
  -exec "$scriptpath/import_owits_hourly_file.sh" "{}" "$mypipe" \;
