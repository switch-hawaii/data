#!/bin/bash
set -e
set -u
set -o pipefail

echo ""
read -p "Are you sure you want to add data to the owits annual data file? " -n 1
if [[ ! $REPLY =~ ^[Yy]$ ]]
then
    exit 1
fi

scriptpath=$(unset CDPATH && cd "$(dirname "$0")" && echo $PWD)
datapath="/Volumes/LaCie/OWITS_DATA"
mypipe="/tmp/owits_pipe"

#echo ""
#echo "Creating owits.annual table..."
#mysql -e '
#create database if not exists owits;
#use owits;
#drop table if exists annual;
#create table annual
#  (grid char(1), i smallint, j smallint,
#  dswrf float, 
#  speed10 float, speed50 float, speed80 float, speed100 float, speed200 float, 
#  dir100 float
#  );
#'

if [[ ! -p $mypipe ]]; then
    mkfifo -m 0666 $mypipe
fi
find "$datapath/F" "$datapath/G" -type f \
  -exec "$scriptpath/average_and_import_owits_file.sh" "{}" "$mypipe" \;

# read in the georeference data
mysql -e'
  use owits;
  drop table if exists cells;
  create table cells (grid char(1), i smallint, j smallint, lat float, lon float);
  alter table cells add index gjill (grid, i, j, lat, lon);
'

echo "loading cell descriptions into owits.cells"
for GRID in D E F G
do
  mysql -v -v -e'
    load data local infile 
      "'"${datapath}/../OWITS/${GRID}_Georef.csv"'"
    into table owits.cells
    fields terminated by "," lines terminated by "\r\n"
    ignore 1 lines
    (i, j, lat, lon)
    set grid="'$GRID'";
'
done