#!/bin/bash
set -e
set -u
set -o pipefail

# this script averages and inserts one owits data file into the mysql owits.annual table
# the first argument is the name of the file to read
# the second argument is the name of a fifo (named pipe) to use to pass the data to mysql
# (this should be created in advance)

# you can call this script repeatedly with something like this:
# mkfifo /tmp/owits_pipe
# find . -type f -exec import_owits_file.sh "{}" /tmp/owits_pipe \;

if [ $# -ne 2 ]
then
  echo "Usage: `basename $0` {data file} {existing fifo}"
  exit 65
fi

source_file=$1
pipe=$2

filename=`basename $source_file`  # e.g., 0005_0005.HAWAII.E.txt
i=${filename:0:4}
j=${filename:5:4}
grid=${filename:17:1}

#echo "Importing i=$i, j=$j, grid=$grid from $source_file"
#exit

case $source_file in 
*.txt ) 
    echo "Loading hourly averages for cell ($grid, $i, $j) from file $source_file"
    # stuff the file into the pipe as it is
    cat $source_file \
      | awk -F, '{for(i=3; i<=NF; i++) sum[i]+=$i } 
                 (NR%6)==0 {
                   printf("%s,%s,", $1,$2); 
                   for( i=3; i<=NF; i++) {
                     printf("%.2f%s", sum[i] / 6, (i==NF) ? "\n" : FS);
                     sum[i]=0;
                   }
                 }' \
      > $pipe &
    ;; 
*.txt.gz ) 
    echo "Loading hourly averages for cell ($grid, $i, $j) from compressed file $source_file"
    # uncompress the file into the pipe
    gzip --stdout -d $source_file \
      | awk -F, '{for(i=3; i<=NF; i++) sum[i]+=$i } 
                 (NR%6)==0 {
                   printf("%s,%s,", $1,$2); 
                   for( i=3; i<=NF; i++) {
                     printf("%.2f%s", sum[i] / 6, (i==NF) ? "\n" : FS);
                     sum[i]=0;
                   }
                 }' \
      > $pipe &
    ;; 
*)
    echo "skipping file $source_file (unknown format)"
    exit
esac

# got a file to read

mysql <<SQLCOMMANDS
load data local infile "$pipe"
  into table owits.hourly
  fields terminated by "," lines terminated by "\n"
  (@DATE, @TIME, 
    temp_sfc, press_sfc, precip, hum_2m,
    dswrf, dlwrf,
    temp10, speed10, dir10, 
    temp50, speed50, dir50, 
    temp80, speed80, dir80, 
    temp100, speed100, dir100, 
    temp200, speed200, dir200)
  set grid="$grid", i=$i, j=$j, 
    datetime_utc=date_sub(str_to_date(concat(@DATE, " ", @TIME), "%Y%m%d %H%i"), interval 50 minute);
SQLCOMMANDS


