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
    echo "Averaging and loading cell ($grid, $i, $j) from file $source_file"
    # stuff the file into the pipe as it is
    cat $source_file \
      | awk -F, '{for(i=1; i<=NF; i++) sum[i]+=$i } END { for( i=1; i<=NF; i++) printf("%.2f%s", sum[i] / FNR, (i==NF) ? "\n" : FS) }' \
      > $pipe &
    ;; 
*.txt.gz ) 
    echo "Averaging and loading cell ($grid, $i, $j) from compressed file $source_file"
    # uncompress the file into the pipe
    gzip --stdout -d $source_file \
      | awk -F, '{for(i=1; i<=NF; i++) sum[i]+=$i } END { for( i=1; i<=NF; i++) printf("%.2f%s", sum[i] / FNR, (i==NF) ? "\n" : FS) }' \
      > $pipe &
    ;; 
*)
    echo "skipping file $source_file (unknown format)"
    exit
esac

# got a file to read

mysql <<SQLCOMMANDS
load data local infile "$pipe"
  into table owits.annual
  fields terminated by "," lines terminated by "\n"
  (@DATE, @TIME, @TSFC, @PSFC, @PCP, @Q2M, dswrf, @DLWRF, @T10, speed10, @W10, @T50, speed50, @W50, @T80, speed80, @W80, @T100, speed100, dir100, @T200, speed200, @W200)
  set grid="$grid", i=$i, j=$j;
SQLCOMMANDS
