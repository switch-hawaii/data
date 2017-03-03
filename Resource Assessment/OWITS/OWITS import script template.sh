#!/bin/bash

# this script inserts one OWITS data file into the postgresql database
# you can call it repeatedly with something like this:
# find . -type f -exec ../fripp/test.sh "{}" \;
# (called from a bash prompt in the OWITS_DATA directory)

#!/bin/bash

csv_file=$1
filename=`basename $csv_file`
i=${filename:0:4}
j=${filename:5:4}
grid=${filename:17:1}

echo "loading grid $grid, i=$i, j=$j from file $csv_file"

psql OWITS <<SQLCOMMANDS
create temporary table timport like data10min;
copy timport ( cell list ) from $csv_file with (format CSV, HEADER);
update timport set grid_id="$grid", i=$i, j=$j;
insert into data10min select * from timport;
SQLCOMMANDS
