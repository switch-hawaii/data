#%%#################################
"""
Screening rules for utility-scale solar projects:
- only zoned for agriculture or country
- no Class A agricultural land
- not golf courses
- not within 50 meters of street centerlines (also filters out rural businesses
  and neighborhoods)
- not steeper than 10% slope
- allowed patch is at least 60 meters in all directions

Steps to implement these rules (eventually preparing allowed_solar.tif geotiff
and solar_cluster_land_class_nsrdb_grid_final.csv text file):

May 2021:

Download the following shape files from
http://planning.hawaii.gov/gis/download-gis-data-expanded/: 
- "Zoning, C&C of Honolulu" (cty_zoning_oah.shp.zip) 
- "Land Study Bureau Classification (LSB)" (lsb.shp.zip)
- "Roads – C&C of Honolulu" (streets_oah.shp.zip)
- "Golf Courses" (golf_courses.shp.zip)

Download 1/3 arc-second elevation data for Oahu area (2 1-degree-square gpkg tiles):
- ned.usgs.gov (redirects) > 3D Elevation Program > Access Data > Get GIS Data > Cloud Browse
  (gets you to http://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/)
- Then Elevation > 13 > TIFF
- Then n22w158 > n22w158.tif or n22w158 > n22w159.tif
    - Note: .gpkg files are available this way, and are 6x smaller, but gpkg is not recognized by rasterio
- Or just go to http://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/13/TIFF/n22w158/
  and http://prd-tnm.s3.amazonaws.com/index.html?prefix=StagedProducts/Elevation/13/TIFF/n22w159/
Alternative option (TIFs):
- go to http://nationalmap.gov/ (redirects to https://www.usgs.gov/core-science-systems/national-geospatial-program/national-map)
- click [Let's make MY MAP!](https://apps.nationalmap.gov/viewer) under "The National Map Viewer"
- click "[Data Download](https://apps.nationalmap.gov/downloader/#/)"
- click the check box next to "Elevation Products (3DEP)", then make sure "1/3 arc-second DEM" and "1 x 1 degree" are selected
- optionally click Show Map Index > Map Indices - 1 Degree (upper right corner)
- zoom in close on Oahu (so two Oahu tiles are in scope)
- click on "Search Products" (upper left corner). This opens the products tab with matching products.
- hover over each product to see which ones are relevant; map will darken if they are in the map
- right-click on "Download TIF" next to relevant tiles ("USGS 13 arc-second n22w158 1 x 1 degree" and 
  "USGS 13 arc-second n22w159 1 x 1 degree") and download
Another option (TIFs):
- Go to https://apps.nationalmap.gov/viewer/ (various paths to get there, e.g., ned.usgs.gov (redirect) > "Let's make MY MAP!")
- Optional, seems to have no effect on results:
    - wait for top toolbar to load
    - click on layers icon in top toolbar (opens layers pane on right)
    - click on 3DEP Elevation - Slope Map in layers pane
        - I can't see any way to actually download the slope data; the options below just 
          get standard elevation data
- zoom close on Oahu
- click small gray triangle at bottom of window to "Open Attribute Table" (pane of available data at the bottom)
- find the two "USGS_13_*" products, scroll right, right-click on URL to download TIFs

Midway through this work, draw polygons around clusters of solar-suitable cells;
assign each a unique ID; save as solar_clusters.shp.zip

"""

"""
API Notes:

This uses GRASS instead of geopandas, rasterio gdal and rasterstats because this
provides all the required behavior in a single environment, e.g., intersecting
vector layers, mosaicing rasters, calculating zonal statistics,
buffering/growing raster regions (haven't found a way to do that without GRASS).

This uses grass.script instead of grass.pygrass because it is possible to
discover and test commands and arguments for grass.script in the GUI grass
environment, then copy them over to the script.

This runs as an external script instead of running as a module inside the GRASS
GUI because that makes it possible to run it a few lines at a time in a Jupyter
environment (e.g., Visual Studio Code) then check the results in GRASS and
figure out settings for the next step. GRASS normally runs Python scripts
directly (see
https://grasswiki.osgeo.org/wiki/Working_with_GRASS_without_starting_it_explicitly),
but the internal editor is limited and it is not easy to run a few lines at a
time.

Before running this, you must install a copy of GRASS GIS and also install
grass-session (pip install grass-session, not available on conda) in your
working environment.

see https://grasswiki.osgeo.org/wiki/Working_with_GRASS_without_starting_it_explicitly
for info on accessing GRASS features from an external program

NOTE: once the first block of code (setup) is run, all the grass commands are in
the search path, so you can test them in the jupyter window (lead with !) or
get help on each one by running, e.g. !v.overlay (grass will give a list of arguments).
"""


##################################
# setup

config = {
    # dir where temporary grass files will be created
    # (you can use the "Add existing or create new database" button
    # on the Data tab in the GRASS GUI to open this dir and see files as
    # they are created)
    # NOTE: this will be deleted if it exists already
    'grass_data_dir': 'grassdata',
    # directory where elevation files are stored
    'dem_dir': 'USGS',
    # crs used for elevation files (change this if USGS changes it).
    # note: if the dem_crs switches to non-geographic coordinates, we
    # also need to change the average latitude and longitude calculations 
    # later
    'dem_crs': 'EPSG:4269',  # NAD83, geographic
    # directory where planning shape files are stored
    'planning_dir': 'planning.hawaii.gov',
    # CRS to use when distances need to be in meters (UTM 4N for Hawaii)
    'utm_crs': 'EPSG:26904',  # NAD83, UTM 4N
    'census_dir': 'census',
}

import glob, pathlib, shutil, shlex, re, time, os

# set some standard environment variables for the way we will use GRASS
os.environ.update(dict(
    GRASS_OVERWRITE='1',
    GRASS_COMPRESS_NULLS='1',
    GRASS_MESSAGE_FORMAT='plain',
    GRASS_VERBOSE='0',
))

# set GRASSBIN environment variable, needed by grass_session
if not os.getenv('GRASSBIN'):
    grass_bin_file = sorted(
        d 
        for d in (
            glob.glob('/Applications/GRASS-*.app/Contents/Resources/bin/grass*')
            + glob.glob('/Applications/GRASS-*.app/Contents/MacOS/grass*')
        ) if not d.endswith('.sh')  # avoid shell scripts with similar names
    )[-1]
    os.environ['GRASSBIN'] = grass_bin_file

# unset PROJ_LIB environment variable (sometimes set by other GIS tools),
# which can confuse GRASS's projection lookup
os.environ.pop('PROJ_LIB', '')

from grass_session import Session
import grass.script as grass

# enable display of error messages on Jupyter
grass.set_capture_stderr(True)

# include {os.getpid()} if you need to run on multiple machines
# without conflicts
temp_map = f'temp1'
temp_map_2 = f'temp2'

# helper functions to run commands
# def run(cmd, **kwargs):
#     """
#     Run cmd without setting --overwrite flag.
#     This is useful for commands that modify environment
#     """
#     return run_command(cmd, **kwargs)

def process_backspace(string):
    while '\x08' in string:
        string = re.sub('[^\x08]\x08', '', string)
        string = re.sub('^\x08', '', string) # drop any from start of string
    return string

def run(*args, **kwargs):   
    """
    Run grass command. This is equivalent to grass.read_command, except that it
    reports any warnings that were sent to stderr.
    """
    print("==========================================")
    print('Running', shlex.join(grass.make_command(*args, **kwargs)))

    # based on read2_command at https://grasswiki.osgeo.org/wiki/GRASS_Python_Scripting_Library
    ps = grass.start_command(*args, stdout=grass.PIPE, stderr=grass.PIPE, **kwargs)
    stdout, stderr = ps.communicate()
    returncode = ps.poll()

    # show any warnings that occurred
    print(stderr.decode())

    if returncode:
        raise RuntimeError('GRASS returned non-zero result')

    return stdout.decode()  # return any stdout info

def planning_file(stem):
    """ Return the name of a .shp.zip file in the planning dir. """
    return os.path.join(config['planning_dir'], stem + '.shp.zip')

# open dem location (in case next cell is skipped)
session = Session()
if os.path.exists(config['grass_data_dir']):
    session.open(gisdb=config['grass_data_dir'], location='dem')

#%%###########################
# delete and recreate dem and utm locations
# Note: r.in.gdal could be used to create a location from one of the raster
# files, but we would first have to create a dummy location with an arbitrary
# projection, because grass_session doesn't seem to have a way to create a
# database without also creating a location.
if os.path.exists(config['grass_data_dir']):
    print("==============================================================")
    print(f"WARNING: deleting existing files from: {config['grass_data_dir']}")
    print("Press ctrl-c or Jupyter 'stop' button (🔲) in the next 10 seconds to interrupt.")
    print("==============================================================")
    for i in range(10, 0, -1):
        print(f"{i}...", end="")
        time.sleep(1)
    print("continuing.")
    shutil.rmtree(config['grass_data_dir'])
session.open(gisdb=config['grass_data_dir'], location='dem', create_opts=config['dem_crs'])
# in principle we should close sessions before switching, but grass_session
# sets path variables only when creating a session (not during open), then 
# clears them from .close(), so then we can't find GRASS again. So we just
# open new sessions as needed, which just leaves a few config files in the 
# temp dir.
session.open(gisdb=config['grass_data_dir'], location='utm', create_opts=config['utm_crs'])
# open dem location for main workflow
session.open(gisdb=config['grass_data_dir'], location='dem')

#%%#################################
# Import and merge 1-degree elevation tiles, then calculate slope

# import all usgs*.tif files in dem_dir
el_files = glob.glob(os.path.join(config['dem_dir'], 'USGS_*.tif'))
el_maps = []
for f in el_files:
    map = pathlib.Path(f).stem  # filename without path or extension
    el_maps.append(map)
    run(
        'r.in.gdal', 
        input=f, 
        output=map, 
        overwrite=True
    )
    # run(f"r.in.gdal input={f} output={map}")

# set resolution and extent to match elevation files
run('g.region', raster=el_maps)
# run(f"g.region raster={','.join(el_maps)}")

# Crop region more narrowly (optional, found by inspection of elevation files in GRASS GUI).
# This clips the overall elevation tile calculated below.
# TODO: maybe use Oahu shape from counties shapefile to set region
run('g.region', n="21:45:00N", s="21:14:30N", e='157:38:30W', w="158:19:00W")

# merge into one tile (using region's resolution and extent)
run('r.patch', overwrite=True, input=el_maps, output='elevation')

# calculate slope for each cell
run(
    'r.slope.aspect', 
    elevation='elevation', slope='slope', 
    format='percent',
    overwrite=True
)

#%%#################################
# Import/define various vector layers that will be useful later
# zoning
# load_zone
# land_class
# golf_courses
# streets
# street_buffers

# TODO: update these queries to cover multiple islands, not just Oahu

# Note: we use a snap of 1e-4 (10 meters) during import, which does a pretty
# good job of cleaning up topological errors. There are a lot in these shape
# files, and without the snap, regions tend to drop out of some of the files,
# e.g., the zoning file, as they are processed. However, some errors are still
# left. These errors and possible solutions are discussed below the load_zone 
# creation step below. For now, we make do with rough but workable topologies.

# zoning
run(
    'v.import',
    input=planning_file('cty_zoning_oah'),
    output='zoning',
    snap=1e-4,
)
run('v.db.addcolumn', map='zoning', columns='load_zone VARCHAR')
run('v.db.update', map='zoning', column='load_zone', value='Oahu')

# load_zone (all land on each island); this assumes
# the zoning file has all islands
run('v.dissolve', input='zoning', output='load_zone', column='load_zone')

# check topology (run these in the console in the GRASS gui)
# v.info load_zone
# d.vect load_zone
# d.vect load_zone type=centroid color=red size=10
# g.gui.vdigit zoning
#    in digitizer, zoom in on problematic areas; change settings to display all vertices;
#    try moving some vertices and lines to see how the topology is arranged (see grass 
#    help for how this works, but the main idea is click with left button to pickup
#    feature, then click with right button to put it down.
# Areas seem to be basically right, but the vector file seems to have a lot of
# areas with no centroids and a lot of islands. The on-screen display also shows
# a lot of stray internal boundaries, which seem to be areas that were
# overlapping in the original zoning file. e.g., near the southern edge of Oahu,
# there is a curved border between resort and business zoning, that was
# digitized with vertices every 3 meter (roughly) on one side and every 0.3
# meter (roughly) on the other side. Snapping can't fix these. In the future we
# could try v.clean to prune these, then re-snap.

# TODO: clean up shapefile topologies a little more after import, so the 
# dissolve comes out cleaner.
# One option would be to use ogr2ogr, as suggested here:
# http://osgeo-org.1560.x6.nabble.com/large-shapefile-not-importing-properly-with-v-import-td5339299.html
# ogr2ogr should be in the same path as r.import, etc. e.g., run 
# f"ogr2ogr -t_srs {config['dem_crs']} {planning_file('cty_zoning_oah_ogr2ogr')} {planning_file('cty_zoning_oah')}"
# But ogr2ogr may be less effective than v.import snap=1e-4 that we use already
# A better option might be to arbitrarily remove any overlapping areas, following advice
# here:
# https://grasswiki.osgeo.org/wiki/Vector_topology_cleaning
# https://grass.osgeo.org/grass78/manuals/v.clean.html
# https://grasswiki.osgeo.org/wiki/Vector_Overlapping_Areas
# https://gis.stackexchange.com/a/170825
# The last one seems about right: each feature has multiple categories assigned
# on each layer (layers seem to be have many-to-many relations between feature IDs 
# and cats), and you need to use v.edit tool=catdel to delete one of them,
# for which you need to specify the feature id and the cat (in a where clause).
# Is there some way to do this via db.execute? Then dissolve on the remaining 
# category and it should be pretty clean.
# Another option could be to use v.clean to prune extra points, which might make
# snapping more effective, then maybe go in to inspect the remaining problems,
# using the methods shown in v.clean.

# land_class
# Shows land_class and load_zone for all land in the study area.
# This is based on Land Study Bureau classifications, but only 
# labels class A, B or C or "other" (everything else).
run(
    'v.import',
    input=planning_file('lsb'),
    output='lsb',
    snap=1e-4,
)
# limit to only the classes that could be restricted
run(
    'v.extract',
    input='lsb',
    output=temp_map,
    where="type IN ('A', 'B', 'C')"
)
# Dissolve each class together, so we only get one polygon per load_zone -
# land_class combination later. This works temporarily, but ends up getting
# broken back into the constituent areas during the overlay with load_zone,
# so it ends up having no effect. (See note below.)
run(
    'v.dissolve',
    input=temp_map,
    output=temp_map_2,
    column='type'
)
# add in "other" land
run(
    'v.overlay',
    ainput=temp_map_2,
    binput='load_zone',
    output='land_class',
    operator='or',
)
# tidy up columns
run('v.db.dropcolumn', map='land_class', columns='a_cat,b_cat')
run('v.db.renamecolumn', map='land_class', column='a_type,land_class')
run('v.db.renamecolumn', map='land_class', column='b_load_zone,load_zone')
# fill in missing/unrestricted land_class as 'other'
run(
    'v.db.update', 
    map='land_class', 
    column='land_class', 
    where='land_class IS NULL', 
    value='other'
)
# should only have a dozen or so rows, but actually has thousands
# !db.describe land_class
# !v.db.select land_class | head -n 20

# TODO: (if possible) Get code above to dissolve common land classes in each
# zone together, so there is a single polygon per load_zone-land_class
# combination. Then we could export the table directly in the next step, without
# needing to group by load_zone and land_class, as we currently do. This would
# also avoid the need to aggregate when creating solar_project later. In QGIS or
# ArcGIS, this would be done by the dissolve step shown above. But in GRASS GIS,
# doing an overlay with previously dissolved layers ends up breaking them back
# apart into their constituent areas (e.g., multiple islands, or places with
# weird internal topology left over from the import). We could in principle
# dissolve the final land_class map on load_zone and land_class at the end, but
# that's tricky to do and there's not much point, since it will all get broken
# up again when this is overlaid with nsrdb_grid later. So we just aggregate in
# SQL before exporting the data. (To dissolve on two columns, we'd probably need
# to define a field lzlc = load_zone || '_' | land_class, then dissolve on lzlc,
# then use join back to the original map using lclz, to read in the original
# attributes, which are discarded during v.dissolve.)

# Save total area in each class in each zone for use in Switch
run('v.to.db', map='land_class', option='area', column='area')
run(
    'db.select',
    sql="""
        SELECT load_zone, land_class, SUM(area) AS area 
        FROM land_class 
        WHERE load_zone='Oahu'
        GROUP BY 1, 2
    """,
    output='load_zone_ag_land_area.csv',
    separator=','
)

# golf courses
run(
    'v.import',
    input=planning_file('golf_courses'),
    output='golf_courses',
    snap=1e-4,
)

# remove 50 meter buffer from centerline of all streets
# note: this removes most neighborhoods with homes from the dataset; another
# option would be to remove all census blocks with population density (1000000 *
# population / area) > 100 people per km2, but that includes some blocks with
# large unpopulated parts

#### Switch to 'utm' location so we can buffer in meters
session.open(gisdb=config['grass_data_dir'], location='utm')
run(
    'v.import',
    input=planning_file('streets_oah'),
    output='streets',
    snap=1e-4,
)
# Next step takes surprisingly long (~8 minutes)
# Maybe it could be sped up by pre-dissolving the lines, but
# I can't see a way to do that in GRASS, and it's unlikely to
# make a difference anyway, because GRASS seems to always 
# dissolve post-buffer, so presumably it already dissolves
# pre-buffer if that would help.
now = time.strftime('%H:%M', time.localtime())
print(f"The next step takes about 8 minutes (starting at {now}):")
start = time.time()
run(
    'v.buffer',
    input='streets',
    output='streets_buffered',
    distance=50,  # in map units, which are meters for this crs
    overwrite=True
)
print(f'elapsed time: {(time.time()-start)/60} minutes')
#### Switch back to 'dem' location
session.open(gisdb=config['grass_data_dir'], location='dem')

# copy streets into main location (only buffers are needed for 
# analysis, but actual streets are useful for map making)
run(
    'v.proj',
    location='utm',
    input='streets_buffered',
    output='streets_buffered',
    overwrite=True
)
run('v.proj', location='utm', input='streets', output='streets', overwrite=True)

#%%#################################
# Combine various vector layers to define allowed_land_use map
# (used to knock out some of the allowed solar raster and possibly
# shown on maps of the area)

# only consider certain land use zoning
run(
    'v.extract', 
    flags='dt',
    input='zoning',
    where="ZONING_LAB in ('AG-1', 'AG-2', 'COUNTRY')",
    output='allowed_zoning',
    new=0,   # put them all in the same category so they will dissolve together
)
# make extra map to display the cutouts
run(
    'v.extract', 
    flags='d',
    input='zoning',
    where="ZONING_LAB NOT IN ('AG-1', 'AG-2', 'COUNTRY')",
    output='unallowed_zoning',
)

# drop unallowed land classes (only A)
# TODO: maybe handle this in the land limitations code instead of just 
# dropping it from the dataset
run(
    'v.extract', 
    flags='dt',
    input='land_class',
    where="land_class == 'A'",
    output='unallowed_land_class',
)
run(
    'v.overlay',
    ainput='allowed_zoning',
    binput='unallowed_land_class',
    operator='not',
    output='allowed_zoning_land_class',
)

# drop golf courses
run(
    'v.overlay',
    ainput='allowed_zoning_land_class',
    binput='golf_courses',
    operator='not',
    output='allowed_zoning_land_class_golf',
)

# drop areas around streets
run(
    'v.overlay',
    ainput='allowed_zoning_land_class_golf',
    binput='streets_buffered',
    operator='not',
    output='allowed_land_use',  # final output
)

#%%############################
# Create rasters of allowed regions, filtering for land use and small parcel sizes

# Clip slope raster using allowed_land_use polygon
run('r.mask', vector='allowed_land_use')
# this calculation will be restricted to the masked area (https://gis.stackexchange.com/a/359248)
run('r.mapcalc', expression='slope_allowed_land_use=slope')
# turn off mask
run('r.mask', flags='r', vector='allowed_land_use')

# Identify land with slopes below 15 percent and below 30 percent
# Later we will identify all suitable land in each, then subtract the 
# below-15 from the below-30, so small fragments of below-15 will be 
# omitted from the 15 layer, but can still be used in the 30 layer if 
# they are near other below-30 land.
slopes = [15, 30]
for s in slopes:
    run(
        'r.mapcalc', 
        expression=f"slope_below_{s}_with_fragments = if(slope_allowed_land_use <= {s}, 1, null())"
    )
    run('r.colors', map=f'slope_below_{s}_with_fragments', color='greens')
# shrink allowed boundaries by 3 pixels (30 meters) to clear out small parcels (10s)
for s in slopes:
    run(
        'r.grow', 
        input=f'slope_below_{s}_with_fragments', 
        output=f'slope_below_{s}_shrunk',
        radius=-3.01,  # a little more than 3 to make sure we hit the edges
        metric='euclidean',
    )
    run('r.colors', map=f'slope_below_{s}_shrunk', color='reds')
# expand allowed areas back to original boundaries, plus a little to
# make sure we get adjacent corners
for s in slopes:
    run(
        'r.grow', 
        input=f'slope_below_{s}_shrunk', 
        output=f'slope_below_{s}_shrunk_grown',
        radius=4.75,
        metric='euclidean',
    )
    run('r.colors', map=f'slope_below_{s}_shrunk_grown', color='blues')

# Use areas from the original raster that are within the regrown area (both non-zero, non-null)
for s in slopes:
    run(
        'r.mapcalc', 
        expression=f"allowed_00_to_{s} = if(slope_below_{s}_with_fragments * slope_below_{s}_shrunk_grown, 1, null())",
    )
    run('r.colors', map=f'allowed_00_to_{s}', color='oranges')

# identify allowed land with slopes between 15 and 30, which may be fragmented
# as long as it is adjacent to allowed land in the 0 to 15 range
run(
    'r.mapcalc',
    expression='allowed_15_to_30 = if(if(allowed_00_to_30) && isnull(allowed_00_to_15), 1, null())'
)


#%%#############################
# Calculate solar area and mean latitude and longitude for each steepness band
# in each combination of solar cluster, land class and NSRDB cell. Then export
# these in solar_project.csv for use by the import_data scripts.

#####
# make grid matching nsrdb grid (0.04 degrees horizontal and vertical)
# Grid extent: -158.32, -157.64, 21.23, 21.75 
# (cell centers: -158.3, -157.66, 21.25, 21.73)
# CRS: EPSG:4269 (NAD 83)
if config["dem_crs"] != "EPSG:4269":
    raise NotImplementedError(
        f"This program needs to be updated to use EPSG:4269 "
        f"to create the NSRDB grid; dem_crs={config['dem_crs']}."
    )
run(
    'v.mkgrid',
    map='nsrdb_grid',
    position='coor',
    grid=(13, 17),                 # number of rows, columns
    coordinates=(-158.32, 21.23),  # lower left lon, lat (outer edge of cell)
    box=(0.04, 0.04),              # width, height of each box
    type='area',
)
# record NSRDB location for use later and drop unneeded columns
run('v.to.db', map='nsrdb_grid', option='coor', columns='lon,lat')
run('v.db.dropcolumn', map='nsrdb_grid', columns='row,col,rown,coln')

######
# read solar cluster definitions (created manually beforehand by drawing
# polygons around clusters of solar-suitable cells, then assigning each a 
# unique ID; and saving as solar_clusters.shp)
run('v.import', input='solar_clusters.shp.zip', output='solar_clusters')
# assign load_zone (TODO: fill in from load_zone map)
run('v.db.addcolumn', map='solar_clusters', columns='load_zone VARCHAR')
run('v.db.update', map='solar_clusters', column='load_zone', value='Oahu')

# NOTE: We could just report the total solar area in each combination of
# load_zone and nsrdb_grid cell, then cluster those automatically in the
# utility-scale import code. But this approach, using manual clusters, works
# reasonably well.

######
# break solar clusters into land classes
run(  
    'v.overlay',
    ainput='solar_clusters',
    binput='land_class',
    operator='and',
    output='solar_cluster_land_class',
)
# rename the useful columns and drop others
run(
    'v.db.dropcolumn', 
    map='solar_cluster_land_class', 
    columns='a_cat,b_cat,b_load_zone,b_area'
)
for pair in [
    'a_cluster_id,cluster_id', 
    'a_load_zone,load_zone', 
    'b_land_class,land_class'
]:
    run('v.db.renamecolumn', map='solar_cluster_land_class', column=pair)

######
# Make solar_cluster_land_class_nsrdb_grid 
# (intersection of solar_cluster, land_class and nsrdb_grid)
run(
    'v.overlay',
    ainput='solar_cluster_land_class',
    binput='nsrdb_grid',
    output='solar_cluster_land_class_nsrdb_grid',
    operator='and'
)

# rename a_ and b_ columns for later use or drop if unneeded
rename_cols = [
    'a_cluster_id,cluster_id',
    'a_load_zone,load_zone',
    'a_land_class,land_class',
    'b_cat,nsrdb_id', 
    'b_lat,nsrdb_lat', 
    'b_lon,nsrdb_lon', 
]
for pair in rename_cols:
    run('v.db.renamecolumn', map='solar_cluster_land_class_nsrdb_grid', column=pair)
run('v.db.dropcolumn', map='solar_cluster_land_class_nsrdb_grid', columns='a_cat')

#####
# create rasters showing latitude, longitude and coverage for all allowed solar areas

bands = ['00_to_15', '15_to_30']
for b in bands:
    # note: x() and y() will be in degrees if dem_crs uses geographic coordinates
    run('r.mapcalc', expression=f'allowed_{b}_lat = allowed_{b} * y()')
    run('r.mapcalc', expression=f'allowed_{b}_lon = allowed_{b} * x()')
# coverage shows 1 for allowed and 0 for non-allowed cells for whole region
for b in bands:
    run('r.mapcalc', expression=f'allowed_{b}_covg = if(isnull(allowed_{b}), 0, 1)')

#####
# Calculate solar area in each polygon of solar_cluster_land_class_nsrdb_grid
# (will be in square meters)
# We store data for different slope classes in different columns instead of
# different rows, which would be more normal. This is necessary because there's
# only one polygon for each region. Later, we split the data for each slope class 
# into different rows in a database table, just before it is exported. 
for b in bands:
    for stat in ['lat', 'lon', 'covg']:
        run(
            'v.rast.stats',
            map='solar_cluster_land_class_nsrdb_grid',
            raster=f'allowed_{b}_{stat}',
            method='average',
            column_prefix=f'{stat}_{b}'
        )
        run(
            'v.db.renamecolumn', 
            map='solar_cluster_land_class_nsrdb_grid', 
            column=f'{stat}_{b}_average,solar_{b}_{stat}'
        )
    # calculate solar area as region area * percent coverage
    run(  # dummy field with full area
        'v.to.db', 
        map='solar_cluster_land_class_nsrdb_grid', 
        option='area', 
        columns=f'solar_{b}_area',
    )
    run( # convert dummy field to area available for solar
        'v.db.update', 
        map='solar_cluster_land_class_nsrdb_grid', 
        column=f'solar_{b}_area',
        query_column=f'solar_{b}_area * solar_{b}_covg',
    )

# Note:
# For some reason solar_cluster_land_class_nsrdb_grid has some area features
# with no centroid or cat or attributes. This gives a warning when 
# exporting, but seems OK to ignore. I can't find any way to identify
# these features (not sure they even have any spatial extent).

# If you run next two lines (our main code path), then it shows 2784 areas and 
# 2183 centroids, and reports a warning that features without categories were skipped.
# !v.info map=solar_cluster_land_class_nsrdb_grid
# !v.out.ogr --o input=solar_cluster_land_class_nsrdb_grid output=solar_cluster_land_class_nsrdb_grid.csv format=CSV type=area
# (using type=centroid gives similar results but more warnings)

# If you run the next three lines, then the export still reports features
# without categories were skipped.
# !v.category --o input=solar_cluster_land_class_nsrdb_grid layer=-1 option=add cat=1000 output=scng_fixed type=area
# !v.info map=scng_fixed
# !v.out.ogr --o input=scng_fixed output=scng_fixed.csv format=CSV type=area

# If you run the next three lines (adding centroids), then it shows centroids
# for each area and you no longer get a warning about a feature with a missing
# category, but you get an informational message that features without
# attributes were written and the .csv output file ends up with rows with only a
# cat (presumably from the centroid) and no other attributes. 
# !v.centroids --o input=solar_cluster_land_class_nsrdb_grid option=add cat=1000 output=scng_fixed
# !v.info map=scng_fixed
# !v.out.ogr --o input=scng_fixed output=scng_fixed.csv format=CSV type=area

# If you then try to add more attributes to help find this feature, e.g., with
# the two lines below, then those columns get added to the output for most rows,
# but some rows still get exported with only a cat and no other attributes.
# run('v.to.db', map='scng_fixed', option='coor', columns='lon,lat')
# run('v.to.db', map='scng_fixed', option='area', columns='area')

######
# Make solar_project table and solar_project.csv with the final output. 
# This splits slope bands into different rows and aggregates all areas
# with the same load_zone, land_class and nsrdb_grid.
run('db.execute', sql='DROP TABLE IF EXISTS solar_project')
run(
    'db.execute',
    sql=f"""
        CREATE TABLE solar_project (
            siteid INTEGER PRIMARY KEY AUTOINCREMENT,
            cluster_id INT,
            load_zone VARCHAR,
            land_class VARCHAR,
            slope_class VARCHAR,
            nsrdb_id INT,
            nsrdb_lat REAL,
            nsrdb_lon REAL,
            solar_lat REAL,
            solar_lon REAL,
            solar_area REAL
        )
    """
)
for b in bands:
    run(
        'db.execute', 
        sql=f"""
            INSERT INTO solar_project (
                cluster_id, load_zone, land_class, slope_class, 
                nsrdb_id, nsrdb_lat, nsrdb_lon,
                solar_lat, solar_lon, solar_area
            )
            SELECT cluster_id, load_zone, land_class, '{b}' as slope_class,
                nsrdb_id, nsrdb_lat, nsrdb_lon,
                sum(solar_{b}_lat*solar_{b}_area)/sum(solar_{b}_area) AS solar_lat, 
                sum(solar_{b}_lon*solar_{b}_area)/sum(solar_{b}_area) AS solar_lon, 
                sum(solar_{b}_area) AS solar_area
            FROM solar_cluster_land_class_nsrdb_grid
            WHERE load_zone = 'Oahu'
            GROUP BY 1, 2, 3, 4, 5, 6, 7
            HAVING solar_area > 0
    """
)
run(
    'db.select',
    sql="SELECT * FROM solar_project",
    output='solar_project.csv',
    separator=','
)
print("============================================")
print("Script complete. Outputs are in solar_project.csv and load_zone_ag_land_area.csv")
print(f"You can delete the {config['grass_data_dir']} directory if not needed for map-making.")

# Output finished. solar_cluster_land_class_nsrdb_grid.csv contains total area
# available for solar in each combination of cluster, land class and nsrdb grid 
# cell. There are also separate columns for low-slope and high-slope land. These will later be aggregated 
# by cluster, land class and slope class, after calculating solar production for
# each grid cell. 

# Note: prior to May 2021, this was not split by land class or slope class. For
# a little while in May 2021, we did external reporting of the amount of land
# in each class that was in each project (single cluster, all land classes). 
# However, we stopped doing that because it became difficult to maintain with
# multiple slope classes as well, and because it was not completely sound, 
# since it implicitly assumed all land classes in the project had the same ratio 
# in each grid cell that made up the project. That may not be true, e.g., worse
# land classes could be in a cell that is steeper with more clouds. So now we
# just define a project for each of these combos, and total land use for each 
# project.

# TODO: update solar import script to 
# - read load_zone_ag_land_area.csv
# - use solar_project instead of solar_cluster_nsrdb_grid_final.csv
# - create separate projects for each cluster, slope class and land class.
# - save this information in the database and use higher costs for the steeper land (probably a completely different technology 
# category so get_scenario_data can exclude steeper projects if needed). Update
# scenario_data to export land_class data for each project. Add switch_hawaii
# code to enforce limits on each land class. Rerun scenarios with different limits
# on each class, and now just do intermediate limits, rather than trying to 
# readjust unlimited-use scenarios to minimize land use in restricted categories,
# without changing costs, since that isn't really possible now. (Or maybe try to 
# minimize use of restricted land while keeping costs fairly close? Probably better
# just to try different limits and report what happens.)
