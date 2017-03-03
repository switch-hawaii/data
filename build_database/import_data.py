#!/usr/bin/env python

"""construct postgresql backend database for SWITCH-Hawaii.
Data is pulled into this database from various sources (Excel files,
GIS results, NSRDB/OWITS data files), and then switch_mod.hawaii.scenario_data
can be used to construct any scenario using the accumulated data.
"""

import sys, csv, datetime, os, collections
from textwrap import dedent
import numpy as np
import pandas as pd
import sqlalchemy
import shared_tables, tracking_pv, util
from util import execute, executemany, pg_host, switch_db

try:
    import openpyxl
except ImportError:
    print "This script requires the openpyxl module to access the data in Microsoft Excel files."
    print "Please execute 'sudo pip install openpyxl' or 'pip install openpyxl' (Windows)."
    raise

db_engine = sqlalchemy.create_engine('postgresql://' + pg_host + '/' + switch_db)

# data files used by multiple functions
data_directory = '..'
psip_data_file = os.path.join(data_directory, 'Generator Info', 'PSIP 2016-12 generator data.xlsx')
# load scenario ID corresponding to this file
load_scen_id = 'PSIP_2016_12'

# technologies that are added to the generator_info and project tables by
# onshore_wind(), offshore_wind(), tracking_pv.tracking_pv() and tracking_pv.distributed_pv(),
# not by generator_info
renewable_techs = ['DistPV', 'CentralTrackingPV', 'CentralFixedPV', 'OnshoreWind', 'OffshoreWind']

def main():
    # run all the import scripts (or at least the ones we want)
    # ev_adoption()
    # fuel_costs()
    # # fuel_costs_no_biofuel()  # obsolete
    # energy_source_properties()
    # rps_timeseries()
    # system_load()
    # generator_info()
    # tracking_pv.tracking_pv()
    # tracking_pv.distributed_pv()
    # onshore_wind()
    # offshore_wind()
    shared_tables.calculate_interconnect_costs()

    pass

def data_dir(*path):
    return os.path.join(data_directory, *path)

def get_workbook(xlsx_file):
    return openpyxl.load_workbook(xlsx_file, data_only=True, read_only=True)

def get_table_from_xlsx(xlsx_file, named_range, transpose=False):
    wb = get_workbook(xlsx_file)
    full_range = wb.get_named_range(named_range)
    # note: named range should be a simple rectangular region;
    # if it contains more than one region we ignore all but the first
    d1 = full_range.destinations[0]
    ws = d1[0]
    region = d1[1]
    data = list(tuple(c.value for c in r) for r in ws[region])
    if transpose:
        data = zip(*data)
    head = data.pop(0)  # take out the header row
    data = zip(*data)   # switch from row to column orientation
    # make a dictionary, with one column for each element of header row
    return dict(zip(head, data))

def get_named_region(xlsx_file, named_range):
    # get a single rectangular region from the specified file in the source data directory
    wb = get_workbook(xlsx_file)
    full_range = wb.get_named_range(named_range)
    if full_range is None:
        raise ValueError(
            'Range "{}" not found in workbook "{}".'.format(named_range, xlsx_file))
    if len(full_range.destinations) > 1:
        raise ValueError(
            'Range "{}" in workbook "{}" contains more than one region.'.format(named_range, xlsx_file))
    ws, region = full_range.destinations[0]
    return ws[region]

def data_frame_from_xlsx(xlsx_file, named_range):
    region = get_named_region(xlsx_file, named_range)
    return pd.DataFrame([cell.value for cell in row] for row in region)

def get_named_cell_from_xlsx(xlsx_file, named_range):
    region = get_named_region(xlsx_file, named_range)
    if isinstance(region, collections.Iterable):
        raise ValueError(
            'Range "{}" in workbook "{}" does not refer to an individual cell.'.format(
                named_range, xlsx_file))
    return region.value

#########################
# rps timeseries (reusing version from before the server crashed)
def rps_timeseries():
    execute("""
        delete from study_date where time_sample='rps';
        delete from study_hour where time_sample='rps';
    """)

    with open('timeseries_rps.tab','r') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            dt = str(r["TIMESERIES"])[2:8]
            execute("""
                INSERT INTO study_date
                    (time_sample, period, study_date,
                    month_of_year, date,
                    hours_in_sample,
                    ts_num_tps, ts_duration_of_tp, ts_scale_to_period)
                VALUES
                    (%s, %s, %s,
                    %s, %s,
                    %s,
                    %s, %s, %s);
            """,
                ('rps', r["ts_period"], r["TIMESERIES"],
                int(dt[2:4]), "20"+dt[0:2]+"-"+dt[2:4]+"-"+dt[4:6],
                float(r["ts_duration_of_tp"])*float(r["ts_scale_to_period"]),
                r["ts_num_tps"], r["ts_duration_of_tp"], r["ts_scale_to_period"])
            )

    with open('timepoints_rps.tab','r') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            i += 1
            sys.stdout.write('row: {}\r'.format(i))
            sys.stdout.flush()
            t = r["timestamp"][5:]
            dt = str(r["timeseries"])[2:8]
            execute("""
                INSERT INTO study_hour
                    (time_sample, study_date, study_hour,
                    hour_of_day, date_time)
                VALUES
                    (%s, %s, %s,
                    %s,
                    cast(%s as timestamp with time zone));
            """,
                ('rps', r["timeseries"], r["timepoint_id"],
                int(t[6:8]),
                "20" + dt[0:2] + '-' + t[:5] + ' ' + t[6:] + ':00-10')
            )

    execute("""
        delete from study_periods where time_sample='rps';
        insert into study_periods
            select distinct time_sample, period from study_date
            where time_sample='rps'
            order by 2;

        -- mini RPS study (even months, even hours) for reasonable results in relatively short time

        drop table if exists tdoublemonths;
        create temporary table tdoublemonths
            (month_of_year smallint primary key, days_in_month smallint);
        insert into tdoublemonths values
          (1, 59), (2, 59), (3, 61), (4, 61), (5, 61), (6, 61),
          (7, 62), (8, 62), (9, 61), (10, 61), (11, 61), (12, 61);

        delete from study_date where time_sample='rps_mini';
        insert into study_date
            (period, study_date, month_of_year, date,
            hours_in_sample, time_sample, ts_num_tps, ts_duration_of_tp, ts_scale_to_period)
            select period, study_date, d.month_of_year, date,
                0.0 as hours_in_sample,
                'rps_mini' as time_sample,
                12 as ts_num_tps, 2.0 as ts_duration_of_tp,
                case when ts_scale_to_period < 100 then 2*%(years_per_period)s
                    else (days_in_month-2)*%(years_per_period)s end as ts_scale_to_period
            from study_date d join tdoublemonths m using (month_of_year)
            where time_sample = 'rps' and mod(month_of_year, 2)=0
            order by 1, 3, 5 desc, 4;

        delete from study_hour where time_sample='rps_mini';
        insert into study_hour (study_date, study_hour, hour_of_day, date_time, time_sample)
          select h.study_date, study_hour, hour_of_day, date_time,
            'rps_mini' as time_sample
          from study_hour h join study_date d on (d.time_sample='rps_mini' and d.study_date=h.study_date)
          where h.time_sample = 'rps' and mod(hour_of_day, 2)=0
          order by period, month_of_year, hours_in_sample desc, hour_of_day;

        delete from study_periods where time_sample='rps_mini';
        insert into study_periods
            select distinct time_sample, period from study_date
            where time_sample='rps_mini'
            order by 2;


    """, dict(years_per_period=8))


    print "Created rps and rps_mini time samples."


#########################
# ev adoption
def ev_adoption():
    # identify pairs of (ev_scen_id, HECO scenario name):
    ev_adoption_scenarios=(
        'Business as Usual', # straight line to 4.3% by 2045
        'No Burning Desire', # 2013 IRP, 17.5% by 2045
        'Stuck in the Middle', # 2013 IRP, a.k.a. 'Moved by Passion', 35.2% by 2045
        'Blazing a Bold Frontier', # 2013 IRP, 70.1% by 2045
        'PSIP 2016-12', # about 55% by 2045
        'Full Adoption', # 100% by 2045
        'Half Adoption', # 50% by 2045
        'Flat 2016', # 0.5% all the way through
    )
    # get the EV adoption curves from an Excel workbook
    # uses logistic curves fitted to HECO IRP 2013 Appendix E-10, p. E-113,
    # as well as VMT data from DBEDT Economic Databook
    # and vehicle registration rates from DBEDT monthly energy spreadsheet
    ev_adoption_curves = get_table_from_xlsx(
        data_dir("EV Adoption", "EV projections.xlsx"),
        named_range='ev_data'
    )

    # create the ev_adoption table
    execute("""
        DROP TABLE IF EXISTS ev_adoption;
        CREATE TABLE ev_adoption (
            load_zone varchar(40),
            year int,
            ev_scenario varchar(40),
            ev_share float,
            ice_miles_per_gallon float,
            ev_miles_per_kwh float,
            ev_extra_cost_per_vehicle_year float,
            n_all_vehicles float,
            vmt_per_vehicle float
        );
    """)

    # insert data into the ev_adoption table
    n_rows = len(ev_adoption_curves['Year'])
    for ev_scenario in ev_adoption_scenarios:
        executemany(
            "INSERT INTO ev_adoption VALUES ({})".format(','.join(["%s"]*9)),
            zip(
                ['Oahu']*n_rows, ev_adoption_curves['Year'], [ev_scenario]*n_rows,
                ev_adoption_curves[ev_scenario], # % adoption
                ev_adoption_curves["ICE miles per gallon"],
                ev_adoption_curves["EV miles per kWh"],
                ev_adoption_curves["EV extra cost per vehicle per year"],
                ev_adoption_curves["number of vehicles"],
                ev_adoption_curves["VMT per vehicle"],
            )
        )

    # add Ulupono data series
    # NOTE: for this series, we are only interested in the EVs, so we model them as if they
    # were the whole fleet. That way, when we report total costs, they don't include any ICE costs.
    uev_scenario = 'Coffman - Ulupono'

    uev = data_frame_from_xlsx(
        data_dir('Ulupono', 'Project EV Electricity Demand_Coffman Reference.xlsx'),
        'ev_data'
    ).T.set_index(0).T
    uev['load_zone'] = 'Oahu'
    uev=uev.rename(columns={"Year": "year"}).set_index(['load_zone', 'year'])

    uev_final = pd.DataFrame(dict(
    ev_scenario=uev_scenario,
    ev_share=1.0,
    ice_miles_per_gallon=30,    # arbitary value, not used
    ev_miles_per_kwh=
        (uev["# of EV's on Road"].values * uev["VMT per Vehicle"].values).sum(axis=1)
        / uev["Electricity (GWh)"].sum(axis=1) / 1e6,
    ev_extra_cost_per_vehicle_year=0.0,
    n_all_vehicles=uev["# of EV's on Road"].sum(axis=1),
    vmt_per_vehicle=
        (uev["# of EV's on Road"].values * uev["VMT per Vehicle"].values).sum(axis=1)
        / uev["# of EV's on Road"].sum(axis=1)
    ))
    # verify it matches the spreadsheet:
    print (
        "The following values should match the energy consumption in " 
        + data_dir('Ulupono', 'Project EV Electricity Demand_Coffman Reference.xlsx')
    )
    print uev_final['n_all_vehicles']*uev_final['ev_share']*uev_final['vmt_per_vehicle']/uev_final['ev_miles_per_kwh']/1e6

    # drop existing records
    execute("""
        DELETE FROM ev_adoption WHERE ev_scenario=%s;
    """, (uev_scenario,))
    uev_final.to_sql('ev_adoption', db_engine, if_exists='append')

    # set n_all_vecicles = total number of EVs each year (cars and trucks)
    # set efficiency and vmt per vehicle based on this
    # set ice efficiency to some arbitrary number (e.g., 30 mpg)

    print "Created ev_adoption table."

    # create the ev_hourly_charge_profile table (simple business-as-usual charging profile,
    # given as hourly weights)
    # see /Users/matthias/Dropbox/Research/shared/Paritosh/M.S Thesis Paritosh/Data Used In Thesis/calculate BAU charging.ipynb
    execute("""
        DROP TABLE IF EXISTS ev_hourly_charge_profile;
        CREATE TABLE ev_hourly_charge_profile (
            hour_of_day smallint,
            charge_weight float
        );
    """)
    with open(data_dir('EV Adoption', 'ev_hourly_charge_profile.tsv')) as f:
        profile = [r.split("\t") for r in f.read().splitlines()][1:] # skip headers

    executemany(
        "INSERT INTO ev_hourly_charge_profile (hour_of_day, charge_weight) VALUES (%s, %s);",
        profile
    )
    print "Created ev_hourly_charge_profile table."



def fuel_costs():
    # create the fuel_costs table if needed
    execute("""
        CREATE TABLE IF NOT EXISTS fuel_costs (
            load_zone varchar(40),
            year int,
            base_year int,
            fuel_type varchar(30),
            price_mmbtu float,
            fixed_cost float,
            max_avail_at_cost float,
            fuel_scen_id varchar(40),
            tier varchar(20),
            max_age int
        );
        ALTER TABLE fuel_costs OWNER TO admin;
    """)

    # TODO: add fixed_cost and max_avail_at_cost for EIA-based forecasts

    def eia_dir(*path):
        return data_dir('EIA-based fuel cost forecasts', *path)

    # Oahu fuel price forecasts, derived from EIA
    # import_eia_fuel_costs(eia_dir("HECO fuel cost forecasts.xlsx"), 'EIA_ref')
    # import_eia_fuel_costs(eia_dir("HECO fuel cost forecasts_low.xlsx"), 'EIA_low')
    # import_eia_fuel_costs(eia_dir("HECO fuel cost forecasts_high.xlsx"), 'EIA_high')
    # import_eia_fuel_costs(eia_dir("HECO fuel cost forecasts_LNG_pegged_to_oil.xlsx"), 'EIA_lng_oil_peg')
    # import_eia_fuel_costs(eia_dir("HECO fuel cost forecasts_high_LNG_pegged_to_oil.xlsx"), 'EIA_high_lng_oil_peg')

    # Oahu hedged fuel costs and equivalent unheged costs from HECO 
    # (note: we use these instead of the PSIP Fuel Price Forecasts workbook because 
    # these adjust to 2016 dollars and include LNG with various durations)
    # import_hedged_fuel_costs(eia_dir("hedged fuel prices.xlsx"), tag='hedged')

    hedged_fuel_scen_id = 'hedged_2016_11_22'
    standard_fuel_scen_id = 'unhedged_2016_11_22'
    import_hedged_fuel_costs(eia_dir("hedged fuel prices 2016-11-22.xlsx"), tag=hedged_fuel_scen_id)
    import_hedged_fuel_costs(eia_dir("unhedged fuel prices 2016-11-22.xlsx"), tag=standard_fuel_scen_id)

    # import_psip_fuel_costs(data_dir("HECO Plans/PSIP-WebDAV/Resource Assumptions/PSIP Fuel Price Forecasts for HE 2016-06-27 regressions.xlsx"))

    # flat fuel price based on 2017 prices in 'unhedged_2016_11_22' 
    execute("""
        CREATE TEMPORARY TABLE tfuelcosts AS
            SELECT * FROM fuel_costs WHERE fuel_scen_id=%s;
        UPDATE TFUELCOSTS a
            SET fuel_scen_id='flat_2016', price_mmbtu=b.price_mmbtu
            FROM tfuelcosts b
            WHERE b.year=2016 AND b.load_zone=a.load_zone AND b.fuel_type=a.fuel_type AND b.tier=a.tier;
        INSERT INTO fuel_costs SELECT * FROM tfuelcosts;
        DROP TABLE tfuelcosts;
    """, (standard_fuel_scen_id,))

def import_eia_fuel_costs(file, fuel_scen_id):

    # get the forecasts from an Excel workbook
    # Based on various sources, cited in the workbook, extended to 2050
    fuel_forecast = get_table_from_xlsx(file, named_range='Adjusted_EIA_Forecast', transpose=True)

    # note: all the EIA spreadsheets use a base year of 2013
    base_year = 2013

    # remove any existing records
    execute("""
        DELETE FROM fuel_costs WHERE fuel_scen_id=%s;
    """, (fuel_scen_id,))

    # take out the list of years, so the dictionary just has one entry for each fuel
    years = fuel_forecast.pop('Year')

    # insert data into the fuel_costs table
    n_rows = len(years)
    for f in fuel_forecast:
        ft=f.split(", ")
        fuel = ft[0]
        tier = ft[1] if len(ft) >= 2 else 'base'
        executemany("""
            INSERT INTO fuel_costs (load_zone, year, base_year, fuel_type, price_mmbtu, fuel_scen_id, tier)
            VALUES (%s, %s, %s, %s, %s, %s, %s)""",
            zip(['Oahu']*n_rows,
                years,
                [base_year]*n_rows,
                [fuel]*n_rows,
                fuel_forecast[f],
                [fuel_scen_id]*n_rows,
                [tier]*n_rows
            )
        )

    print "Added EIA-derived forecast (fuel_scen_id={}) to fuel_costs table.".format(fuel_scen_id)


def import_hedged_fuel_costs(file, tag='hedged'):

    prices = data_frame_from_xlsx(file, named_range='fuel_prices')
    prices = prices.set_index(0)
    prices.index.name = 'year'
    prices = prices.T.set_index(['fuel_type', 'tier']).T.astype(float)
    # switch to one row per value, and assign a name to the value
    prices = pd.DataFrame({'price_mmbtu': prices.stack(['fuel_type', 'tier'])})
    prices['load_zone'] = 'Oahu'
    prices['base_year'] = get_named_cell_from_xlsx(file, named_range='base_year')

    tiers = data_frame_from_xlsx(file, named_range='tier_properties')
    # Transpose, set row and column labels, and convert to floating point (converting None to NaN)
    tiers = tiers.set_index(0).T.set_index(['fuel_type', 'tier']).astype(float)

    # fixed prices vary depending on the finance term; terms are pre-specified in this region
    fixed_costs = data_frame_from_xlsx(file, named_range='tier_fixed_costs')
    # use the first column as indexes (mostly to get column names), then set column headers
    fixed_costs = fixed_costs.set_index(0).T.set_index(['fuel_type', 'tier']).T
    # drop unneeded row for current finance term (we only want the values from the data table below that)
    fixed_costs = fixed_costs.iloc[1:]
    # give the index a name
    fixed_costs.index.name = 'term'
    # convert to row-wise format, give the fixed_cost column a name, and convert the indexes to columns
    fixed_costs = pd.DataFrame({'fixed_cost': fixed_costs.unstack()}).reset_index()
    # add a fuel_scen_id
    fixed_costs['fuel_scen_id'] = tag
    # use the term column as the maximum age for each tier with non-zero fixed costs
    limited_life = fixed_costs['fixed_cost'] > 0
    fixed_costs.loc[limited_life, 'max_age'] = fixed_costs.loc[limited_life, 'term']
    del fixed_costs['term']

    # remove duplicate rows (we don't need multiple rows with multiple ages for the $0 cost tiers)
    # also restore the indexes, to enable joining later
    fixed_costs = fixed_costs.drop_duplicates().set_index(['fuel_type', 'tier'])
    # merge the columns into the tiers table (adding all fuel_scen_id's and max_age's)
    tiers = tiers.join(fixed_costs)

    # merge the columns into the prices table (have to drop the year index to make this work)
    prices = prices.reset_index('year').join(tiers)

    # add the project lifespan into the tier id (have to convert tier index to a column to do this,
    # so might as well remove all indexes)
    prices = prices.reset_index()
    limited_life = prices['fixed_cost'] > 0
    prices.loc[limited_life, 'tier'] += '_' + prices.loc[limited_life, 'max_age'].astype(int).astype(str)

    # restore the indexes and sort the table
    prices = prices.set_index(['fuel_scen_id', 'year', 'fuel_type', 'tier']).sort_index()

    # remove any existing records
    execute("DELETE FROM fuel_costs WHERE fuel_scen_id LIKE %s;", (tag,))

    prices.to_sql('fuel_costs', db_engine, if_exists='append')

    print "Added hedged prices (fuel_scen_id = {}) to fuel_costs table.".format(list(prices.index.levels[0]))


def import_psip_fuel_costs(file):
    # TODO: change this to do a more complete treatment of LNG options and coal
    # (not important immediately because we're just using this for "greenfield"
    # analysis of demand response)

    file = data_dir+"/HECO Plans/PSIP-WebDAV/Resource Assumptions/PSIP Fuel Price Forecasts for HE 2016-06-27 regressions.xlsx"
    fuel_scen_id = 'PSIP_2016_09'

    prices = data_frame_from_xlsx(file, named_range='real_fuel_prices').T.set_index(0).T
    year = data_frame_from_xlsx(file, named_range='years')

    prices = prices.set_index(year[0])
    prices = prices.astype(float)
    # drop unneeded columns and rename the remaining natural gas column
    del prices['ULSD']
    del prices['HECO LNG commodity']
    del prices['HECO LNG delivered']
    prices.rename(columns={'HG LNG delivered': 'LNG'}, inplace=True)

    # switch to one row per value, and assign a name to the value
    prices = pd.DataFrame({'price_mmbtu': prices.stack()})
    prices.index.rename(['year', 'fuel_type'], inplace=True)

    prices['load_zone'] = 'Oahu'
    prices['fuel_scen_id'] = fuel_scen_id
    prices['tier'] = 'base'
    prices['fixed_cost'] = 0
    prices['base_year'] = get_named_cell_from_xlsx(file, named_range='base_year')

    # remove any existing records
    execute("DELETE FROM fuel_costs WHERE fuel_scen_id like %s;", (fuel_scen_id,))

    prices.to_sql('fuel_costs', db_engine, if_exists='append')

    # reuse existing solid fuel data
    execute("""
        INSERT INTO fuel_costs
            SELECT load_zone, year, fuel_type, price_mmbtu,
                %s as fuel_scen_id, tier, fixed_cost, max_avail_at_cost, base_year
            FROM fuel_costs
            WHERE fuel_scen_id = 'EIA_ref' AND fuel_type in ('Coal', 'Pellet-Biomass');
    """, (fuel_scen_id,))

    # convert LNG to a 'bulk' tier and lookup relevant data
    execute("""
        UPDATE fuel_costs AS a
          SET tier = b.tier, fixed_cost = b.fixed_cost, max_avail_at_cost = b.max_avail_at_cost
          FROM fuel_costs b
          WHERE a.fuel_scen_id = %s AND a.fuel_type = 'LNG' AND a.tier = 'base'
              AND a.fuel_type = b.fuel_type
              AND b.fuel_scen_id = 'hedged_20' AND b.tier = 'bulk';
    """, (fuel_scen_id,))

    # # ULSD is not in the energy_source_properties database and isn't used in current scenarios
    # execute("DELETE FROM fuel_costs WHERE fuel_type = 'ULSD' AND fuel_scen_id = %s;", (fuel_scen_id,))

    print "Added PSIP prices (fuel_scen_id = {}) to fuel_costs table.".format(fuel_scen_id)



#########################
# Fuel properties, maintained manually in the Excel forecast workbook
def energy_source_properties():
    properties = get_table_from_xlsx(
        data_dir+"/EIA-based fuel cost forecasts/HECO fuel cost forecasts.xlsx",
        named_range='Fuel_Properties'
    )

    # create the fuel_properties table if needed
    execute("""
        CREATE TABLE IF NOT EXISTS energy_source_properties (
            energy_source VARCHAR(30) PRIMARY KEY,      -- name of the fuel
            fuel_rank DECIMAL(4, 2),           -- usually 1-5, but may be decimal, e.g., 1.5
            rps_eligible SMALLINT,             -- 0 or 1
            co2_intensity FLOAT                -- tCO2 per MMBtu
        );
    """)

    # create a temporary table to hold the data before aggregating by fuel type
    execute("""
        DROP TABLE IF EXISTS t_energy_source_properties;
        CREATE TEMPORARY TABLE t_energy_source_properties (LIKE energy_source_properties);
    """)

    # insert data into the energy_source_properties table
    executemany("""
        INSERT INTO t_energy_source_properties (energy_source, fuel_rank, rps_eligible, co2_intensity)
        VALUES (%s, %s, %s, %s)""",
        zip(
            [f.split(', ')[0] for f in properties['Fuel']],
            properties['Rank'],
            properties['RPS Eligible'],
            [i/1000.0 for i in properties['kg CO2 per MMbtu']],
        )
    )

    # move the data into the main energy_source_properties table
    execute("""
        DELETE FROM energy_source_properties;
        INSERT INTO energy_source_properties SELECT DISTINCT * FROM t_energy_source_properties;
        DROP TABLE t_energy_source_properties;
    """)

    print "Created energy_source_properties table."

def fuel_costs_no_biofuel():
    """Create no-biofuel fuel cost scenarios"""
    # note: these are not used anymore; the same effect can be achieved by setting
    # '--biofuel-limit 0'
    execute("""
        DELETE FROM fuel_costs WHERE fuel_scen_id LIKE 'EIA_%_no_biofuel';
        DROP TABLE IF EXISTS t_fuel_costs_no_biofuel;
        CREATE TABLE t_fuel_costs_no_biofuel AS
        SELECT c.*
            FROM fuel_costs c JOIN energy_source_properties p ON c.fuel_type = p.energy_source
            WHERE rps_eligible = 0 AND fuel_scen_id LIKE 'EIA_%';
        UPDATE t_fuel_costs_no_biofuel SET fuel_scen_id = fuel_scen_id || '_no_biofuel';
        INSERT INTO fuel_costs SELECT * FROM t_fuel_costs_no_biofuel;
        DROP TABLE t_fuel_costs_no_biofuel;
    """)

def onshore_wind():
    """Import old onshore wind data into newer tables."""
    # TODO: write code to create these records directly from OWITS data and GIS files
    # and also store location and interconnect distance
    execute("""
        delete from cap_factor
            where project_id in
                (select project_id from project where technology = 'OnshoreWind');
        delete from project where technology='OnshoreWind';
        insert into project
            (load_zone, technology, site, orientation, max_capacity, connect_distance_km)
            select load_zone, technology, concat('OnWind_', site),
                orientation, max_capacity, connect_length_km
            from max_capacity_pre_2016_06_21
                left join connect_cost_pre_2016_06_21 using (load_zone, technology, site, orientation)
                where technology='OnshoreWind';
        insert into cap_factor (project_id, date_time, cap_factor)
            select project_id, date_time, cap_factor
                from cap_factor_pre_2016_06_21 c
                    join project p on
                        (p.load_zone=c.load_zone and p.technology=c.technology and
                        p.site = concat('OnWind_', c.site) and p.orientation=c.orientation)
                    where c.technology='OnshoreWind';
    """)
    # also collect data for existing wind farms, stored in "existing_plants" tables
    execute("""
        insert into project
            (load_zone, technology, site, orientation, max_capacity, connect_distance_km)
            select load_zone, 'OnshoreWind' AS technology, CONCAT('OnWind_', plant_name) AS site,
                'na' AS orientation, peak_mw AS max_capacity, 0 AS connect_length_km
            from existing_plants
            where aer_fuel_code='WND';
        -- set latitude and longitude for existing wind farms so they can be matched during
        -- import of existing project profiles in new_generator_info() later
        update project set latitude=21.682, longitude=-157.976 where site='OnWind_Kahuku';
        update project set latitude=21.620, longitude=-158.048 where site='OnWind_Kawailoa';
        insert into cap_factor (project_id, date_time, cap_factor)
            select p.project_id, date_time, cap_factor
            from existing_plants e join existing_plants_cap_factor c using (project_id)
                join project p on (
                    p.load_zone=c.load_zone and
                    p.technology='OnshoreWind' and c.technology='Wind'
                    and p.site=CONCAT('OnWind_', plant_name))
            order by 1, 2;
    """)

def offshore_wind():
    """Import capacity factor for offshore wind farms. This is calculated as the
    average of three proposed offshore sites to get approximately the right amount
    for diversified offshore wind. (It might be better just to model them as three
    separate projects.)

    Note: The 2016 PSIP used hourly output possibly for 2014, from an existing
    wind farm on the Big Island with a capacity factor of 42%. We don't use this
    because it's the wrong profile for offshore Oahu, and especially because it
    has inconsistent timing with our other weather and load data so it would
    create an artificial appearance of diversity (strong winds when Oahu actually
    has windless/sunless days).
    """
    # approximate locations for the centers of three proposed wind farms
    # were found on 2016-04-07 by inspecting the
    # "Atlantic and Pacific OCS Wind Planning Areas and Wind Energy Areas"
    # shapefile from http://www.boem.gov/Renewable-Energy-GIS-Data/
    # (http://www.boem.gov/uploadedFiles/BOEM/Renewable_Energy_Program/Mapping_and_Data/Wind_Planning_Areas.zip)
    locs = np.array([[21.656, -158.572], [21.096, -157.987], [20.969, -157.799]])
    cells = np.array(list(execute("select i, j, lat, lon from cell where grid_id = 'E';")))
    cell_lat_lon = cells[:,-2:]
    # this makes one row for each site, one col for each cell, showing distance in degrees**2
    dist2 = ((locs[:,np.newaxis,:] - cell_lat_lon[np.newaxis,:,:])**2).sum(axis=2)
    match_cells = dist2.argmin(axis=1)
    turbine_cells = cells[match_cells]

    # normalized power curve for generic offshore wind turbine from http://www.nrel.gov/docs/fy14osti/61714.pdf p. 5,
    # with operating range extended to 30 m/s like Repower 6 M shown on p. 4.
    power_curve = np.array(zip(
        range(32),
        [0] * 4 + [0.0281, 0.074, 0.1373, 0.2266, 0.3443, 0.4908, 0.6623, 0.815, 0.9179, 0.9798]
        + [1] * 17 + [0]
    ))

    hourly = pd.DataFrame(
        list(execute(
            """
                select grid_id, i, j, complete_time_stamp as date_time, s100
                    from hourly_average
                    where grid_id='E' and (i, j) in %(cells)s;""",
            {"cells": tuple(tuple(c for c in r) for r in turbine_cells[:,:2].astype(int))}
        )),
        columns=['grid_id', 'i', 'j', 'date_time', 's100'],
    )
    hourly.set_index(['date_time', 'grid_id', 'i', 'j'], inplace=True)

    hourly['p100'] = np.interp(hourly['s100'].values, power_curve[:,0], power_curve[:,1])
    hourly = hourly.unstack(level=['grid_id', 'i', 'j'])
    # use mean of three sites as hourly output; also derate same as we do for land sites (from IRP 2013)
    hourly['power'] = hourly['p100'].mean(axis=1)*0.8747

    # delete any old OffshoreWind records from cap_factor
    execute("""
        delete from cap_factor
            where project_id in
                (select project_id from project where load_zone = 'Oahu' and technology = 'OffshoreWind');
    """)

    # add the new project to the project table
    execute("""
        delete from project where technology = 'OffshoreWind' and load_zone = 'Oahu';
        insert into project
            (load_zone, technology, site, orientation, max_capacity)
            values ('Oahu', 'OffshoreWind', 'OffWind', 'na', 800);
    """)
    # retrieve the project_id for the new project
    project_id = execute("""
        select project_id from project where load_zone = 'Oahu' and technology = 'OffshoreWind';
    """).next()[0]

    # put the power data into cap_factor
    out_data = zip([project_id]*len(hourly), hourly.index.astype(datetime.datetime), hourly['power'].values)

    # TODO: consider removing and restoring indexes before this
    executemany("""
        insert into cap_factor (project_id, date_time, cap_factor)
        values (%s, %s, %s);
    """, out_data)

    # note: we don't add latitude, longitude or interconnect_id (and cost) because we don't
    # have project-specific connection costs for them. So they will automatically use
    # the generic connection cost from generator_info (assigned later).
    # That happens to be zero in this case since the connection cost is included in the overnight cost.

def generator_info():
    # note: for now, these must always be run in this sequence, because
    # new_generator_info() creates the generator_info and part_load_fuel_consumption
    # tables, and then existing_generator_info() appends to them (without clearing
    # out any prior records)
    shared_tables.create_table('project')
    shared_tables.create_table('generator_info')
    new_generator_info()
    existing_generator_info()

def new_generator_info():
    """Read data from 'PSIP 2016-09 generator data.xlsx and
    store it in generator_info and generator_costs_by_year.

    This spreadsheet is based on HECO's 2016-09-07 PSIP assumptions, shown in their workplan
    and stored in a separate spreadsheet on 2016-06-16."""

    base_year = get_named_cell_from_xlsx(psip_data_file, 'o_m_base_year')

    gen_info = data_frame_from_xlsx(psip_data_file, 'technology_info')
    # set column headers and row indexes (index in the dataframe become index in the table)
    gen_info = gen_info.T.set_index(0).T.set_index('technology')
    gen_info.rename(
        columns={
            'fixed_o_m_per_kw_year': 'fixed_o_m',
            'variable_o_m_per_mwh': 'variable_o_m',
            'full_load_heat_rate': 'heat_rate'
        }, inplace=True)
    # convert from cost per MWh to cost per kWh
    gen_info['variable_o_m'] *= 0.001
    # convert unit_size = NA or 0 to NaN
    gen_info['unit_size'] = gen_info['unit_size'].where(
        (gen_info['unit_size'] != "NA") & (gen_info['unit_size'] != 0)
    )
    # report base_year for inflation calculations later
    gen_info['base_year'] = base_year

    # convert all columns except fuel to numeric values (some of these were imported
    # as objects due to invalid values, which have now been removed)
    # gen_info.convert_objects() does this nicely, but is deprecated.
    for c in gen_info.columns:
        if c != 'fuel':
            gen_info[c] = pd.to_numeric(gen_info[c])

    # remove all existing technology definitions (they will all be re-created
    # here and in existing_generator_info)
    execute('DELETE FROM generator_info;')
    # store data
    gen_info.to_sql('generator_info', db_engine, if_exists='append')

    # load gen capital cost info
    gen_costs = data_frame_from_xlsx(psip_data_file, 'technology_costs')
    inflation_rate = get_named_cell_from_xlsx(psip_data_file, 'inflation_rate')

    gen_costs = gen_costs.T.set_index(0).T
    gen_costs.columns.name = 'technology'
    # drop info rows (now 0-6)
    gen_costs = gen_costs.iloc[7:]
    # rename the first column
    gen_costs = gen_costs.rename(columns={gen_costs.columns[0]: 'year'})
    # set year index
    gen_costs = gen_costs.set_index('year')
    # convert N/A to nan
    gen_costs = gen_costs.where((gen_costs != ' N/A ') & (gen_costs != 'N/A')).astype(float)

    # select only the technologies that are in gen_info
    gen_costs = gen_costs[gen_info.index]

    # convert to real dollars in base year
    gen_costs = gen_costs.multiply(
        (1 + inflation_rate) ** (base_year - gen_costs.index.values),
        axis='rows'
    )

    # switch to stacked orientation (one row per year, tech), but keep as a DataFrame
    gen_costs = pd.DataFrame({'capital_cost_per_kw': gen_costs.stack()})
    # record the base year to allow adjustment to other years later
    gen_costs['base_year'] = base_year
    gen_costs['cap_cost_scen_id'] = 'psip_1609'
    gen_costs.to_sql('generator_costs_by_year', db_engine, if_exists='replace')

    # make a new version with 2017 costs all the way through
    execute("""
        INSERT INTO generator_costs_by_year 
            (year, technology, capital_cost_per_kw, base_year, cap_cost_scen_id)
            SELECT 
                a.year, a.technology, b.capital_cost_per_kw, b.base_year, 
                'psip_1609_flat' AS cap_cost_scen_id
            FROM generator_costs_by_year a 
                JOIN generator_costs_by_year b USING (cap_cost_scen_id, technology)
            WHERE a.cap_cost_scen_id = 'psip_1609' AND b.year = 2017;
    """)
    # or maybe:
    # gen_costs.loc[:, 'capital_cost_per_kw'] = gen_costs.loc[2017, 'capital_cost_per_kw']
    # gen_costs.to_sql('generator_costs_by_year', db_engine, if_exists='append')

    # import part-load heat rates
    gen_fuel_cons = data_frame_from_xlsx(psip_data_file, 'part_load_fuel_consumption')
    gen_fuel_cons = gen_fuel_cons.T.set_index(0).T
    gen_fuel_cons = gen_fuel_cons.rename(columns={
        'load level (MW)': 'output_mw',
        'fuel consumption (MMBtu/h)': 'fuel_consumption_mmbtu_per_h',
    }).set_index('technology')
    gen_fuel_cons.to_sql('part_load_fuel_consumption', db_engine, if_exists='replace')

    #############
    # import definitions for non-renewable/non-resource-limited projects
    # We could just construct these in scenario_data.py from the generator_info entries,
    # except that we need to specify a maximum capacity for each project to support
    # the RPS calculation (which uses that in a big-M constraint that allocates output
    # among fuels)
    project_info = data_frame_from_xlsx(psip_data_file, 'non_renewable_project_info') \
        .T.set_index(0) \
        .T.set_index(['load_zone', 'technology', 'site', 'orientation'])
    # remove any non-renewable project definitions from project table
    execute(
        'DELETE FROM project WHERE technology NOT IN %s',
        [tuple(renewable_techs)]
    )
    # insert non-renewable project definitions into project table
    project_info.to_sql('project', db_engine, if_exists='append')


def existing_generator_info():
    """copy data from 'Data/Generator Info/Existing Plant Data.xlsx' into
    generator_info, part_load_fuel_consumption, project and proj_existing_builds
    """
    gen_info_file = data_dir('Generator Info', 'Existing Plant Data.xlsx')

    ################
    # create generator technology definitions for non-renewable projects
    # (renewable ones were already added via the tracking_pv code)
    gen_info = data_frame_from_xlsx(gen_info_file, 'non_renewable_technology_info') \
        .T.set_index(0).T.set_index('technology')

    # convert from cost per MWh to cost per kWh
    gen_info['variable_o_m'] *= 0.001

    # add some fields
    gen_info['min_vintage_year'] = gen_info['build_year']
    gen_info['max_age_years'] = gen_info['retirement year'] - gen_info['build_year']
    gen_info['base_year'] = get_named_cell_from_xlsx(gen_info_file, named_range='base_year')

    # keep only basic generator info (dropping project-related fields)
    gen_info = gen_info.loc[:, 'unit_size':]

    # store generator info
    # note: the table should have been emptied and recreated by new_generator_info()
    # before calling this function
    gen_info.to_sql('generator_info', db_engine, if_exists='append')

    ################
    # create heat rate curves
    heat_rate_curves = data_frame_from_xlsx(gen_info_file, 'heat_rate_info')
    # place dummy values in the first level of the index; otherwise NaNs match any slice
    heat_rate_curves.loc[0, :] = heat_rate_curves.loc[0, :].fillna('x')
    # create the column index
    heat_rate_curves = heat_rate_curves.T.set_index([0, 1]).T
    # create the row index
    heat_rate_curves = heat_rate_curves.set_index(('x','technology'))
    heat_rate_curves.index.names=['technology']
    # get heat-rate specific info
    heat_rate_curves = heat_rate_curves[['PP', 'FC']].rename(
        columns={'PP':'output_mw', 'FC': 'fuel_consumption_mmbtu_per_h'}
    ).astype(float)

    # switch to database format
    heat_rate_curves = heat_rate_curves.stack()[['output_mw', 'fuel_consumption_mmbtu_per_h']]
    # don't use min/1/2/3/max labels
    heat_rate_curves.index = heat_rate_curves.index.droplevel(1)
    # sort rows appropriately (only matters for display)
    heat_rate_curves = heat_rate_curves.reset_index()
    heat_rate_curves = heat_rate_curves.sort_values(['technology', 'output_mw', 'fuel_consumption_mmbtu_per_h'])
    heat_rate_curves = heat_rate_curves.set_index('technology')
    # drop blank entries and treat the rest as floating point
    heat_rate_curves = heat_rate_curves.astype(float).dropna(axis=0, subset=['output_mw'])

    # store in database (should already be emptied and created by new_generator_info())
    heat_rate_curves.to_sql('part_load_fuel_consumption', db_engine, if_exists='append')

    ################
    # add existing non-renewable projects to project table
    # (renewable ones were already added via the wind and tracking_pv code)
    projects = data_frame_from_xlsx(gen_info_file, 'project_info').T.set_index(0).T
    # filter to include only non-renewable projects
    projects = projects.ix[~projects['technology'].isin(renewable_techs), :]

    # create columns not provided in xlsx file
    # don't allow new instances of existing thermal projects (these could
    # be allowed by creating a max_capacity column in the spreadsheet and specifying
    # future costs for those technologies elsewhere)
    projects['max_capacity'] = projects['proj_existing_cap']
    projects['interconnect_id'] = None
    projects['connect_distance_km'] = 0.0
    projects['connect_cost_per_mw'] = 0.0

    # set index and choose correct columns for database
    projects = projects.set_index(['load_zone', 'technology', 'site', 'orientation'])
    projects = projects[[
        'max_capacity', 'latitude', 'longitude', 'interconnect_id',
        'connect_distance_km', 'connect_cost_per_mw'
    ]]
    # insert non-renewable project definitions into project table
    projects.to_sql('project', db_engine, if_exists='append')

    ##############
    # create proj_existing_builds, holding construction dates and costs for existing projects;
    # this assigns solar projects to suitable resource tranches (which reduces the amount
    # of those tranches available for new construction)
    proj_build = data_frame_from_xlsx(gen_info_file, 'project_info').T.set_index(0).T

    # assign existing utility-scale renewable energy projects to the nearest site
    near_query = """
        select
            site, orientation,
            ((latitude-%(latitude)s)^2+(longitude - %(longitude)s)^2)^0.5 as dist
        from project
        where load_zone=%(load_zone)s and technology=%(technology)s
        order by 3
        limit 1;
    """
    for row in proj_build.itertuples():
        if row.technology in renewable_techs and row.technology != 'DistPV':
            # find the nearest project and assign this capacity to that
            # note: row is a namedtuple; vars() converts it to a dict (also works in 2.7):
            # https://docs.python.org/3.3/library/collections.html#collections.somenamedtuple._asdict
            nearest = pd.read_sql(sql=near_query, con=db_engine, params=vars(row))
            proj_build.loc[row.Index, ['site', 'orientation']] = nearest.loc[0, ['site', 'orientation']]

    # replace single DistPV project with several projects spread among the
    # better-than-average resources within that zone (e.g., south-facing roofs)

    # remove the DistPV rows from proj_build and keep them for further reference
    proj_build_dist_pv = proj_build[proj_build['technology']=='DistPV']
    proj_build = proj_build.drop(proj_build_dist_pv.index)

    # get a list of all better-than-average solar sites in each zone
    dist_pv_tranche_query = """
        WITH site_cap_factor AS (
            SELECT
                load_zone, technology, site, orientation, max_capacity as site_capacity,
                AVG(cap_factor) AS cap_factor
            FROM project JOIN cap_factor USING (project_id)
            WHERE technology='DistPV'
            GROUP BY 1, 2, 3, 4, 5
        ), zone_cap_factor AS (
            SELECT
                load_zone, technology,
                sum(cap_factor*site_capacity)/SUM(site_capacity) AS zone_cap_factor
                FROM site_cap_factor
                GROUP BY 1, 2
        ), good_sites as (
            SELECT *
            FROM site_cap_factor s JOIN zone_cap_factor z USING (load_zone, technology)
            WHERE cap_factor >= zone_cap_factor
        ), zone_good_capacity AS (
            SELECT load_zone, technology, SUM(site_capacity) AS zone_good_capacity
            FROM good_sites
            GROUP BY 1, 2
        )
        SELECT *
        FROM good_sites JOIN zone_good_capacity USING (load_zone, technology);
    """
    dist_pv_tranches = pd.read_sql(dist_pv_tranche_query, con=db_engine)

    # pair project templates with tranches based on load_zone and technology
    # (but not site or orientation)
    new_rows = (
        proj_build_dist_pv.drop(['site', 'orientation'], axis=1) \
        .merge(dist_pv_tranches, on=['load_zone', 'technology'], how='left')
    )
    # allocate existing capacity among tranches
    new_rows['proj_existing_cap'] = (
        new_rows['proj_existing_cap'] * new_rows['site_capacity'] / new_rows['zone_good_capacity']
    )
    # append matching columns to proj_build
    proj_build = proj_build.append(new_rows.reindex(columns=proj_build.columns))

    # lookup project_id's for existing projects
    proj_id = pd.read_sql(
        'SELECT project_id, load_zone, technology, site, orientation FROM project;',
        con=db_engine
    )
    proj_build = proj_build.merge(proj_id, how='left')
    proj_unmatched = proj_build[proj_build['project_id'].isnull()]
    if proj_unmatched.shape[0] > 0:
        print "="*70
        print "WARNING: The following existing projects were not found in the project table:"
        print proj_unmatched
        print 'See "{}" for details.'.format(gen_info_file)
        print "="*70

    # create/replace proj_existing_builds table (with appropriate formats for columns)
    proj_build['build_year'] = proj_build['build_year'].astype(int)
    proj_build = proj_build.set_index(['project_id', 'build_year'])
    proj_build = proj_build[['proj_existing_cap', 'proj_overnight_cost', 'proj_fixed_om']].astype(float)

    proj_build.to_sql('proj_existing_builds', db_engine, if_exists='replace')

    # make sure no projects are over-allocated
    # (may also prompt an error or infeasibility in SWITCH later)
    excess_allocation = pd.read_sql(
        """
            SELECT
                project.project_id, load_zone, technology, site, orientation, max_capacity,
                sum(proj_existing_cap) as proj_existing_cap
            FROM project JOIN proj_existing_builds USING (project_id)
            GROUP BY 1, 2, 3, 4, 5, 6
            HAVING sum(proj_existing_cap) > max_capacity;
        """,
        con=db_engine
    ).set_index('project_id')
    if excess_allocation.shape[0] > 0:
        print "="*70
        print "WARNING: The following projects have existing capacity greater than"
        print "the maximum possible capacity:"
        print excess_allocation
        print 'See "{}" for details.'.format(gen_info_file)
        print "="*70


def system_load():
    # TODO: extend to other load zones by adding more rows to the 
    # 'sales_forecast' region of the psip_data_file

    # get historical peak and average loads
    hist = pd.read_sql(
        sql="""
            SELECT
                load_zone, EXTRACT(year FROM date_time) as year_hist,
                MAX(system_load) as peak_hist, AVG(system_load) as avg_hist
            FROM system_load
            GROUP BY 1, 2;
        """,
        con=db_engine
    )
    # forecast peak and energy
    fore = data_frame_from_xlsx(psip_data_file, 'sales_forecast')
    fore = fore.T.set_index(0).T
    fore = fore.rename(columns={'year': 'year_fore'})
    # calculate scale factors for system_load_scale table
    sls = pd.merge(hist, fore, on='load_zone')
    sls['load_scen_id'] = load_scen_id
    sls['peak_fore'] = sls['underlying forecast (MW)'] + sls['energy efficiency (MW)']
    sls['avg_fore'] = (sls['underlying forecast (GWh)'] + sls['energy efficiency (GWh)'])/8.76
    sls['scale'] = (sls['peak_fore'] - sls['avg_fore']) / (sls['peak_hist'] - sls['avg_hist'])
    sls['offset'] = sls['peak_fore'] - sls['scale'] * sls['peak_hist']

    # put into standard order, drop unneeded columns, convert to the right types for the database
    db_columns = [
        'load_zone', 'load_scen_id', 'year_hist', 'year_fore',
        'peak_hist', 'peak_fore', 'avg_hist', 'avg_fore', 'scale', 'offset'
    ]
    system_load_scale = pd.DataFrame()
    for c in db_columns:
        if c in ['load_zone', 'load_scen_id']:
            system_load_scale[c] = sls[c].astype(str)
        elif c in ['year_hist', 'year_fore']:
            system_load_scale[c] = sls[c].astype(int)
        else:
            system_load_scale[c] = sls[c].astype(float)
    system_load_scale.set_index(db_columns[:4], inplace=True)
    # store data
    execute("DELETE FROM system_load_scale WHERE load_scen_id=%s;", (load_scen_id,))
    system_load_scale.to_sql('system_load_scale', db_engine, if_exists='append')
    
    # create another forecast with peak and average loads from 2007, carried through to the future
    execute("""
        CREATE TEMPORARY TABLE tsls AS 
            SELECT * FROM system_load_scale WHERE load_scen_id=%s;
        UPDATE tsls a 
            SET peak_fore=b.peak_hist, avg_fore=b.avg_hist, load_scen_id='flat_2007'
            FROM tsls b
            WHERE b.year_hist=2007 and b.year_fore=2045;
        UPDATE tsls 
            SET scale = (peak_fore - avg_fore) / (peak_hist - avg_hist);
        UPDATE tsls 
            SET "offset" = peak_fore - scale * peak_hist;
        INSERT INTO system_load_scale SELECT * FROM tsls;
        DROP TABLE tsls;
    """, (load_scen_id,))
        
def interconnect():
    # also see database/build_database/shared_tables.py for code to fill in
    # project.interconnect_id, project.connect_distance_km and project.connect_cost_per_mw
    # based on this table
    # note: we could eventually add interconnect-specific connection costs here,
    # to be used instead of generic project interconnection costs; in that case
    # the code in shared_tables.calculate_interconnect_costs() would also need
    # to be updated
    execute("""
        DROP TABLE IF EXISTS interconnect;
        CREATE TABLE interconnect (
            interconnect_id integer PRIMARY KEY NOT NULL,
            county text,
            latitude float,
            longitude float
        );
        ALTER TABLE interconnect OWNER TO admin;
        -- At some point interconnect was filled in with the equivalent of the
        -- following command. The original code is missing, but these appear to be
        -- the population-weighted centers of each county.
        INSERT INTO interconnect (interconnect_id, county, latitude, longitude) VALUES
            (1, 'Honolulu', 21.372464, -157.913673),
            (2, 'Hawaii', 19.672837, -155.421895),
            (3, 'Maui', 20.863747, -156.493816),
            (4, 'Kauai', 22.021022, -159.442112),
            (5, 'Kalawao', 21.188495, -156.979972);
    """)


if __name__ == "__main__":
    main()
