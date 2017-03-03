from util import execute

connect_cost_per_mw_km = 1000.0     # TODO: specify this somewhere better

queries = {}

# note: some of the project data could go in a separate site table,
# but we keep it all in one place for now for simplicity
# ??? can we add existing projects to this table too? (no reason not to;
# may need to add a flag indicating whether more capacity can be built in
# each project.)
# note: we assume max_capacity indicates the max amount of each technology
# if that is the only thing built at this site; if multiple projects
# are built on the same site, we require sum(Build[site, tech]/max_capacity[site, tech]) <= 1.
# note: we use double precision instead of real to avoid rounding errors when comparing
# max_capacity to proj_existing_builds (which ends up with double precision) and generally
# to maintain consistency throughout the work
queries[("project", "create_table")] = """
    CREATE TABLE IF NOT EXISTS project (
        project_id SERIAL PRIMARY KEY,
        load_zone VARCHAR(20),
        technology VARCHAR(50),
        site VARCHAR(20),
        orientation VARCHAR(5),
        max_capacity DOUBLE PRECISION,
        latitude DOUBLE PRECISION,
        longitude DOUBLE PRECISION,
        interconnect_id INT,
        connect_distance_km DOUBLE PRECISION,
        connect_cost_per_mw DOUBLE PRECISION
    );
    ALTER TABLE project OWNER TO admin;
"""
queries[("cap_factor", "create_table")] = """
    CREATE TABLE IF NOT EXISTS cap_factor (
        project_id INT NOT NULL,
        date_time TIMESTAMP WITH TIME ZONE,
        cap_factor REAL
    );
    ALTER TABLE cap_factor OWNER TO admin;
"""
queries[("generator_info", "create_table")] = """
    CREATE TABLE IF NOT EXISTS generator_info (
        technology VARCHAR(50),
        min_vintage_year INT,
        unit_size DOUBLE PRECISION,
        connect_cost_per_kw_generic DOUBLE PRECISION,
        fixed_o_m DOUBLE PRECISION,
        variable_o_m DOUBLE PRECISION,
        fuel VARCHAR(20),
        heat_rate DOUBLE PRECISION,
        max_age_years INT,
        forced_outage_rate DOUBLE PRECISION,
        scheduled_outage_rate DOUBLE PRECISION,
        intermittent INT,
        resource_limited INT,
        distributed INT,
        baseload INT,
        must_run INT,
        non_cycling INT,
        cogen INT,
        base_year INT
    );
    ALTER TABLE generator_info OWNER TO admin;
"""

# queries[("cap_factor", "create_indexes")] = """
#     DO $$
#     BEGIN
#         BEGIN
#             ALTER TABLE cap_factor
#                 ADD CONSTRAINT pt PRIMARY KEY (project_id, date_time),
#                 ADD CONSTRAINT tp UNIQUE (date_time, project_id)
#         EXCEPTION
#             WHEN duplicate_object THEN NULL; -- ignore if index exists already
#         END;
#     END $$;
# """
queries[("cap_factor", "create_indexes")] = """
    ALTER TABLE cap_factor
        ADD CONSTRAINT pt PRIMARY KEY (project_id, date_time),
        ADD CONSTRAINT tp UNIQUE (date_time, project_id)
"""
queries[("cap_factor", "drop_indexes")] = """
    ALTER TABLE cap_factor 
        DROP CONSTRAINT IF EXISTS pt,
        DROP CONSTRAINT IF EXISTS tp
"""

def create_table(table):
    execute(queries[(table, "create_table")])

def create_indexes(table):
    if (table, "create_indexes") in queries:
        execute(queries[(table, "create_indexes")])

def drop_indexes(table):
    if (table, "drop_indexes") in queries:
        execute(queries[(table, "drop_indexes")])

def calculate_interconnect_costs():
    """Choose closest interconnect location to each project, and calculate distance to it.
    Also calculate connect_cost_per_mw based on distance and generic connection cost for
    each technology.
    note: this could eventually be updated to use interconnect-specific costs, where
    provided, instead of generic project interconnect costs; in that case, code that
    creates the interconnect table in import_data.py would need to be updated.
    """
    execute("""
        WITH distances as (
            select p.project_id, i.interconnect_id, 
                -- haversine distance formula, radius of earth = 6371 km
                2 *  6371 * sqrt(
                    pow(sin(radians((i.latitude - p.latitude)/2)), 2) 
                    + cos(radians(p.latitude)) * cos(radians(i.latitude)) 
                        * pow(sin(radians((i.longitude - p.longitude)/2)), 2))
                as distance
                from project p, interconnect i
                where p.latitude is not null and p.longitude is not null
        ), closest as (
            select project_id, min(distance) as distance 
                from distances group by 1
        ), neighbor as (
            select c.project_id, d.interconnect_id, c.distance
                from closest c join distances d using (project_id, distance)
            -- note, this may return multiple interconnects with the same distance
            -- but that is rare, and one will be chosen arbitrarily in the update query
        )
        update project p
            set interconnect_id = n.interconnect_id, 
                connect_distance_km = n.distance
            from neighbor n
            where n.project_id = p.project_id;
    """)
    execute("""
        update project p
            set connect_cost_per_mw = 
                1000 * coalesce(connect_cost_per_kw_generic, 0)
                    + %(connect_cost_per_mw_km)s * coalesce(connect_distance_km, 0)
            from generator_info g
            where g.technology=p.technology;
    """, dict(connect_cost_per_mw_km=connect_cost_per_mw_km))
