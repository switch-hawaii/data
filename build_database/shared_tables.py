from __future__ import absolute_import
from util import execute

queries = {}

# note: some of the project data could go in a separate site table,
# but we keep it all in one place for now for simplicity
# ??? can we add existing projects to this table too? (no reason not to;
# may need to add a flag indicating whether more capacity can be built in
# each project.)
# note: we assume gen_capacity_limit_mw indicates the max amount of each technology
# if that is the only thing built at this site; if multiple projects
# are built on the same site, we require sum(Build[site, tech]/gen_capacity_limit_mw[site, tech]) <= 1.
# note: we use double precision instead of real to avoid rounding errors when comparing
# gen_capacity_limit_mw to gen_build_predetermined (which ends up with double precision) and generally
# to maintain consistency throughout the work
queries[("projects", "create_table")] = """
    CREATE TABLE IF NOT EXISTS projects (
        project_id SERIAL PRIMARY KEY,
        load_zone VARCHAR(20),
        technology VARCHAR(50),
        site VARCHAR(20),
        orientation VARCHAR(5),
        gen_capacity_limit_mw DOUBLE PRECISION,
        latitude DOUBLE PRECISION,
        longitude DOUBLE PRECISION,
        interconnect_id INT,
        connect_distance_km DOUBLE PRECISION,
        spur_line_cost_per_mw DOUBLE PRECISION
    );
"""
queries[("variable_capacity_factors", "create_table")] = """
    CREATE TABLE IF NOT EXISTS variable_capacity_factors (
        project_id INT NOT NULL,
        date_time TIMESTAMP WITH TIME ZONE,
        cap_factor REAL
    );
"""
queries[("generator_info", "create_table")] = """
    CREATE TABLE IF NOT EXISTS generator_info (
        tech_scenario VARCHAR(30) NOT NULL,
        technology VARCHAR(50),
        min_vintage_year INT,
        gen_unit_size DOUBLE PRECISION,
        gen_min_build_capacity DOUBLE PRECISION,
        substation_cost_per_kw DOUBLE PRECISION,
        variable_o_m DOUBLE PRECISION,
        gen_energy_source VARCHAR(20),
        gen_full_load_heat_rate DOUBLE PRECISION,
        gen_max_age INT,
        gen_forced_outage_rate DOUBLE PRECISION,
        gen_scheduled_outage_rate DOUBLE PRECISION,
        gen_is_variable INT,
        resource_limited INT,
        distributed INT,
        gen_is_baseload INT,
        must_run INT,
        non_cycling INT,
        gen_is_cogen INT,
        gen_min_uptime DOUBLE PRECISION,  -- hours
        gen_min_downtime DOUBLE PRECISION,  -- hours
        gen_startup_fuel DOUBLE PRECISION,  -- MMBtu per MW
        base_year INT,
        gen_storage_efficiency DOUBLE PRECISION,
        gen_storage_energy_to_power_ratio DOUBLE PRECISION,
        gen_storage_max_cycles_per_year DOUBLE PRECISION
    );
"""

# queries[("variable_capacity_factors", "create_indexes")] = """
#     DO $$
#     BEGIN
#         BEGIN
#             ALTER TABLE variable_capacity_factors
#                 ADD CONSTRAINT pt PRIMARY KEY (project_id, date_time),
#                 ADD CONSTRAINT tp UNIQUE (date_time, project_id)
#         EXCEPTION
#             WHEN duplicate_object THEN NULL; -- ignore if index exists already
#         END;
#     END $$;
# """

# note: if this reports 'relation "pt" already exists', it probably means an index
# named pt is already attached to an old (renamed) version of variable_capacity_factors. That can
# be viewed via "select * from pg_indexes where indexname='pt';" and renamed via
# "alter index pt rename to pt_2018_07_23;"
queries[("variable_capacity_factors", "create_indexes")] = """
    ALTER TABLE variable_capacity_factors
        ADD CONSTRAINT pt PRIMARY KEY (project_id, date_time),
        ADD CONSTRAINT tp UNIQUE (date_time, project_id);
"""
queries[("variable_capacity_factors", "drop_indexes")] = """
    ALTER TABLE variable_capacity_factors
        DROP CONSTRAINT IF EXISTS pt,
        DROP CONSTRAINT IF EXISTS tp;
"""

queries[("periods", "create_table")] = """
    CREATE TABLE IF NOT EXISTS periods (
        time_sample character varying(40) NOT NULL,
        period bigint NOT NULL,
        period_end integer
    );
    ALTER TABLE periods
        DROP CONSTRAINT IF EXISTS periods_pkey,
        ADD CONSTRAINT periods_pkey PRIMARY KEY (time_sample, period);
"""
queries[("timeseries", "create_table")] = """
    CREATE TABLE timeseries (
        period bigint,
        timeseries bigint NOT NULL,
        month_of_year integer,
        date date,
        hours_in_sample double precision,
        time_sample character varying(40) NOT NULL,
        ts_num_tps integer,
        ts_duration_of_tp double precision,
        ts_scale_to_period double precision
    );
    ALTER TABLE timeseries
        DROP CONSTRAINT IF EXISTS timeseries_pkey,
        ADD CONSTRAINT timeseries_pkey PRIMARY KEY (time_sample, timeseries);
"""
queries[("timepoints", "create_table")] = """
    CREATE TABLE timepoints (
        timeseries bigint NOT NULL,
        timepoint bigint NOT NULL,
        hour_of_day integer,
        date_time timestamp with time zone NOT NULL,
        time_sample character varying(40)
    );
    ALTER TABLE timepoints
        DROP CONSTRAINT IF EXISTS timepoints_pkey,
        ADD CONSTRAINT timepoints_pkey PRIMARY KEY
            (time_sample, timeseries, timepoint);
"""


def create_table(table):
    execute(queries[(table, "create_table")])

def create_indexes(table):
    if (table, "create_indexes") in queries:
        execute(queries[(table, "create_indexes")])

def drop_indexes(table):
    if (table, "drop_indexes") in queries:
        execute(queries[(table, "drop_indexes")])

def create_database():
    # can also get current definitions via pg_dump -s > postgresql.sql
    print("""
        I created a self-signed certificate and private key (server.crt and
        server.key) in the postgresql data directory following instructions at
        http://www.postgresql.org/docs/9.3/static/ssl-tcp.html . Then I added
        "ssl=on" to postgresql.conf and changed the "host" entries in
        pg_hba.conf to "hostssl" to only allow ssl connections (per
        http://www.postgresql.org/docs/9.3/static/ssl-tcp.html and
        http://www.postgresql.org/docs/9.3/static/auth-pg-hba-conf.html). Then I
        restarted the server. So now the postgresql server on
        redr.eng.hawaii.edu should only be accepting ssl connections, and
        should be requiring MD5 encryption of passwords.

        I tested the new setup with my existing client software (pgAdmin3, psql
        and psycopg2, as used by get_scenario_data.py), and it seems to be
        working fine.

        By default in postgresql: any user can access any database (connect,
        maybe more); this is not shown as a command applied to the database in
        pgadmin3, it just seems to be the default access level when new schemas
        are created (including public), the command "GRANT ALL ON SCHEMA public
        TO public;" is automatically run. this means any user can connect to any
        database

        Note: all of the commands below can be done by a superuser account; some
        can be done after SET ROLE switch_hawaii_owner.

        To manage access to the databases, run the following commands for each
        database (and any new databases), then give all allowed users explicit
        access to the database (further below). These examples work for
        switch_hawaii. Some info can be found here.

        REVOKE ALL ON DATABASE switch_hawaii FROM public;   -- all accounts have full access by default
        REVOKE ALL ON SCHEMA public FROM PUBLIC;    -- public is granted full access when schema is created

        Create read-only and owner access roles for the database:

        CREATE ROLE switch_hawaii_owner;
        GRANT ALL ON DATABASE switch_hawaii TO switch_hawaii_owner;
        GRANT ALL ON SCHEMA public TO switch_hawaii_owner;
        GRANT ALL ON ALL TABLES IN SCHEMA public TO switch_hawaii_owner;

        CREATE ROLE switch_hawaii_reader;
        GRANT CONNECT ON DATABASE switch_hawaii TO switch_hawaii_reader;
        GRANT USAGE ON SCHEMA public TO switch_hawaii_reader;
        GRANT SELECT ON ALL TABLES IN SCHEMA public TO switch_hawaii_reader;

        Then you can create owners and read-only users as follows:
        CREATE ROLE <user> LOGIN PASSWORD '<password>;
        GRANT switch_hawaii_owner TO <user>;
        GRANT switch_hawaii_reader TO <user>;

        -- The lines below are one way to give reading and quasi-owning privileges to the right groups when new tables
        -- are created. However, they have to be run once for each account that may create tables, which is a hassle,
        -- and even then, new tables remain owned by whoever created them, so other user members of the xxx_owner
        -- group can't drop or alter them. So it is better just to use the event trigger shown below.
        -- SET ROLE <user>;
        -- ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT ALL ON TABLES TO switch_hawaii_owner;
        -- ALTER DEFAULT PRIVILEGES IN SCHEMA public GRANT SELECT ON TABLES TO switch_hawaii_reader;

        -- This event trigger makes sure every new table is owned by the group rather than the user, so any other group
        -- member can drop or alter the table. They also enable giving read access to the read-only users.
        -- This is based on code at https://blog.hagander.net/setting-owner-at-create-table-237/.

        CREATE OR REPLACE FUNCTION trg_create_set_owner()
         RETURNS event_trigger
         LANGUAGE plpgsql
        AS $$
        DECLARE
          obj record;
        BEGIN
          FOR obj IN SELECT tablename FROM pg_tables WHERE tableowner = current_user LOOP
          -- postgresql 9.5 or later:
          -- SELECT table_name FROM pg_event_trigger_ddl_commands() WHERE command_tag='CREATE TABLE'
            EXECUTE format('GRANT SELECT ON %s TO switch_hawaii_reader', obj.tablename);
            EXECUTE format('ALTER TABLE %s OWNER TO switch_hawaii_owner', obj.tablename);
          END LOOP;
        END;
        $$;

        DROP EVENT TRIGGER trg_create_set_owner;
        CREATE EVENT TRIGGER trg_create_set_owner
         ON ddl_command_end
         WHEN tag IN ('CREATE TABLE', 'CREATE TABLE AS')
         EXECUTE PROCEDURE trg_create_set_owner();

        Then users should install software to access the server (psycopg2 and/or
        sqlalchemy for Python, RPostgres for R, PgAdmin3/PgAdmin4/psql for
        direct access). Then they should create a ~/.pgpass file with a line
        like redr.eng.hawaii.edu:5432:*:<user name>:<password> (and chmod
        0600).
    """)
