from textwrap import dedent

switch_db = 'switch'
pg_host = 'redr.eng.hawaii.edu'

try:
    import psycopg2
except ImportError:
    print "This module requires the psycopg2 module to access the postgresql database."
    print "Please execute 'sudo pip install psycopg2' or 'pip install psycopg2' (Windows)."
    raise

try:
    # note: the connection gets created when the module loads and never gets closed 
    # (until presumably python exits)
    con = psycopg2.connect(database=switch_db, host=pg_host, sslmode='require')

    # note: we don't autocommit because it makes executemany() very slow; 
    # instead we call con.commit() after each query
    con.autocommit = False

    # note: con and cur stay open until the module goes out of scope
    cur = con.cursor()    

except psycopg2.OperationalError:
    print dedent("""
        ############################################################################################
        Error while connecting to {db} database on postgresql server {server}.
        Please ensure that your user name on the local system is the same as your postgresql user 
        name or there is a local PGUSER environment variable set with your postgresql user name.
        There should also be a line like "*:*:{db}:<username>:<password>" in ~/.pgpass or 
        %APPDATA%\postgresql\pgpass.conf (Windows). On Unix systems, .pgpass should be chmod 0600.
        See http://www.postgresql.org/docs/9.3/static/libpq-pgpass.html for more details.
        ############################################################################################
        """.format(db=switch_db, server=pg_host))
    raise


def execute(query, *args, **kwargs):
    return _execute(query, False, *args, **kwargs)

def executemany(query, *args, **kwargs):
    return _execute(query, True, *args, **kwargs)

def _execute(query, many, *args, **kwargs):
    q = dedent(query)
    func = cur.executemany if many else cur.execute
    print q
    try:
        func(q, *args, **kwargs)
        con.commit()
        return cur
    except:
        con.rollback()
        raise
    

