#!/usr/bin/python
import csv
import datetime
import os
import sys
import psycopg2
import time


try:

    con = psycopg2.connect(database='switch', host='switch.eng.hawaii.edu', port='5432', user='deepakc_super', password='myPassword')
    completeFilePath = "C:/Users/HP USER/Downloads/form714-database (1)/Form 714 Export/Part 3 Schedule 2 - Planning Area Hourly Demand.csv"

    cur = con.cursor()
    
    cur.execute("""DROP TABLE IF EXISTS "timport"""")

    query = """CREATE TABLE "timport"(load_zone text, date_time timestamp with time zone, rowdata text)"""
    cur.execute(query)

    f = open(completeFilePath, 'r')
    query = """COPY "timport"(rowdata) FROM STDIN WITH CSV HEADER DELIMITER AS '\t' """
    cur.copy_expert(query, f)
    con.commit()

    
    
    if f:
        f.close()    
    
finally:
    if con:
        con.close()
    
