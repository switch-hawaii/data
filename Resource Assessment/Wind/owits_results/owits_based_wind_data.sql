-- note: these are clusters and production profiles that Paritosh Das prepared as
-- documented in http://evtc.fsec.ucf.edu/publications/documents/HI-13-16.pdf.
-- It is not clear whether we have code to reproduce this anywhere, or if it was
-- done ad hoc and/or Paritosh lost it.
-- The commands below were used to copy these tables from the redr switch database;
-- there is code in the import_data.py script for that to create these from somewhat
-- older tables also in the switch database.

\copy (select * from project where technology='OnshoreWind' order by 1) to wind_project.csv with CSV header
\copy (select c.* from project p natural join cap_factor c where technology='OnshoreWind' order by 1, 2) to wind_cap_factor.csv with CSV header
