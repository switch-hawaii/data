This repository contains data and code used to create the data warehouse for Switch
models of the Oahu power system.

The main script is `build_database/import_data.py`. This will pull in data from
many input files and directories and create a postgresql database called
`switch_hawaii`. See `build_database/import_data.py` for a description of the
configuration files you need to connect to the database.

Most of the inputs for `import_data.py` are included in this repository, but
some large input files need to be downloaded from public sources. If you get a
"file not found" error, look for a "data sources.txt" file, download script or
similar file in the relevant directory, and follow the instructions to download
those files. There are also some files titled `steps to create ... .txt` that
describe steps to follow with GIS software to perform land use screening. The
results from this screening are already stored in the repository for use by
`import_data.py`, so you do not need to re-run them. But you can review those
instructions to see exactly how the screening was done, and you can modify them
and follow the new instructions if you want to change the screening rules.

After the `switch_hawaii` database is constructed, you can use
`switch_model.hawaii.scenario_data` (part of the main Switch software
distribution) to extract data for the particular dates and cost scenarios needed
to run an individual model. See `get_scenario_data.py` scripts in various model
repositories on https://github.com/switch-hawaii/ for examples of how to do
this.

These are some of the important input files used to create the data warehouse:

- `Generator Info/Existing Plant Data.xlsx`
  - data describing the capabilities of existing power plants
- `Generator Info/PSIP 2016-12 ATB 2020 generator data.xlsx`
  - data describing new renewable projects that could be developed, as well as
    HECO sales forecasts
  - note that HECO DER forecasts and near-term construction forecasts are
    shown in `switch_model.hawaii.heco_outlook_2020_08` and
    `switch_model.hawaii.heco_outlook_2020_08` in the main Switch repository
- `EV Adoption/IGP adoption forecast 2020-05-13.xlsx`
  - EV adoption projections based on Hawaiian Electric IGP docket, used for
    Ulupono scenarios in the PBR docket in 2020
- `EV Adoption/EV projections with buses.xlsx`
  - EV adoption projections used for academic studies with Switch
- `EV Adoption/EoT Roadmap Nonmanaged Charging in 2030.xlsx`
  - business-as-usual EV charging shapes used to define the
    `EoT_2018_avg` schedule used for Ulupono scenarios in the PBR docket in 2020
- `EV Adoption/ev_hourly_charge_profile.tsv`
  - business-as-usual EV charging shapes to define the `das_2015` charging
    schedule used for some earlier studies
- `Generator Info/build_database/import_data.py`
  - imports all data into `switch_hawaii` data warehouse
- `build_database/solar_resources.py`
  - contains functions to calculate hourly performance of rooftop and
    utility-scale solar (called by `import_data.py`)
