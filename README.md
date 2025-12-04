Romania GHG Dashboard (Streamlit)
=================================

This project extracts greenhouse gas (GHG) emission data for Romania (national, regional, and cities) from the Excel and NetCDF files in this repository and renders multiple interactive dashboards using Streamlit.

What it does
------------
- Preprocesses all relevant files and stores only Romania-related data in a separate SQLite DB at `data/ro_ghg.db`.
- Provides dashboards for:
  - National totals over time (EDGAR)
  - NUTS2 regional totals (with ranking + timeseries)
  - City totals and sectoral timeseries (UCDB)
  - Spatial grid map over Romania (EDGAR CO₂ 2024 totals)

Prerequisites
-------------
- Use the provided virtual environment at `env/` (Conda-based). The setup commands below assume it exists and is active via explicit path usage.

Preprocessing
-------------
Run the preprocessing script to create/update the SQLite database. It scans the repository for the known datasets and extracts Romania-only subsets.

```
env/bin/python scripts/preprocess.py
```

The database will be created at `data/ro_ghg.db`.

Run the app
-----------
Launch the Streamlit app after preprocessing:

```
env/bin/streamlit run app/streamlit_app.py
```

Notes
-----
- The app reads data only from `data/ro_ghg.db`; if you change/replace source files, re-run the preprocessing step.
- City totals come from `allcountries_onlycities_summary.xlsx`; sectoral city timeseries come from `EDGAR/EDGAR_emiss_on_UCDB_2024.xlsx` (only `GWP_100_AR5_GHG` metrics are stored to keep the DB light).
- National and regional (NUTS2) series are extracted from the EDGAR annual totals and NUTS2 datasets.
- The grid map uses the EDGAR 2024 CO₂ totals NetCDF and plots a 2D heatmap over Romania.
