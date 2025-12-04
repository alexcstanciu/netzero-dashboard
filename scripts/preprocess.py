#!/usr/bin/env python3
"""
Preprocess GHG datasets to extract only Romania (national + cities) and
store them into a local SQLite database for the Streamlit dashboards.

Inputs (expected paths within the repo):
 - allcountries_summary.xlsx
 - allcountries_onlycities_summary.xlsx
 - EDGAR/EDGAR_emiss_on_UCDB_2024.xlsx (sheet: EDGAR_emiss_on_UCDB_2024)
 - EDGAR/Annual_Totals/EDGAR_AR5_GHG_1970_2024/EDGAR_AR5_GHG_1970_2024.xlsx (sheet: TOTALS BY COUNTRY)
 - EDGAR/EDGAR_2025_total_GHG_AR5_NUTS2_by_country_and_sector_1990-2024/EDGAR_2025_total_GHG_AR5_NUTS2_by_country_and sector_1990-2024.xlsx (sheets: TOTALS by NUTS2, TOTALS by NUTS2 and Sector)
 - EDGAR/Annual_gridmaps/EDGAR_2025_GHG_CO2_2024_TOTALS_emi_nc/EDGAR_2025_GHG_CO2_2024_TOTALS_emi.nc

Output:
 - data/ro_ghg.db (SQLite database)

Requires: pandas, openpyxl, xarray, netCDF4 or h5netcdf
"""
from __future__ import annotations

import os
import sqlite3
from pathlib import Path
from typing import List, Tuple

import pandas as pd

# Optional imports for NetCDF; allow script to run even if missing
try:
    import xarray as xr  # type: ignore
except Exception:
    xr = None  # type: ignore
try:
    import geopandas as gpd  # type: ignore
except Exception:
    gpd = None  # type: ignore
try:
    from shapely.geometry import Point, Polygon  # type: ignore
except Exception:
    Point = None  # type: ignore
    Polygon = None  # type: ignore


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
DB_PATH = DATA_DIR / "ro_ghg.db"


def ensure_db() -> sqlite3.Connection:
    DATA_DIR.mkdir(exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    # Improve write performance
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    return conn


def normalize_country_name(s: str | None) -> str:
    return (s or "").strip().lower()


def process_allcountries_summaries(conn: sqlite3.Connection) -> None:
    # Cities-only file
    cities_path = ROOT / "allcountries_onlycities_summary.xlsx"
    if cities_path.exists():
        dfc = pd.read_excel(cities_path, sheet_name=0)
        dfc.columns = [str(c).strip() for c in dfc.columns]
        if "Country" in dfc.columns:
            mask = dfc["Country"].astype(str).str.lower().eq("romania")
            dfc_ro = dfc.loc[mask].copy()
            # Unit normalization: "Total (t CO2)" is tonnes; add kt and Mt for convenience
            if "Total (t CO2)" in dfc_ro.columns:
                dfc_ro["total_tonnes"] = pd.to_numeric(dfc_ro["Total (t CO2)"], errors="coerce")
                dfc_ro["total_kt"] = dfc_ro["total_tonnes"] / 1_000.0
                dfc_ro["total_mt"] = dfc_ro["total_tonnes"] / 1_000_000.0
            # Round estimated population to whole persons (nullable integer)
            if "Est. Population" in dfc_ro.columns:
                pop = pd.to_numeric(dfc_ro["Est. Population"], errors="coerce")
                dfc_ro["Est. Population"] = pop.round().astype("Int64")
            dfc_ro.to_sql("cities_only_summary", conn, if_exists="replace", index=False)
            print(f"Stored {len(dfc_ro)} rows to cities_only_summary")
        else:
            print("[warn] cities summary: 'Country' column not found; skipped")
    else:
        print("[info] allcountries_onlycities_summary.xlsx not found; skipping")

    # All regions (may include national + admin levels)
    reg_path = ROOT / "allcountries_summary.xlsx"
    if reg_path.exists():
        dfr = pd.read_excel(reg_path, sheet_name=0)
        dfr.columns = [str(c).strip() for c in dfr.columns]
        if "Country" in dfr.columns:
            mask = dfr["Country"].astype(str).str.lower().eq("romania")
            dfr_ro = dfr.loc[mask].copy()
            if "Total (t CO2)" in dfr_ro.columns:
                dfr_ro["total_tonnes"] = pd.to_numeric(dfr_ro["Total (t CO2)"], errors="coerce")
                dfr_ro["total_kt"] = dfr_ro["total_tonnes"] / 1_000.0
                dfr_ro["total_mt"] = dfr_ro["total_tonnes"] / 1_000_000.0
            dfr_ro.to_sql("regions_summary", conn, if_exists="replace", index=False)
            print(f"Stored {len(dfr_ro)} rows to regions_summary")
        else:
            print("[warn] allcountries_summary: 'Country' column not found; skipped")
    else:
        print("[info] allcountries_summary.xlsx not found; skipping")


def process_edgar_ucdb(conn: sqlite3.Connection) -> None:
    path = ROOT / "EDGAR/EDGAR_emiss_on_UCDB_2024.xlsx"
    if not path.exists():
        print("[info] EDGAR_emiss_on_UCDB_2024.xlsx not found; skipping")
        return
    # Read only the main sheet
    df = pd.read_excel(path, sheet_name="EDGAR_emiss_on_UCDB_2024")
    # Filter to Romania
    if "UC_country" not in df.columns:
        print("[warn] UCDB: 'UC_country' column not found; skipped")
        return
    df = df[df["UC_country"].astype(str).str.strip().str.lower().eq("romania")].copy()
    # Keep only GHG (AR5) columns for compactness; format: EMI_GWP_100_AR5_GHG_<Sector>_<Year>
    ghg_cols = [c for c in df.columns if isinstance(c, str) and c.startswith("EMI_GWP_100_AR5_GHG_")]
    base_cols = [c for c in ["ID_UC_G0", "UC_name", "UC_country"] if c in df.columns]
    df_small = df[base_cols + ghg_cols].copy()
    # Unpivot to long format
    long = df_small.melt(id_vars=base_cols, var_name="metric", value_name="value_t")
    # Parse sector and year from metric
    def parse_metric(s: str) -> Tuple[str, int | None]:
        # EMI_GWP_100_AR5_GHG_<Sector>_<Year> -> parts length 7
        parts = s.split("_")
        sector = parts[5] if len(parts) >= 7 else "Unknown"
        year = int(parts[6]) if len(parts) >= 7 and parts[6].isdigit() else None
        return sector, year

    parsed = long["metric"].astype(str).apply(parse_metric)
    long["sector"] = parsed.apply(lambda x: x[0])
    long["year"] = parsed.apply(lambda x: x[1])
    long.drop(columns=["metric"], inplace=True)
    # Some rows may not parse; drop missing years
    long = long.dropna(subset=["year"]).copy()
    long["year"] = long["year"].astype(int)
    # Convert tonnes to kt for readability
    long["value_kt"] = pd.to_numeric(long["value_t"], errors="coerce") / 1_000.0
    long.to_sql("edgar_ucdb_cities_ghg", conn, if_exists="replace", index=False)
    print(f"Stored {len(long)} rows to edgar_ucdb_cities_ghg")


def process_edgar_country_totals(conn: sqlite3.Connection) -> None:
    path = ROOT / "EDGAR/Annual_Totals/EDGAR_AR5_GHG_1970_2024/EDGAR_AR5_GHG_1970_2024.xlsx"
    sheet = "TOTALS BY COUNTRY"
    if not path.exists():
        print("[info] EDGAR_AR5_GHG_1970_2024.xlsx not found; skipping country totals")
        return
    df = pd.read_excel(path, sheet_name=sheet)
    # Find data rows having expected column names present; after header preamble there is a row with Name, Substance, Y_1970... etc
    # Pandas should have recognized headers, so look for 'Name' and 'Substance' presence
    if not {"Name", "Substance"}.issubset(df.columns):
        # Attempt to re-read without header and locate header row
        raw = pd.read_excel(path, sheet_name=sheet, header=None)
        header_idx = raw.index[raw.apply(lambda r: r.astype(str).str.contains("Country_code_A3").any(), axis=1)].tolist()
        if header_idx:
            df = pd.read_excel(path, sheet_name=sheet, header=header_idx[0])
        else:
            # Fallback: find first row where first col equals 'IPCC_annex'
            header_idx = raw.index[raw.iloc[:, 0].astype(str).eq("IPCC_annex")].tolist()
            if header_idx:
                df = pd.read_excel(path, sheet_name=sheet, header=header_idx[0])
    # Filter Romania + GHG compound
    if not {"Name", "Substance"}.issubset(df.columns):
        print("[warn] Country totals sheet structure unexpected; skipped")
        return
    rom = df[(df["Name"].astype(str).str.lower().eq("romania")) & (df["Substance"].astype(str).str.contains("GWP_100_AR5_GHG"))].copy()
    if rom.empty:
        print("[warn] No Romania rows found in country totals; skipped")
        return
    # Keep Y_XXXX columns
    ycols = [c for c in rom.columns if isinstance(c, str) and c.startswith("Y_")]
    out = rom.melt(id_vars=["Name", "Substance"], value_vars=ycols, var_name="year_label", value_name="value_gg")
    out["year"] = out["year_label"].str.replace("Y_", "", regex=False).astype(int)
    out.drop(columns=["year_label"], inplace=True)
    # Convert Gg to Mt (1 Gg = 0.001 Mt)
    out["value_mtco2e"] = pd.to_numeric(out["value_gg"], errors="coerce") / 1_000.0
    out.rename(columns={"Name": "country", "Substance": "compound"}, inplace=True)
    out.to_sql("edgar_country_totals", conn, if_exists="replace", index=False)
    print(f"Stored {len(out)} rows to edgar_country_totals")


def process_edgar_nuts2(conn: sqlite3.Connection) -> None:
    path = ROOT / "EDGAR/EDGAR_2025_total_GHG_AR5_NUTS2_by_country_and_sector_1990-2024/EDGAR_2025_total_GHG_AR5_NUTS2_by_country_and sector_1990-2024.xlsx"
    sheet = "TOTALS by NUTS2"
    if not path.exists():
        print("[info] NUTS2 workbook not found; skipping")
        return
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
    # Identify data rows: they start with 'GWP_100_AR5_GHG' in col0
    rows = raw[raw.iloc[:, 0].astype(str).str.contains("GWP_100_AR5_GHG", na=False)].copy()
    if rows.empty:
        print("[warn] No NUTS2 data rows detected; skipped")
        return
    # Columns: [Compound, A3, Country, NUTS2_code, NUTS2_name, y1990...y2024]
    # Detect start year from the metadata block (search for 'Start year: <n>')
    start_year = 1990
    end_year = 2024
    try:
        # search metadata area
        meta = raw.iloc[:20, :]
        stacked = meta.astype(str).stack()
        sy = pd.to_numeric(stacked.str.extract(r"Start year:\s*(\d{4})", expand=False).dropna(), errors="coerce").dropna()
        ey = pd.to_numeric(stacked.str.extract(r"End year:\s*(\d{4})", expand=False).dropna(), errors="coerce").dropna()
        if not sy.empty:
            start_year = int(sy.iloc[0])
        if not ey.empty:
            end_year = int(ey.iloc[0])
    except Exception:
        pass
    years = list(range(start_year, end_year + 1))
    n_years = len(years)

    # Keep only Romania rows: either country code == 'ROU' or NUTS2 code starts with 'RO'
    rom_rows = rows[(rows.iloc[:, 1].astype(str).eq("ROU")) | (rows.iloc[:, 3].astype(str).str.startswith("RO"))].copy()
    if rom_rows.empty:
        print("[warn] No Romania NUTS2 rows found")
        return
    # Build tidy DataFrame
    rom_rows.columns = list(range(rom_rows.shape[1]))  # numeric cols
    meta_cols = rom_rows.iloc[:, :5].copy()
    values = rom_rows.iloc[:, 5:5 + n_years].copy()
    values.columns = years
    tidy = (
        values.assign(compound=meta_cols[0].values,
                      country_a3=meta_cols[1].values,
                      country=meta_cols[2].values,
                      nuts2_code=meta_cols[3].values,
                      nuts2_name=meta_cols[4].values)
        .melt(id_vars=["compound", "country_a3", "country", "nuts2_code", "nuts2_name"], var_name="year", value_name="value_gg")
    )
    # Convert Gg to kt (1 Gg = 1 kt); also provide Mt
    tidy["value_ktco2e"] = pd.to_numeric(tidy["value_gg"], errors="coerce")
    tidy["value_mtco2e"] = tidy["value_ktco2e"] / 1_000.0
    tidy.to_sql("edgar_nuts2_totals", conn, if_exists="replace", index=False)
    print(f"Stored {len(tidy)} rows to edgar_nuts2_totals")


## EDGAR-FOOD processing removed per user request.


def process_edgar_grid_romania(conn: sqlite3.Connection) -> None:
    if xr is None:
        print("[info] xarray not available; skipping NetCDF grid processing")
        return
    # File path for emissions (tonnes)
    nc_path = ROOT / "EDGAR/Annual_gridmaps/EDGAR_2025_GHG_CO2_2024_TOTALS_emi_nc/EDGAR_2025_GHG_CO2_2024_TOTALS_emi.nc"
    if not nc_path.exists():
        print("[info] Grid NetCDF not found; skipping grid processing")
        return
    ds = xr.open_dataset(nc_path)
    if "emissions" not in ds.data_vars:
        print("[warn] NetCDF 'emissions' var missing; skipped")
        return
    # Create a bounding box first to limit size; then clip to polygon if geopandas is available
    lat_min, lat_max = 43.4, 48.5
    lon_min, lon_max = 20.0, 30.1
    # Handle ascending/descending coordinate axes safely
    lat_vals = ds["lat"].values
    lon_vals = ds["lon"].values
    lat_slice = slice(lat_min, lat_max) if lat_vals[0] < lat_vals[-1] else slice(lat_max, lat_min)
    lon_slice = slice(lon_min, lon_max) if lon_vals[0] < lon_vals[-1] else slice(lon_max, lon_min)
    sub = ds.sel(lat=lat_slice, lon=lon_slice)

    da = sub["emissions"]  # tonnes
    df = da.to_dataframe(name="emissions_t").reset_index()

    # Clip to Romania polygon. Priority:
    # 1) If user provided a GeoJSON at data/romania_boundary.geojson, use it.
    # 2) Else, use a coarse built-in polygon approximation.
    if Polygon is not None and Point is not None:
        boundary_file = ROOT / "data/romania_boundary.geojson"
        try:
            rom_poly = None
            if boundary_file.exists():
                try:
                    if gpd is not None:
                        gdf = gpd.read_file(boundary_file)
                        geom = gdf.to_crs("EPSG:4326").geometry
                        # Prefer modern GeoPandas API; fallback to shapely if unavailable
                        if hasattr(geom, "union_all"):
                            rom_poly = geom.union_all()
                        else:
                            try:
                                from shapely.ops import unary_union as _unary_union  # type: ignore
                                rom_poly = _unary_union(list(geom))
                            except Exception:
                                # Last resort (may trigger deprecation warning in newer GeoPandas)
                                rom_poly = geom.unary_union
                    else:
                        print("[warn] GeoPandas not available to read GeoJSON; using coarse polygon")
                except Exception as e:
                    print(f"[warn] Failed reading custom boundary {boundary_file}: {e}")
            if rom_poly is None:
                # Coarse hand-crafted polygon around Romania borders (lon, lat)
                rom_pts = [
                    (20.26, 45.9), (21.0, 44.7), (22.7, 43.7), (25.1, 43.7), (27.0, 44.0),
                    (28.2, 45.1), (29.67, 45.42), (28.7, 46.0), (27.5, 47.3), (26.7, 47.9),
                    (25.0, 47.0), (22.8, 46.3), (21.5, 46.2), (20.26, 45.9)
                ]
                rom_poly = Polygon(rom_pts)

            mask = [rom_poly.contains(Point(lon, lat)) for lon, lat in zip(df["lon"], df["lat"])]
            df = df[pd.Series(mask, index=df.index)]
        except Exception as e:
            print(f"[warn] Could not apply Romania polygon clip: {e}")

    # Filter out zeros / NaNs
    df = df[pd.to_numeric(df["emissions_t"], errors="coerce").fillna(0) > 0].copy()
    # Add kt
    df["emissions_kt"] = df["emissions_t"] / 1_000.0
    # Round coordinates to reduce DB size
    df["lat"] = df["lat"].round(3)
    df["lon"] = df["lon"].round(3)
    df.to_sql("edgar_grid_emi_romania", conn, if_exists="replace", index=False)
    print(f"Stored {len(df)} rows to edgar_grid_emi_romania")


def main() -> None:
    conn = ensure_db()
    try:
        process_allcountries_summaries(conn)
        process_edgar_ucdb(conn)
        process_edgar_country_totals(conn)
        process_edgar_nuts2(conn)
        process_edgar_grid_romania(conn)
    finally:
        conn.close()
    print(f"Done. Database at: {DB_PATH}")


if __name__ == "__main__":
    main()
