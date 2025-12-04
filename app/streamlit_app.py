import os
from pathlib import Path

import pandas as pd
import streamlit as st
import sqlite3


ROOT = Path(__file__).resolve().parents[1]
DB_PATH = ROOT / "data/ro_ghg.db"


@st.cache_data(show_spinner=False)
def load_table(name: str) -> pd.DataFrame:
    if not DB_PATH.exists():
        return pd.DataFrame()
    con = sqlite3.connect(DB_PATH)
    try:
        return pd.read_sql_query(f"SELECT * FROM {name}", con)
    finally:
        con.close()


def kpi_card(label: str, value: float | int | str, suffix: str = ""):
    st.metric(label, f"{value}{suffix}")


def page_overview():
    st.subheader("Romania – Total GHG Emissions (EDGAR)")
    df = load_table("edgar_country_totals")
    if df.empty:
        st.info("No country totals in database. Run preprocessing.")
        return
    df = df[df["compound"].astype(str).str.contains("GWP_100_AR5_GHG", na=False)]
    ts = df.copy()
    ts["year"] = pd.to_numeric(ts["year"], errors="coerce")
    ts = ts.dropna(subset=["year"]).sort_values("year")
    if ts.empty:
        st.info("No data to display.")
        return
    st.line_chart(ts.set_index("year")["value_mtco2e"], y_label="MtCO₂e")
    # KPIs
    latest = ts.iloc[-1]
    first = ts.iloc[0]
    delta = latest["value_mtco2e"] - first["value_mtco2e"]
    c1, c2, c3 = st.columns(3)
    with c1:
        kpi_card("Latest year", int(latest["year"]))
    with c2:
        kpi_card("Emissions (latest)", round(float(latest["value_mtco2e"]), 1), " MtCO₂e")
    with c3:
        kpi_card("Change vs first", round(float(delta), 1), " MtCO₂e")


def page_regions():
    st.subheader("NUTS2 Regions – GHG Totals")
    df = load_table("edgar_nuts2_totals")
    if df.empty:
        st.info("No NUTS2 data in database. Run preprocessing.")
        return
    df = df[df["compound"].astype(str).str.contains("GWP_100_AR5_GHG", na=False)]
    # Sidebar controls
    years = sorted(pd.to_numeric(df["year"], errors="coerce").dropna().astype(int).unique().tolist())
    year = st.sidebar.slider("Select year", min_value=int(min(years)), max_value=int(max(years)), value=int(max(years)))
    # Ranking bar chart for selected year
    df["nuts2_name"] = df["nuts2_name"].fillna("(unknown)").astype(str)
    cur = df[df["year"].astype(int).eq(int(year))].sort_values("value_mtco2e", ascending=False)
    st.write(f"Regional totals in {year} (MtCO₂e)")
    st.bar_chart(cur.set_index("nuts2_name")["value_mtco2e"])

    # Timeseries for selected regions
    regions = sorted(set(df["nuts2_name"].dropna().astype(str).tolist()))
    sel = st.multiselect("Regions for timeseries", regions, default=regions[:3] if len(regions) >= 3 else regions)
    if sel:
        sub = df[df["nuts2_name"].isin(sel)].copy()
        sub["year"] = pd.to_numeric(sub["year"], errors="coerce")
        sub = sub.dropna(subset=["year"]).sort_values("year")
        ts = sub.pivot_table(index="year", columns="nuts2_name", values="value_mtco2e")
        st.line_chart(ts, y_label="MtCO₂e")


def page_cities():
    st.subheader("Cities – GHG Data")
    src = st.radio("Data source", ["Summary (totals)", "EDGAR UCDB (sectors)"])

    if src.startswith("Summary"):
        df = load_table("cities_only_summary")
        if df.empty:
            st.info("No cities summary in database. Run preprocessing.")
            return
        view = df[[c for c in ["Region Name", "Est. Population", "total_mt", "CO2 per capita"] if c in df.columns]].copy()
        view = view.rename(columns={"Region Name": "City", "total_mt": "Total (MtCO₂)"})
        st.dataframe(view.sort_values("Total (MtCO₂)", ascending=False), width='stretch')
        st.caption("Source: allcountries_onlycities_summary.xlsx")

    else:
        df = load_table("edgar_ucdb_cities_ghg")
        if df.empty:
            st.info("No UCDB (city-sector) data in database. Run preprocessing.")
            return
        df["UC_name"] = df["UC_name"].fillna("(unknown)").astype(str)
        df["sector"] = df["sector"].fillna("(unknown)").astype(str)
        cities = sorted(set(df["UC_name"].dropna().astype(str).tolist()))
        city = st.selectbox("City", cities)
        sectors = sorted(set(df["sector"].dropna().astype(str).tolist()))
        sel_sec = st.multiselect("Sectors", sectors, default=["Energy", "Transport", "Industry"])  # common sectors
        sub = df[(df["UC_name"].eq(city)) & (df["sector"].isin(sel_sec))]
        if sub.empty:
            st.info("No data for selection.")
            return
        pivot = sub.pivot_table(index="year", columns="sector", values="value_kt", aggfunc="sum").sort_index()
        st.area_chart(pivot, y_label="ktCO₂e")
        st.caption("Source: EDGAR_emiss_on_UCDB_2024.xlsx – GHG (AR5) sectors")




def page_grid_map():
    st.subheader("Spatial Grid – CO₂ Emissions (EDGAR 2024)")
    df = load_table("edgar_grid_emi_romania")
    if df.empty:
        st.info("No grid data in database. Run preprocessing.")
        return
    # Filter small values and control heatmap radius (in screen pixels)
    max_kt = float(df["emissions_kt"].max()) if not df.empty else 1.0
    min_kt = st.slider("Minimum cell emissions (kt)", min_value=0.0, max_value=max_kt, value=min(0.5, max_kt), step=max(0.1, max_kt/1000))
    radius_px = st.slider("Heatmap radius (pixels)", min_value=10, max_value=200, value=60, step=5)
    sub = df[df["emissions_kt"] >= min_kt].copy()

    import pydeck as pdk
    midpoint = [sub["lat"].mean(), sub["lon"].mean()] if not sub.empty else [45.8, 24.9]
    layer = pdk.Layer(
        "HeatmapLayer",
        data=sub,
        get_position="[lon, lat]",
        get_weight="emissions_kt",
        aggregation="SUM",
        radiusPixels=radius_px,
    )
    view_state = pdk.ViewState(latitude=midpoint[0], longitude=midpoint[1], zoom=6.2, pitch=0)
    st.pydeck_chart(pdk.Deck(layers=[layer], initial_view_state=view_state))
    st.caption("Source: EDGAR_2025_GHG_CO2_2024_TOTALS_emi.nc – kt per grid cell; 2D heatmap")


def main():
    st.set_page_config(page_title="Romania GHG Dashboard", layout="wide")
    st.title("Romania GHG Dashboard")
    if not DB_PATH.exists():
        st.warning("Database not found. Run preprocessing script first: `env/bin/python scripts/preprocess.py`")

    tabs = st.tabs(["Overview", "Regions", "Cities", "Grid Map"])  # noqa
    with tabs[0]:
        page_overview()
    with tabs[1]:
        page_regions()
    with tabs[2]:
        page_cities()
    with tabs[3]:
        page_grid_map()


if __name__ == "__main__":
    main()
