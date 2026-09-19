#!/usr/bin/env python3
"""Pre-compute Transmission & Spread context tables for the dashboard.

Reads epi and genomic data from the DuckDB knowledge warehouse and writes:
  - transmission_burden.tsv
  - transmission_spatial.tsv
  - transmission_anomaly.tsv
  - transmission_potential.tsv
  - transmission_outbreak_summary.tsv
  - strain_transmission_profile.tsv

Requires only python + pandas + numpy + duckdb (pgirl_knowledge env).
"""
import argparse
import contextlib
import datetime
import math
import re
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd


def _to_year_month(value) -> str | None:
    """Return 'YYYY-MM' from a datetime-like or decimal-year value."""
    if pd.isna(value):
        return None
    if isinstance(value, (pd.Timestamp, datetime.datetime, datetime.date)):
        return f"{value.year}-{value.month:02d}"
    # decimal year, e.g. 2014.78
    if isinstance(value, (int, float, np.integer, np.floating)):
        year = int(value)
        month = int((value - year) * 12) + 1
        month = max(1, min(12, month))
        return f"{year}-{month:02d}"
    # string/ISO date fallback
    with contextlib.suppress(Exception):
        dt = pd.to_datetime(value)
        return f"{dt.year}-{dt.month:02d}"
    return None


def _to_timestamp(value) -> pd.Timestamp | None:
    """Return a pandas Timestamp from a datetime-like or decimal-year value."""
    if pd.isna(value):
        return None
    if isinstance(value, (pd.Timestamp, datetime.datetime, datetime.date)):
        return pd.Timestamp(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        year = int(value)
        month = int((value - year) * 12) + 1
        month = max(1, min(12, month))
        return pd.Timestamp(f"{year}-{month:02d}-15")
    with contextlib.suppress(Exception):
        return pd.to_datetime(value)
    return None


EPIDEMIOLOGICAL_QUERY = """
SELECT
    r.record_date,
    r.report_date,
    r.reference_date,
    CAST(COALESCE(r.country, gl1.country, gl2.country) AS TEXT) AS country,
    CAST(COALESCE(r.admin1, gl1.admin1, gl2.admin1) AS TEXT) AS admin1,
    CAST(COALESCE(r.admin2, gl1.admin2, gl2.admin2) AS TEXT) AS admin2,
    r.location_id,
    r.cases,
    r.deaths,
    r.suspected,
    r.recovered,
    r.value,
    r.measure,
    r.indicator_label,
    r.case_classification,
    d.dataset_id,
    d.dataset_name,
    CAST(COALESCE(gl1.latitude, gl2.latitude) AS DOUBLE) AS latitude,
    CAST(COALESCE(gl1.longitude, gl2.longitude) AS DOUBLE) AS longitude
FROM epidemiological_records r
JOIN epidemiological_datasets d ON r.dataset_id = d.dataset_id
LEFT JOIN geographic_locations gl1 ON CAST(r.location_id AS BIGINT) = gl1.location_id
LEFT JOIN geographic_locations gl2 ON r.location_code = gl2.location_code
    AND r.location_code_type = gl2.location_code_type
WHERE d.species = ? AND d.dataset_type != 'summary'
"""

SAMPLE_QUERY = """
SELECT
    s.sample_id,
    s.sample_name,
    s.ppx_accession,
    s.insdc_accession,
    s.collection_date,
    s.country,
    s.admin1,
    s.locality,
    COALESCE(s.clade, t.clade) AS clade,
    s.lineage,
    s.outbreak,
    s.host,
    s.is_query,
    gl.latitude,
    gl.longitude,
    t.tip_id,
    t.tip_date,
    t.div,
    t.aa_mutation_count,
    t.nuc_mutation_count
FROM samples s
LEFT JOIN sample_geo_location sgl ON s.sample_id = sgl.sample_id
LEFT JOIN geographic_locations gl ON sgl.location_id = gl.location_id
LEFT JOIN tree_tips t ON s.sample_id = t.sample_id
WHERE s.species = ?
"""


def classify_indicator(text: str) -> str:
    """Classify a measure/indicator string as case, death, suspected, recovered, or other."""
    if not text:
        return "other"
    t = text.lower()
    if "death" in t or "killed" in t or "mortality" in t or "deced" in t:
        return "death"
    if "suspect" in t and "case" in t:
        return "suspected"
    if "recover" in t or "cured" in t:
        return "recovered"
    if "case" in t or "infected" in t or "affected" in t or "human cases" in t:
        return "case"
    return "other"


def to_float(v):
    if v is None:
        return np.nan
    try:
        return float(v)
    except (ValueError, TypeError):
        return np.nan


def parse_epi_records(df: pd.DataFrame) -> pd.DataFrame:
    """Harmonise raw epi records into one daily value per dataset/location."""
    if df.empty:
        return pd.DataFrame(
            columns=[
                "record_date",
                "country",
                "admin1",
                "location_id",
                "dataset_id",
                "latitude",
                "longitude",
                "daily_cases",
                "daily_deaths",
                "daily_suspected",
                "daily_recovered",
            ]
        )

    df = df.copy()
    date_col = df["record_date"].fillna(df["report_date"]).fillna(df["reference_date"])
    df["record_date"] = pd.to_datetime(date_col, errors="coerce")
    df = df.dropna(subset=["record_date"])

    # value from long/pivoted formats
    df["value_f"] = df["value"].apply(to_float)
    df["cases_f"] = df["cases"].apply(to_float)
    df["deaths_f"] = df["deaths"].apply(to_float)
    df["suspected_f"] = df["suspected"].apply(to_float)
    df["recovered_f"] = df["recovered"].apply(to_float)

    # classify each row
    df["itype"] = (df["measure"].fillna("") + " " + df["indicator_label"].fillna("")).apply(
        classify_indicator
    )

    # Fill in missing cases/deaths from value when the indicator matches
    m_case = df["itype"] == "case"
    m_death = df["itype"] == "death"
    m_susp = df["itype"] == "suspected"
    m_recv = df["itype"] == "recovered"

    df.loc[m_case, "cases_f"] = df.loc[m_case, "cases_f"].fillna(df.loc[m_case, "value_f"])
    df.loc[m_death, "deaths_f"] = df.loc[m_death, "deaths_f"].fillna(df.loc[m_death, "value_f"])
    df.loc[m_susp, "suspected_f"] = df.loc[m_susp, "suspected_f"].fillna(df.loc[m_susp, "value_f"])
    df.loc[m_recv, "recovered_f"] = df.loc[m_recv, "recovered_f"].fillna(df.loc[m_recv, "value_f"])

    # Negative reported counts are almost always data corrections; clamp to 0
    for c in ["cases_f", "deaths_f", "suspected_f", "recovered_f"]:
        df[c] = df[c].clip(lower=0)

    # Choose the broadest case/death value per (dataset, country, admin1, date)
    # by taking the maximum among same-day rows of the same type.
    admin_fill = df["admin1"].fillna("National")
    df["admin1"] = admin_fill
    df["country"] = df["country"].fillna("Unknown")

    grouped = (
        df.groupby(["dataset_id", "country", "admin1", "record_date"], as_index=False)
        .agg(
            daily_cases=("cases_f", "max"),
            daily_deaths=("deaths_f", "max"),
            daily_suspected=("suspected_f", "max"),
            daily_recovered=("recovered_f", "max"),
            latitude=("latitude", "first"),
            longitude=("longitude", "first"),
        )
    )
    return grouped


COUNTRY_NAME_MAP = {
    "#adm0 +code": None,
    "COD": "Democratic Republic of the Congo",
    "Zaire": "Democratic Republic of the Congo",
    "Zaire (Democratic Republic of the Congo - DRC)": "Democratic Republic of the Congo",
    "Democratic Republic of the Congo (formerly Zaire)": "Democratic Republic of the Congo",
    "Sudan (South Sudan)": "South Sudan",
    "Guinea 2": "Guinea",
    "Guinée": "Guinea",
    "Liberia 2": "Liberia",
    "England": "United Kingdom",
}


def _fill_country_coords(samples: pd.DataFrame, con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Backfill missing sample lat/lon from the country_coordinates table."""
    try:
        cc = con.execute(
            "SELECT country, latitude, longitude FROM country_coordinates"
        ).df()
    except (duckdb.CatalogException, duckdb.InvalidInputException):
        return samples
    if cc.empty:
        return samples
    samples = samples.merge(cc, on="country", how="left", suffixes=("", "_cc"))
    samples["latitude"] = samples["latitude"].astype(float).fillna(samples["latitude_cc"].astype(float))
    samples["longitude"] = samples["longitude"].astype(float).fillna(samples["longitude_cc"].astype(float))
    samples = samples.drop(columns=["latitude_cc", "longitude_cc"], errors="ignore")
    return samples


def _haversine_km(lat1, lon1, lat2, lon2):
    """Vectorised great-circle distance in km between two lat/lon arrays."""
    r = 6371.0
    lat1 = pd.to_numeric(lat1, errors="coerce").astype(float)
    lat2 = pd.to_numeric(lat2, errors="coerce").astype(float)
    lon1 = pd.to_numeric(lon1, errors="coerce").astype(float)
    lon2 = pd.to_numeric(lon2, errors="coerce").astype(float)
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = np.radians(lat2 - lat1)
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2) ** 2
    return 2 * r * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def build_daily_burden(records: pd.DataFrame) -> pd.DataFrame:
    """Combine datasets to one daily count per country/admin1."""
    if records.empty:
        return records
    records["country"] = records["country"].replace(COUNTRY_NAME_MAP).fillna("Unknown")
    # If both national and sub-national datasets exist, prefer sub-national rows
    # by taking the *maximum* daily value per (country, admin1, date).
    # This avoids double-counting while keeping the broadest reported number.
    daily = (
        records.groupby(["country", "admin1", "record_date"], as_index=False)
        .agg(
            {
                "daily_cases": "max",
                "daily_deaths": "max",
                "daily_suspected": "max",
                "daily_recovered": "max",
                "latitude": "first",
                "longitude": "first",
            }
        )
    )
    daily = daily.sort_values(["country", "admin1", "record_date"]).reset_index(drop=True)
    daily["daily_cases"] = daily["daily_cases"].fillna(0)
    daily["daily_deaths"] = daily["daily_deaths"].fillna(0)
    return daily


def add_new_counts(df: pd.DataFrame) -> pd.DataFrame:
    """Derive new cases/deaths from (assumed) cumulative daily counts."""
    if df.empty:
        return df
    df = df.sort_values(["country", "admin1", "record_date"]).copy()
    df["new_cases"] = df.groupby(["country", "admin1"])["daily_cases"].diff().fillna(df["daily_cases"])
    df["new_deaths"] = df.groupby(["country", "admin1"])["daily_deaths"].diff().fillna(df["daily_deaths"])
    # Corrections can make differences negative
    df["new_cases"] = df["new_cases"].clip(lower=0)
    df["new_deaths"] = df["new_deaths"].clip(lower=0)
    df["cfr_cum"] = df["daily_deaths"] / (df["daily_cases"] + 1e-9)
    df["cfr_new"] = df["new_deaths"] / (df["new_cases"] + 1e-9)
    return df


def weekly_aggregation(df: pd.DataFrame) -> pd.DataFrame:
    """Resample daily counts to ISO weeks (Monday-based)."""
    if df.empty:
        return df
    weekly = (
        df.groupby(["country", "admin1", pd.Grouper(key="record_date", freq="W-MON")])
        .agg(
            cases_cum=("daily_cases", "last"),
            deaths_cum=("daily_deaths", "last"),
            daily_suspected=("daily_suspected", "sum"),
            daily_recovered=("daily_recovered", "sum"),
            new_cases=("new_cases", "sum"),
            new_deaths=("new_deaths", "sum"),
            latitude=("latitude", "first"),
            longitude=("longitude", "first"),
        )
        .reset_index()
    )
    weekly = weekly.rename(columns={"record_date": "week_start"})
    weekly["cfr_cum"] = weekly["deaths_cum"] / (weekly["cases_cum"] + 1e-9)
    weekly["cfr_new"] = weekly["new_deaths"] / (weekly["new_cases"] + 1e-9)
    weekly["week_start"] = weekly["week_start"].dt.date
    return weekly


def mann_kendall(y: np.ndarray) -> tuple:
    """Return (p_value, slope) for Mann-Kendall trend + Theil-Sen slope."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3 or np.all(np.isnan(y)):
        return 1.0, 0.0
    # Theil-Sen slope (median of all pair-wise slopes)
    slopes = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            if not (np.isnan(y[i]) or np.isnan(y[j])) and (j - i) != 0:
                slopes.append((y[j] - y[i]) / (j - i))
    sen_slope = float(np.median(slopes)) if slopes else 0.0

    # Mann-Kendall S
    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            if np.isnan(y[i]) or np.isnan(y[j]):
                continue
            if y[j] > y[i]:
                s += 1
            elif y[j] < y[i]:
                s -= 1
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    if var_s == 0 or s == 0:
        p = 1.0
    else:
        z = (s - 1) / math.sqrt(var_s) if s > 0 else (s + 1) / math.sqrt(var_s)
        # two-sided p-value from erfc
        p = math.erfc(abs(z) / math.sqrt(2))
    return p, sen_slope


def wilson_ci(k: pd.Series, n: pd.Series, z: float = 1.96) -> pd.DataFrame:
    """Wilson score interval for a binomial proportion."""
    n_safe = np.maximum(0, n) + 1e-9
    p = np.clip(k / n_safe, 0, 1)
    denom = 1 + z * z / n_safe
    centre = (p + z * z / (2 * n_safe)) / denom
    half = z * np.sqrt(np.maximum(0, (p * (1 - p) + z * z / (4 * n_safe)) / n_safe)) / denom
    return pd.DataFrame({"cfr_lower": centre - half, "cfr_upper": centre + half})


def add_trend_stats(weekly: pd.DataFrame) -> pd.DataFrame:
    """Add Mann-Kendall / Sen slope and smoothing for each location."""
    if weekly.empty:
        return weekly
    weekly = weekly.sort_values(["country", "admin1", "week_start"]).copy()

    stats = []
    for _, sub in weekly.groupby(["country", "admin1"]):
        vals = sub["new_cases"].values
        p, slope = mann_kendall(vals)
        p_d, slope_d = mann_kendall(sub["new_deaths"].values)
        p_c, slope_c = mann_kendall(sub["cfr_cum"].fillna(0).values)
        sub = sub.copy()
        sub["mk_p_cases"] = p
        sub["sen_slope_cases"] = slope
        sub["mk_p_deaths"] = p_d
        sub["sen_slope_deaths"] = slope_d
        sub["mk_p_cfr"] = p_c
        sub["sen_slope_cfr"] = slope_c
        sub["trend_cases"] = np.where(p < 0.05, np.sign(slope), 0)
        stats.append(sub)
    weekly = pd.concat(stats, ignore_index=True)

    weekly["ma4_new_cases"] = (
        weekly.groupby(["country", "admin1"])["new_cases"]
        .transform(lambda s: s.rolling(window=4, min_periods=1).mean())
    )
    weekly["ma4_new_deaths"] = (
        weekly.groupby(["country", "admin1"])["new_deaths"]
        .transform(lambda s: s.rolling(window=4, min_periods=1).mean())
    )
    weekly["growth_rate"] = (
        weekly.groupby(["country", "admin1"])["new_cases"]
        .transform(lambda s: (s - s.shift(1)) / (s.shift(1).replace(0, np.nan).fillna(s)))
    )

    ci = wilson_ci(weekly["new_deaths"], weekly["new_cases"])
    weekly = pd.concat([weekly, ci], axis=1)
    return weekly


def add_anomaly_stats(weekly: pd.DataFrame) -> pd.DataFrame:
    """Add EWMA and CUSUM alert columns."""
    if weekly.empty:
        return weekly
    weekly = weekly.sort_values(["country", "admin1", "week_start"]).copy()

    frames = []
    for _, sub in weekly.groupby(["country", "admin1"]):
        x = sub["new_cases"].fillna(0).astype(float).values
        sub = sub.copy()
        if len(x) < 2:
            sub["ewma"] = np.nan
            sub["ewma_ucl"] = np.nan
            sub["cusum_pos"] = np.nan
            sub["cusum_neg"] = np.nan
            sub["alert"] = False
            frames.append(sub)
            continue
        # historical baseline = first 4 weeks or all if short
        n_hist = min(4, len(x) - 1)
        mu = float(np.mean(x[:n_hist])) if n_hist > 0 else float(np.mean(x))
        sd = float(np.std(x[:n_hist])) if n_hist > 0 else float(np.std(x))
        if sd == 0 or np.isnan(sd):
            sd = 1.0

        lam = 0.3
        ewma = np.empty_like(x, dtype=float)
        ewma[0] = mu
        for i in range(1, len(x)):
            ewma[i] = lam * x[i] + (1 - lam) * ewma[i - 1]

        c_pos = np.zeros_like(x, dtype=float)
        c_neg = np.zeros_like(x, dtype=float)
        k = 0.5 * sd
        for i in range(1, len(x)):
            c_pos[i] = max(0.0, c_pos[i - 1] + x[i] - (mu + k))
            c_neg[i] = min(0.0, c_neg[i - 1] + x[i] - (mu - k))

        ucl = mu + 2.5 * sd
        sub["ewma"] = ewma
        sub["ewma_ucl"] = ucl
        sub["cusum_pos"] = c_pos
        sub["cusum_neg"] = c_neg
        sub["alert"] = (sub["new_cases"] > ucl) | (c_pos > 4 * sd) | (c_neg < -4 * sd)
        frames.append(sub)
    weekly = pd.concat(frames, ignore_index=True)
    return weekly


def build_spatial(weekly: pd.DataFrame, con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Aggregate latest cumulative counts per admin1 and compute Getis-Ord Gi* z-scores."""
    if weekly.empty:
        return pd.DataFrame(
            columns=[
                "country", "admin1", "latitude", "longitude", "cases",
                "deaths", "cfr", "g_z", "g_p"
            ]
        )
    latest = (
        weekly.sort_values("week_start")
        .groupby(["country", "admin1"], as_index=False)
        .last()[["country", "admin1", "cases_cum", "deaths_cum"]]
        .rename(columns={"cases_cum": "cases", "deaths_cum": "deaths"})
    )
    latest["cfr"] = latest["deaths"] / (latest["cases"] + 1e-9)

    # Use country-level centroids (more complete than admin1 coords)
    try:
        geo = con.execute(
            "SELECT country, AVG(latitude) AS latitude, AVG(longitude) AS longitude "
            "FROM geographic_locations "
            "WHERE country IS NOT NULL AND latitude IS NOT NULL AND longitude IS NOT NULL "
            "GROUP BY country"
        ).df()
        if not geo.empty:
            latest = latest.merge(geo, on="country", how="left")
    except Exception:
        latest["latitude"] = np.nan
        latest["longitude"] = np.nan

    for col in ["latitude", "longitude"]:
        if col not in latest.columns:
            latest[col] = np.nan

    pts = latest.dropna(subset=["latitude", "longitude"]).copy()
    if len(pts) < 2:
        latest["g_z"] = np.nan
        latest["g_p"] = np.nan
        return latest

    coords = pts[["latitude", "longitude"]].to_numpy()
    x = pts["cases"].to_numpy(dtype=float)
    n = len(x)

    # Euclidean distances on lat/lon (sufficient for relative ordering and small areas)
    dist = np.sqrt(((coords[:, None, :] - coords[None, :, :]) ** 2).sum(axis=2))
    # k-nearest neighbours (including self)
    k = min(5, n)
    w = np.zeros_like(dist, dtype=float)
    for i in range(n):
        idx = np.argpartition(dist[i], k - 1)[:k]
        w[i, idx] = 1.0

    # local G* (point included through w_ii = 1)
    total = x.sum()
    g = (w @ x) / (total + 1e-9)
    W = w.sum(axis=1)
    p = W / n
    mean_x = x.mean()
    var_x = np.maximum(0.0, (x ** 2).mean() - mean_x ** 2)
    if var_x <= 0:
        pts["g_z"] = 0.0
        pts["g_p"] = 1.0
    else:
        # variance of the local sum under randomisation
        var = np.maximum(0.0, (W * (n - W) / (n * (n - 1))) * var_x / (total ** 2))
        z = (g - p) / np.sqrt(var + 1e-9)
        pts["g_z"] = z
        pts["g_p"] = np.where(z > 0, 0.5 * np.exp(-(z ** 2) / 2), 1.0)

    # merge back to all locations
    latest = latest.merge(pts[["country", "admin1", "g_z", "g_p"]], on=["country", "admin1"], how="left")
    return latest


def build_anomaly(weekly: pd.DataFrame) -> pd.DataFrame:
    """Return only flagged anomalous weeks."""
    if weekly.empty or "alert" not in weekly.columns:
        return pd.DataFrame()
    flagged = weekly[weekly["alert"]].copy()
    cols = [
        "country",
        "admin1",
        "week_start",
        "new_cases",
        "new_deaths",
        "ewma",
        "ewma_ucl",
        "cusum_pos",
        "cusum_neg",
    ]
    return flagged[[c for c in cols if c in flagged.columns]]


def build_potential(
    con: duckdb.DuckDBPyConnection,
    species: str,
    weekly: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build a raw query-to-closest-genome transmission potential table.

    Each query is linked to its N closest historical background genomes by
    phylogenetic distance (div).  The epidemiological burden of each
    background's country/admin1 and year is looked up and returned for
    downstream transmission potential/burden analysis.
    """
    samples = con.execute(SAMPLE_QUERY, [species]).df()
    if samples.empty:
        return samples
    samples["collection_date"] = pd.to_datetime(samples["collection_date"], errors="coerce")
    samples["tip_date"] = pd.to_datetime(samples["tip_date"].apply(_to_timestamp), errors="coerce")
    samples["country"] = (
        samples["country"]
        .fillna("Unknown")
        .astype(str)
        .replace(["<NA>", "nan", "None"], "Unknown")
        .replace(COUNTRY_NAME_MAP)
        .fillna("Unknown")
    )
    samples["admin1"] = (
        samples["admin1"]
        .astype(str)
        .replace(["<NA>", "nan", "None", ""], "National")
        .fillna("National")
    )
    samples["clade"] = samples["clade"].astype(str).replace(["<NA>", "nan", "None", "unassigned"], np.nan)
    samples["lineage"] = samples["lineage"].astype(str).replace(["<NA>", "nan", "None"], np.nan)
    samples["outbreak"] = samples["outbreak"].astype(str).replace(["<NA>", "nan", "None", "unassigned"], np.nan)
    samples["strain"] = (
        samples["clade"]
        .combine_first(samples["lineage"])
        .combine_first(samples["outbreak"])
    )
    samples["year"] = np.where(
        samples["collection_date"].notna(),
        samples["collection_date"].dt.year,
        samples["tip_date"].dt.year,
    )

    samples = _fill_country_coords(samples, con)

    queries = samples[samples["is_query"].astype(str).str.lower() == "true"].copy()
    if queries.empty:
        return pd.DataFrame()

    background = samples[samples["is_query"].astype(str).str.lower() != "true"].copy()

    epi_admin, epi_country = build_yearly_burden(weekly if weekly is not None else pd.DataFrame())

    def _lookup_burden(country, admin1, year):
        """Look up annual epi burden with admin1 -> country -> nearest-year fallbacks."""
        result = {
            "epi_cases": np.nan,
            "epi_deaths": np.nan,
            "epi_cfr": np.nan,
            "epi_year": np.nan,
            "epi_available": False,
            "epi_match_level": "none",
        }
        if pd.isna(year) or country in (None, "Unknown", "nan", "None"):
            return result

        # 1. Exact admin1 + year match
        if not epi_admin.empty:
            row = epi_admin[
                (epi_admin["country"] == country)
                & (epi_admin["admin1"] == admin1)
                & (epi_admin["year"] == year)
            ]
            if not row.empty:
                r = row.iloc[0]
                result.update(
                    epi_cases=r["epi_cases"],
                    epi_deaths=r["epi_deaths"],
                    epi_cfr=r["epi_cfr"],
                    epi_year=year,
                    epi_available=True,
                    epi_match_level="admin1",
                )
                return result

        if epi_country.empty:
            return result

        # 2. Country + year match
        row = epi_country[
            (epi_country["country"] == country) & (epi_country["year"] == year)
        ]
        if not row.empty:
            r = row.iloc[0]
            result.update(
                epi_cases=r["epi_cases"],
                epi_deaths=r["epi_deaths"],
                epi_cfr=r["epi_cfr"],
                epi_year=year,
                epi_available=True,
                epi_match_level="country",
            )
            return result

        # 3. Nearest available year (within +/-2) at country level
        cand = epi_country[epi_country["country"] == country].copy()
        if not cand.empty:
            cand["year_diff"] = (cand["year"] - year).abs()
            cand = cand[cand["year_diff"] <= 2].sort_values("year_diff")
            if not cand.empty:
                r = cand.iloc[0]
                result.update(
                    epi_cases=r["epi_cases"],
                    epi_deaths=r["epi_deaths"],
                    epi_cfr=r["epi_cfr"],
                    epi_year=int(r["year"]),
                    epi_available=True,
                    epi_match_level="country_nearest_year",
                )
        return result

    rows = []
    for _, q in queries.iterrows():
        q_name = q["sample_name"]
        q_date = q["collection_date"] if not pd.isna(q["collection_date"]) else q["tip_date"]
        q_year = q["year"] if not pd.isna(q["year"]) else None
        q_div = q["div"]

        # No phylogenetic distance => cannot pick close relatives; still emit a query row.
        if background.empty or pd.isna(q_div):
            rows.append(
                {
                    "query_sample": q_name,
                    "query_collection_date": q_date.date() if not pd.isna(q_date) else None,
                    "query_country": q["country"],
                    "query_admin1": q["admin1"],
                    "query_clade": q["clade"],
                    "query_lineage": q["lineage"],
                    "query_outbreak": q["outbreak"],
                    "query_strain": q["strain"],
                    "query_ppx_accession": q["ppx_accession"],
                    "query_div": q_div,
                    "query_latitude": q["latitude"],
                    "query_longitude": q["longitude"],
                    "background_sample": np.nan,
                    "background_div": np.nan,
                    "background_div_diff": np.nan,
                    "background_collection_date": None,
                    "background_country": np.nan,
                    "background_admin1": np.nan,
                    "background_year": np.nan,
                    "background_latitude": np.nan,
                    "background_longitude": np.nan,
                    "background_clade": np.nan,
                    "background_lineage": np.nan,
                    "background_strain": np.nan,
                    "background_epi_year": np.nan,
                    "background_annual_cases": np.nan,
                    "background_annual_deaths": np.nan,
                    "background_annual_cfr": np.nan,
                    "background_epi_available": False,
                    "background_epi_match_level": "none",
                    "potential_label": "No background or divergence information",
                }
            )
            continue

        bg = background.copy()
        bg["div_diff"] = (bg["div"] - q_div).abs()
        closest = bg.sort_values("div_diff").head(10)

        for _, nb in closest.iterrows():
            bg_date = nb["collection_date"] if not pd.isna(nb["collection_date"]) else nb["tip_date"]
            bg_year = nb["year"] if not pd.isna(nb["year"]) else q_year
            if pd.isna(bg_year):
                continue
            b = _lookup_burden(nb["country"], nb["admin1"], int(bg_year))
            rows.append(
                {
                    "query_sample": q_name,
                    "query_collection_date": q_date.date() if not pd.isna(q_date) else None,
                    "query_country": q["country"],
                    "query_admin1": q["admin1"],
                    "query_clade": q["clade"],
                    "query_lineage": q["lineage"],
                    "query_outbreak": q["outbreak"],
                    "query_strain": q["strain"],
                    "query_ppx_accession": q["ppx_accession"],
                    "query_div": q_div,
                    "query_latitude": q["latitude"],
                    "query_longitude": q["longitude"],
                    "background_sample": nb["sample_name"],
                    "background_div": nb["div"],
                    "background_div_diff": nb["div_diff"],
                    "background_collection_date": bg_date.date() if not pd.isna(bg_date) else None,
                    "background_country": nb["country"],
                    "background_admin1": nb["admin1"],
                    "background_year": int(bg_year),
                    "background_latitude": nb["latitude"],
                    "background_longitude": nb["longitude"],
                    "background_clade": nb["clade"],
                    "background_lineage": nb["lineage"],
                    "background_strain": nb["strain"],
                    "background_epi_year": b["epi_year"],
                    "background_annual_cases": float(b["epi_cases"]) if not pd.isna(b["epi_cases"]) else np.nan,
                    "background_annual_deaths": float(b["epi_deaths"]) if not pd.isna(b["epi_deaths"]) else np.nan,
                    "background_annual_cfr": float(b["epi_cfr"]) if not pd.isna(b["epi_cfr"]) else np.nan,
                    "background_epi_available": bool(b["epi_available"]),
                    "background_epi_match_level": b["epi_match_level"],
                    "potential_label": (
                        f"Closest historical genome {nb['sample_name']} "
                        f"(div diff {nb['div_diff']:.4f})"
                    ),
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out["geo_distance_km"] = _haversine_km(
            out["query_latitude"],
            out["query_longitude"],
            out["background_latitude"],
            out["background_longitude"],
        )
        qd = pd.to_datetime(out["query_collection_date"], errors="coerce")
        bd = pd.to_datetime(out["background_collection_date"], errors="coerce")
        out["temporal_distance_days"] = (qd - bd).abs().dt.days
    return out


def build_strain_profile(
    con: duckdb.DuckDBPyConnection,
    species: str,
    weekly: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """One row per clade summarising historical transmission behaviour.

    Answers "how has this strain behaved?": temporal span, geographic spread,
    and the epidemiological burden observed in the place-years where the
    strain circulated (de-duplicated across co-circulating clades).
    """
    cols = [
        "strain",
        "strain_level",
        "lineages",
        "n_genomes",
        "first_seen",
        "last_seen",
        "active_years",
        "n_countries",
        "countries",
        "n_admin1",
        "origin_country",
        "max_spread_km",
        "spread_rate_km_per_year",
        "linked_cases",
        "linked_deaths",
        "mean_cfr",
        "behavior_label",
    ]
    samples = con.execute(SAMPLE_QUERY, [species]).df()
    if samples.empty:
        return pd.DataFrame(columns=cols)

    samples["collection_date"] = pd.to_datetime(samples["collection_date"], errors="coerce")
    samples["tip_date"] = pd.to_datetime(samples["tip_date"].apply(_to_timestamp), errors="coerce")
    samples["country"] = (
        samples["country"]
        .fillna("Unknown")
        .astype(str)
        .replace(["<NA>", "nan", "None"], "Unknown")
        .replace(COUNTRY_NAME_MAP)
        .fillna("Unknown")
    )
    samples["admin1"] = (
        samples["admin1"]
        .astype(str)
        .replace(["<NA>", "nan", "None", ""], "National")
        .fillna("National")
    )
    samples["clade"] = samples["clade"].astype(str).replace(
        ["<NA>", "nan", "None", "unassigned"], np.nan
    )
    samples["lineage"] = samples["lineage"].astype(str).replace(
        ["<NA>", "nan", "None"], np.nan
    )
    samples["outbreak"] = samples["outbreak"].astype(str).replace(
        ["<NA>", "nan", "None", "unassigned"], np.nan
    )
    # Strain = first available of clade -> lineage -> outbreak
    samples["strain"] = (
        samples["clade"]
        .combine_first(samples["lineage"])
        .combine_first(samples["outbreak"])
    )
    samples["strain_level"] = np.select(
        [
            samples["clade"].notna(),
            samples["lineage"].notna(),
            samples["outbreak"].notna(),
        ],
        ["clade", "lineage", "outbreak"],
        default="none",
    )
    samples["sample_date"] = samples["collection_date"].where(
        samples["collection_date"].notna(), samples["tip_date"]
    )
    samples["year"] = samples["sample_date"].dt.year
    samples = _fill_country_coords(samples, con)

    bg = samples[
        (samples["is_query"].astype(str).str.lower() != "true")
        & samples["strain"].notna()
    ].copy()
    if bg.empty:
        return pd.DataFrame(columns=cols)

    epi_admin, epi_country = build_yearly_burden(
        weekly if weekly is not None else pd.DataFrame()
    )

    # Distinct place-years where each strain was observed
    obs = (
        bg[["strain", "country", "admin1", "year"]]
        .dropna(subset=["year"])
        .drop_duplicates()
    )
    obs["year"] = obs["year"].astype(int)
    # Number of distinct strains per place-year (for burden de-duplication)
    obs["n_strains_place_year"] = obs.groupby(["country", "admin1", "year"])[
        "strain"
    ].transform("nunique")

    admin_lookup = {}
    if not epi_admin.empty:
        for _, r in epi_admin.iterrows():
            admin_lookup[(r["country"], r["admin1"], int(r["year"]))] = r
    country_lookup = {}
    if not epi_country.empty:
        for _, r in epi_country.iterrows():
            country_lookup[(r["country"], int(r["year"]))] = r

    def _burden(country, admin1, year):
        r = admin_lookup.get((country, admin1, year))
        if r is None:
            r = country_lookup.get((country, year))
        if r is None:
            return pd.Series(
                {"epi_cases": np.nan, "epi_deaths": np.nan, "epi_cfr": np.nan}
            )
        return r[["epi_cases", "epi_deaths", "epi_cfr"]]

    rows = []
    for strain, g in bg.groupby("strain"):
        dates = g["sample_date"].dropna()
        first_seen = dates.min() if not dates.empty else pd.NaT
        last_seen = dates.max() if not dates.empty else pd.NaT
        active_years = int(g["year"].nunique()) if g["year"].notna().any() else 0
        countries = sorted(g["country"].dropna().unique())
        n_admin1 = int(g.loc[g["admin1"] != "National", "admin1"].nunique())

        origin = np.nan
        if not dates.empty:
            origin = g.loc[g["sample_date"] == first_seen, "country"].iloc[0]

        # Max pairwise distance between distinct sampled locations
        pts = g[["latitude", "longitude"]].dropna().drop_duplicates().values
        max_spread = np.nan
        if len(pts) >= 2:
            i, j = np.triu_indices(len(pts), k=1)
            d = _haversine_km(pts[i, 0], pts[i, 1], pts[j, 0], pts[j, 1])
            max_spread = float(np.nanmax(d)) if len(d) else np.nan
        elif len(pts) == 1:
            max_spread = 0.0

        span_years = (
            (last_seen - first_seen).days / 365.25
            if not dates.empty and len(dates) > 1
            else 0.0
        )
        spread_rate = (
            max_spread / span_years
            if not pd.isna(max_spread) and span_years > 0
            else np.nan
        )

        # Linked burden: each place-year's burden shared across co-circulating strains
        o = obs[obs["strain"] == strain]
        linked_cases = 0.0
        linked_deaths = 0.0
        cfrs = []
        for _, r in o.iterrows():
            b = _burden(r["country"], r["admin1"], int(r["year"]))
            share = 1.0 / r["n_strains_place_year"]
            if not pd.isna(b["epi_cases"]):
                linked_cases += b["epi_cases"] * share
            if not pd.isna(b["epi_deaths"]):
                linked_deaths += b["epi_deaths"] * share
            if not pd.isna(b["epi_cfr"]):
                cfrs.append(b["epi_cfr"])
        mean_cfr = float(np.mean(cfrs)) if cfrs else np.nan

        n_countries = len(countries)
        n_genomes = len(g)
        if n_genomes < 2:
            label = "sporadic"
        elif n_countries >= 3 and span_years >= 2:
            label = "persistent_multi_country"
        elif n_countries >= 2:
            label = "multi_country_spread"
        elif span_years >= 2:
            label = "persistent_local"
        else:
            label = "single_outbreak"

        rows.append(
            {
                "strain": strain,
                "strain_level": g["strain_level"].iloc[0],
                "lineages": _join_unique(g["lineage"]),
                "n_genomes": n_genomes,
                "first_seen": first_seen.date() if not pd.isna(first_seen) else None,
                "last_seen": last_seen.date() if not pd.isna(last_seen) else None,
                "active_years": active_years,
                "n_countries": n_countries,
                "countries": ";".join(countries),
                "n_admin1": n_admin1,
                "origin_country": origin,
                "max_spread_km": max_spread,
                "spread_rate_km_per_year": spread_rate,
                "linked_cases": linked_cases,
                "linked_deaths": linked_deaths,
                "mean_cfr": mean_cfr,
                "behavior_label": label,
            }
        )
    return pd.DataFrame(rows).sort_values("n_genomes", ascending=False)


def build_outbreak_summary(con, species):
    """Summarise historical outbreak totals from yearly-range datasets."""
    query = """
    SELECT
        r.record_date,
        r.country,
        r.admin1,
        r.cases,
        r.deaths,
        r.suspected,
        r.recovered,
        r.indicator_label,
        d.dataset_name AS source_dataset
    FROM epidemiological_records r
    JOIN epidemiological_datasets d ON r.dataset_id = d.dataset_id
    LEFT JOIN geographic_locations gl ON r.location_id = gl.location_id
    WHERE d.species = ? AND d.dataset_type = 'summary'
    """
    df = con.execute(query, [species]).df()
    if df.empty:
        return pd.DataFrame(
            columns=[
                "start_year",
                "country",
                "admin1",
                "cases",
                "deaths",
                "cfr",
                "species_subtype",
                "source_dataset",
            ]
        )
    df = df.copy()
    df["record_date"] = pd.to_datetime(df["record_date"], errors="coerce")
    df = df.dropna(subset=["record_date"])
    df["start_year"] = df["record_date"].dt.year.astype(int)
    df["cfr"] = df["deaths"].astype(float) / df["cases"].astype(float)
    df["cfr"] = df["cfr"].replace([np.inf, -np.inf], np.nan)
    df = df.rename(columns={"indicator_label": "species_subtype"})
    return df[
        [
            "start_year",
            "country",
            "admin1",
            "cases",
            "deaths",
            "cfr",
            "species_subtype",
            "source_dataset",
        ]
    ]


def build_yearly_burden(weekly: pd.DataFrame) -> tuple:
    """Aggregate weekly burden to yearly per country/admin1 and country only."""
    if weekly.empty or "week_start" not in weekly.columns:
        epi_admin = pd.DataFrame(
            columns=["country", "admin1", "year", "epi_cases", "epi_deaths", "epi_trend", "epi_cfr"]
        )
        epi_country = pd.DataFrame(columns=["country", "year", "epi_cases", "epi_deaths", "epi_cfr"])
        return epi_admin, epi_country

    b = weekly[["country", "admin1", "week_start", "new_cases", "new_deaths", "trend_cases"]].copy()
    b["admin1"] = b["admin1"].fillna("National").astype(str)
    b["year"] = pd.to_datetime(b["week_start"], errors="coerce").dt.year
    b = b.dropna(subset=["year"])
    epi_admin = (
        b.groupby(["country", "admin1", "year"])
        .agg(
            epi_cases=("new_cases", "sum"),
            epi_deaths=("new_deaths", "sum"),
            epi_trend=("trend_cases", lambda x: x.mode().iloc[0] if not x.mode().empty else pd.NA),
        )
        .reset_index()
    )
    epi_admin["epi_cfr"] = epi_admin["epi_deaths"] / epi_admin["epi_cases"]
    epi_admin["epi_cfr"] = epi_admin["epi_cfr"].replace([np.inf, -np.inf], np.nan)

    epi_country = (
        b.groupby(["country", "year"])
        .agg(epi_cases=("new_cases", "sum"), epi_deaths=("new_deaths", "sum"))
        .reset_index()
    )
    epi_country["epi_cfr"] = epi_country["epi_deaths"] / epi_country["epi_cases"]
    epi_country["epi_cfr"] = epi_country["epi_cfr"].replace([np.inf, -np.inf], np.nan)
    return epi_admin, epi_country


def _mode_str(series: pd.Series) -> str:
    """Most frequent non-null value in a series ('' if none)."""
    s = series.dropna()
    s = s[~s.isin(["<NA>", "nan", "None", "unassigned", ""])]
    if s.empty:
        return ""
    return str(s.mode().iloc[0])


def _join_unique(series: pd.Series, sep: str = ";") -> str:
    return sep.join(sorted({str(v).strip() for v in series.dropna() if str(v).strip().lower() not in {"", "nan", "<na>", "none"}}))


def build_genome_linkage(
    weekly: pd.DataFrame, con: duckdb.DuckDBPyConnection, species: str
) -> pd.DataFrame:
    """Add genome counts and identifiers to each weekly burden row by month-year and geo."""
    if weekly.empty:
        return weekly

    weekly = weekly.copy()
    weekly["year_month"] = weekly["week_start"].apply(
        lambda d: f"{d.year}-{d.month:02d}" if not pd.isna(d) else None
    )

    samples = con.execute(SAMPLE_QUERY, [species]).df()
    out_cols = [
        "n_genomes",
        "genome_sample_names",
        "genome_clades",
        "genome_lineages",
        "genome_outbreaks",
        "genome_strains",
        "dominant_strain",
        "genome_ppx_accessions",
    ]
    for col in out_cols:
        weekly[col] = 0 if col == "n_genomes" else ""

    if samples.empty:
        weekly = weekly.drop(columns=["year_month"], errors="ignore")
        return weekly

    samples["collection_date"] = pd.to_datetime(samples["collection_date"], errors="coerce")
    samples["country"] = (
        samples["country"]
        .astype(str)
        .replace(["<NA>", "nan", "None"], "Unknown")
        .replace(COUNTRY_NAME_MAP)
        .fillna("Unknown")
    )
    samples["admin1"] = (
        samples["admin1"]
        .astype(str)
        .replace(["<NA>", "nan", "None", ""], "National")
        .fillna("National")
    )
    samples = _fill_country_coords(samples, con)
    samples["sample_name"] = samples["sample_name"].astype(str)
    samples["clade"] = samples["clade"].astype(str).replace(
        ["<NA>", "nan", "None", "unassigned"], np.nan
    )
    samples["lineage"] = samples["lineage"].astype(str).replace(
        ["<NA>", "nan", "None"], np.nan
    )
    samples["outbreak"] = samples["outbreak"].astype(str).replace(
        ["<NA>", "nan", "None", "unassigned"], np.nan
    )
    samples["strain"] = (
        samples["clade"]
        .combine_first(samples["lineage"])
        .combine_first(samples["outbreak"])
    )
    samples["ppx_accession"] = samples["ppx_accession"].astype(str)
    samples["year_month"] = samples.apply(
        lambda r: _to_year_month(r["collection_date"]) if pd.notna(r["collection_date"]) else _to_year_month(r["tip_date"]),
        axis=1,
    )
    samples = samples.dropna(subset=["year_month"])

    if samples.empty:
        weekly = weekly.drop(columns=["year_month"], errors="ignore")
        return weekly

    admin = (
        samples.groupby(["country", "admin1", "year_month"])
        .agg(
            n_genomes=("sample_name", "count"),
            genome_sample_names=("sample_name", _join_unique),
            genome_clades=("clade", _join_unique),
            genome_lineages=("lineage", _join_unique),
            genome_outbreaks=("outbreak", _join_unique),
            genome_strains=("strain", _join_unique),
            dominant_strain=("strain", _mode_str),
            genome_ppx_accessions=("ppx_accession", _join_unique),
        )
        .reset_index()
    )
    country = (
        samples.groupby(["country", "year_month"])
        .agg(
            n_genomes=("sample_name", "count"),
            genome_sample_names=("sample_name", _join_unique),
            genome_clades=("clade", _join_unique),
            genome_lineages=("lineage", _join_unique),
            genome_outbreaks=("outbreak", _join_unique),
            genome_strains=("strain", _join_unique),
            dominant_strain=("strain", _mode_str),
            genome_ppx_accessions=("ppx_accession", _join_unique),
        )
        .reset_index()
    )

    # Merge admin-level first: short column names for admin values
    weekly = weekly.merge(admin, on=["country", "admin1", "year_month"], how="left", suffixes=("", "_admin"))
    # Merge country-level second: country values get the _country suffix
    weekly = weekly.merge(country, on=["country", "year_month"], how="left", suffixes=("", "_country"))

    m_nat = weekly["admin1"] == "National"
    n_admin = pd.to_numeric(weekly["n_genomes_admin"], errors="coerce").fillna(0).astype(int)
    n_country = pd.to_numeric(weekly["n_genomes_country"], errors="coerce").fillna(0).astype(int)
    weekly["n_genomes"] = n_admin
    weekly.loc[m_nat, "n_genomes"] = n_country[m_nat].values
    weekly["n_genomes"] = pd.to_numeric(weekly["n_genomes"], errors="coerce").fillna(0).astype(int)

    for col in [
        "genome_sample_names",
        "genome_clades",
        "genome_lineages",
        "genome_outbreaks",
        "genome_strains",
        "dominant_strain",
        "genome_ppx_accessions",
    ]:
        admin_col = weekly[f"{col}_admin"].astype("object").fillna("").astype(str)
        country_col = weekly[f"{col}_country"].astype("object").fillna("").astype(str)
        weekly[col] = admin_col
        weekly.loc[m_nat, col] = country_col[m_nat].values
        weekly[col] = weekly[col].fillna("").astype(str)

    weekly = weekly.drop(
        columns=[
            c
            for c in weekly.columns
            if c == "year_month" or c.endswith("_admin") or c.endswith("_country")
        ],
        errors="ignore",
    )
    return weekly


def write_tsv(df: pd.DataFrame, path: Path):
    if df.empty:
        # Write empty file with header if columns exist
        df.head(0).to_csv(path, sep="\t", index=False)
    else:
        df.to_csv(path, sep="\t", index=False, float_format="%.6g")


def main():
    parser = argparse.ArgumentParser(description="Pre-compute Transmission & Spread context tables")
    parser.add_argument("--db", default="results/knowledge_warehouse/knowledge_warehouse.duckdb")
    parser.add_argument("--species", default="ebov")
    parser.add_argument("--outdir", default="results/transmission_context")
    args = parser.parse_args()

    outdir = Path(args.outdir) / args.species
    outdir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(args.db, read_only=True)

    # 1. Burden
    epi = con.execute(EPIDEMIOLOGICAL_QUERY, [args.species]).df()
    records = parse_epi_records(epi)
    daily = build_daily_burden(records)
    daily = add_new_counts(daily)
    weekly = weekly_aggregation(daily)
    weekly = add_trend_stats(weekly)
    weekly = add_anomaly_stats(weekly)

    weekly = weekly.rename(
        columns={
            "daily_cases": "weekly_cases_sum",
            "daily_deaths": "weekly_deaths_sum",
            "daily_suspected": "weekly_suspected_sum",
            "daily_recovered": "weekly_recovered_sum",
        }
    )
    burden_cols = [
        "country",
        "admin1",
        "week_start",
        "cases_cum",
        "deaths_cum",
        "new_cases",
        "new_deaths",
        "cfr_cum",
        "cfr_new",
        "cfr_lower",
        "cfr_upper",
        "ma4_new_cases",
        "ma4_new_deaths",
        "growth_rate",
        "mk_p_cases",
        "sen_slope_cases",
        "trend_cases",
        "mk_p_deaths",
        "sen_slope_deaths",
        "mk_p_cfr",
        "sen_slope_cfr",
        "ewma",
        "ewma_ucl",
        "cusum_pos",
        "cusum_neg",
        "alert",
        "n_genomes",
        "genome_sample_names",
        "genome_clades",
        "genome_lineages",
        "genome_outbreaks",
        "genome_strains",
        "dominant_strain",
        "genome_ppx_accessions",
    ]
    # Link each weekly burden row to genomes from the same month-year and location
    weekly = build_genome_linkage(weekly, con, args.species)
    # Stable schema: keep burden_cols order, add any missing as empty
    weekly = weekly.reindex(columns=burden_cols)
    write_tsv(weekly, outdir / "transmission_burden.tsv")

    # 2. Spatial
    spatial = build_spatial(weekly if not weekly.empty else pd.DataFrame(), con)
    write_tsv(spatial, outdir / "transmission_spatial.tsv")

    # 3. Anomaly
    anomaly = build_anomaly(weekly)
    write_tsv(anomaly, outdir / "transmission_anomaly.tsv")

    # 4. Potential
    potential = build_potential(con, args.species, weekly if not weekly.empty else pd.DataFrame())
    write_tsv(potential, outdir / "transmission_potential.tsv")

    # 5. Outbreak summary
    summary = build_outbreak_summary(con, args.species)
    write_tsv(summary, outdir / "transmission_outbreak_summary.tsv")

    # 6. Strain transmission profile (per-clade historical behaviour)
    strain = build_strain_profile(
        con, args.species, weekly if not weekly.empty else pd.DataFrame()
    )
    write_tsv(strain, outdir / "strain_transmission_profile.tsv")

    con.close()
    print(f"Wrote Transmission & Spread context to {outdir}")


if __name__ == "__main__":
    main()
