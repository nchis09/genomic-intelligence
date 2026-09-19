#!/usr/bin/env python3
"""Populate country-level latitude/longitude coordinates.

Default data source is a CSV of country centroids downloaded from
`https://github.com/gavinr/world-countries-centroids` and stored at
`database/country_centroids.csv`.

If a local shapefile and geopandas are available, use `--shapefile` to compute
polygon centroids from it instead.  If neither source is available, the script
can also derive centroids from existing sub-national `geographic_locations` rows.

Usage:
    python bin/populate_country_coordinates.py \
        --db results/knowledge_warehouse/knowledge_warehouse.duckdb \
        --outdir results/transmission_context
"""

import argparse
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd


CSV_COUNTRY_MAP = {
    "Congo DRC": "Democratic Republic of the Congo",
    "Congo": "Republic of the Congo",
    "United States": "United States of America",
    "United Kingdom": "United Kingdom",
    "England": "United Kingdom",
    "Russian Federation": "Russia",
}

DB_COUNTRY_MAP = {
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


def _centroids_from_shapefile(shapefile: Path) -> pd.DataFrame:
    """Read a world shapefile with geopandas and return country centroids."""
    try:
        import geopandas as gpd
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError(
            "geopandas is not installed. Install it or omit --shapefile "
            "to derive centroids from the DuckDB geographic_locations table."
        ) from exc

    gdf = gpd.read_file(shapefile)
    # Heuristic column names used by Natural Earth / GADM shapefiles
    country_col = None
    for col in ["ADMIN", "SOVEREIGNT", "NAME", "NAME_EN", "COUNTRY"]:
        if col in gdf.columns:
            country_col = col
            break
    if country_col is None:
        raise ValueError(f"Could not find a country name column in {shapefile}; columns: {list(gdf.columns)}")

    gdf["longitude"] = gdf.geometry.centroid.x
    gdf["latitude"] = gdf.geometry.centroid.y
    return gdf[[country_col, "latitude", "longitude"]].rename(columns={country_col: "country"})


def _centroids_from_csv(csv_path: Path) -> pd.DataFrame:
    """Read a country centroids CSV and map names to the warehouse's canonical names."""
    df = pd.read_csv(csv_path)
    df = df.rename(columns={"COUNTRY": "country", "latitude": "latitude", "longitude": "longitude"})
    if "longitude" not in df.columns or "latitude" not in df.columns:
        raise ValueError(f"CSV must contain 'longitude' and 'latitude' columns: {csv_path}")
    df["country"] = df["country"].replace(CSV_COUNTRY_MAP)
    return df[["country", "latitude", "longitude"]]


def _centroids_from_db(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Derive country centroids from existing sub-national geographic_locations."""
    df = con.execute(
        """
        SELECT
            country,
            AVG(latitude) AS latitude,
            AVG(longitude) AS longitude
        FROM geographic_locations
        WHERE latitude IS NOT NULL
          AND longitude IS NOT NULL
        GROUP BY country
        """
    ).df()
    return df


def _country_name_map(con: duckdb.DuckDBPyConnection) -> dict:
    """Return a map from raw country variants to a canonical name if one exists."""
    rows = con.execute(
        "SELECT DISTINCT country FROM geographic_locations WHERE country IS NOT NULL"
    ).df()["country"].tolist()
    return {r: r for r in rows if r}


def main():
    parser = argparse.ArgumentParser(description="Populate country-level coordinates")
    parser.add_argument("--db", required=True, help="Path to the DuckDB knowledge warehouse")
    parser.add_argument("--outdir", required=True, help="Output directory for the TSV")
    parser.add_argument("--shapefile", help="Path to a world shapefile (requires geopandas)")
    parser.add_argument("--centroids-csv", help="Path to a country centroids CSV")
    parser.add_argument("--update-db", action="store_true", help="Update geographic_locations country-level rows")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(args.db)

    if args.shapefile:
        centroids = _centroids_from_shapefile(Path(args.shapefile))
    elif args.centroids_csv:
        centroids = _centroids_from_csv(Path(args.centroids_csv))
    else:
        default_csv = Path(__file__).parent.parent / "database" / "country_centroids.csv"
        if default_csv.exists():
            print(f"Using default centroids CSV: {default_csv}")
            centroids = _centroids_from_csv(default_csv)
        else:
            print("No --shapefile or --centroids-csv given; deriving centroids from geographic_locations sub-national rows.")
            centroids = _centroids_from_db(con)

    if centroids.empty:
        print("No country centroids could be derived.", file=sys.stderr)
        con.close()
        sys.exit(1)

    # Keep only one row per country using the mean of any duplicate entries
    # Normalise to the warehouse's canonical country names so the TSV can be joined
    centroids["country"] = centroids["country"].replace(DB_COUNTRY_MAP)
    centroids = (
        centroids.groupby("country", as_index=False)
        .agg({"latitude": "mean", "longitude": "mean"})
    )
    centroids = centroids.replace([np.inf, -np.inf], np.nan).dropna()

    tsv_path = outdir / "country_coordinates.tsv"
    centroids.to_csv(tsv_path, sep="\t", index=False, float_format="%.6f")
    print(f"Wrote {len(centroids)} country centroids to {tsv_path}")

    if args.update_db:
        con.execute("CREATE TABLE IF NOT EXISTS country_coordinates (country TEXT, latitude DOUBLE, longitude DOUBLE)")
        con.execute("DELETE FROM country_coordinates")
        con.register("centroids_df", centroids)
        con.execute("INSERT INTO country_coordinates SELECT * FROM centroids_df")
        con.unregister("centroids_df")

        # Also backfill country-level rows in geographic_locations where lat/lon are NULL
        existing = con.execute(
            "SELECT DISTINCT country FROM geographic_locations WHERE admin1 IS NULL"
        ).df()["country"].tolist()
        for _, row in centroids.iterrows():
            if row["country"] in existing:
                con.execute(
                    """
                    UPDATE geographic_locations
                    SET latitude = ?, longitude = ?
                    WHERE country = ?
                      AND admin1 IS NULL
                    """,
                    [row["latitude"], row["longitude"], row["country"]],
                )
        print("Updated country_coordinates table and backfilled geographic_locations country-level rows.")

    con.close()


if __name__ == "__main__":
    main()
