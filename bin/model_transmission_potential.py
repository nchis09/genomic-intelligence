#!/usr/bin/env python3
"""Model transmission potential from the raw query-to-background linkage table.

Algorithm: spatiotemporal kernel (gravity-style) model.

    potential_score = sum_i [ w_i * background_annual_cases_i ]
    w_i = exp(-alpha * div_diff_i)          # genetic proximity
        * exp(-beta  * geo_km_i)            # geographic proximity
        * exp(-gamma * |days_i|)            # temporal proximity

Weights are normalised per query.  The weighted-case score is min-max
normalised across queries to [0, 1] and labelled low/medium/high by
tertiles.  Queries whose neighbours have no linked epi burden are
labelled ``no_epi_link``.

Usage:
    python bin/model_transmission_potential.py \
        --input for_validation/transmission_context/ebov/transmission_potential.tsv \
        --outdir for_validation/transmission_context/ebov
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

SCORE_COLUMNS = [
    "query_sample",
    "query_collection_date",
    "query_country",
    "query_admin1",
    "n_neighbors",
    "n_neighbors_with_epi",
    "mean_div_diff",
    "min_geo_km",
    "potential_cases_weighted",
    "potential_deaths_weighted",
    "potential_cfr_weighted",
    "transmission_potential_score",
    "severity_score",
    "likely_source_countries",
    "predicted_behavior",
    "risk_label",
]


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


def _risk_label(score: pd.Series) -> pd.Series:
    """Tertile-based risk labels for a score column."""
    labels = pd.Series("unknown", index=score.index, dtype=object)
    valid = score.dropna()
    if len(valid) < 3:
        return labels
    q1, q2 = valid.quantile([1 / 3, 2 / 3])
    # Tied scores (e.g. all-zero when no epi links) collapse the tertile edges;
    # dedupe bins and shrink the label set to match so pd.cut doesn't fail.
    edges = sorted({-np.inf, q1, q2, np.inf})
    n_bins = len(edges) - 1
    use_labels = {3: ["low", "medium", "high"], 2: ["low", "high"], 1: ["medium"]}[n_bins]
    labels.loc[score.notna()] = pd.cut(
        score[score.notna()],
        bins=edges,
        labels=use_labels,
    ).astype(str)
    return labels


def model_potential(
    df: pd.DataFrame,
    alpha: float,
    beta: float,
    gamma: float,
    strain_profile: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Compute per-query transmission potential scores from the raw link table."""
    if df.empty:
        return pd.DataFrame(columns=SCORE_COLUMNS)

    df = df.copy()

    # --- distances ---------------------------------------------------------
    if "geo_distance_km" not in df.columns:
        df["geo_distance_km"] = _haversine_km(
            df["query_latitude"],
            df["query_longitude"],
            df["background_latitude"],
            df["background_longitude"],
        )
    if "temporal_distance_days" not in df.columns:
        qd = pd.to_datetime(df["query_collection_date"], errors="coerce")
        bd = pd.to_datetime(df["background_collection_date"], errors="coerce")
        df["temporal_distance_days"] = (qd - bd).abs().dt.days

    # --- kernel weights ------------------------------------------------------
    # Missing distances are treated as "far" (max observed) so weight -> ~0.
    div = df["background_div_diff"].fillna(df["background_div_diff"].max() or 1.0)
    geo = df["geo_distance_km"].fillna(df["geo_distance_km"].max() or 20000.0)
    tim = df["temporal_distance_days"].fillna(df["temporal_distance_days"].max() or 3650.0)

    df["w_genetic"] = np.exp(-alpha * div)
    df["w_geo"] = np.exp(-beta * geo)
    df["w_time"] = np.exp(-gamma * tim)
    df["weight"] = df["w_genetic"] * df["w_geo"] * df["w_time"]

    wsum = df.groupby("query_sample")["weight"].transform("sum")
    df["weight_norm"] = np.where(wsum > 0, df["weight"] / wsum, 0.0)

    # --- weighted burden -----------------------------------------------------
    cases = df["background_annual_cases"].fillna(0)
    deaths = df["background_annual_deaths"].fillna(0)
    cfr = df["background_annual_cfr"].fillna(0)

    df["weighted_cases"] = df["weight_norm"] * cases
    df["weighted_deaths"] = df["weight_norm"] * deaths
    df["weighted_cfr"] = df["weight_norm"] * cfr

    if "background_epi_available" in df.columns:
        epi_flag = df["background_epi_available"].astype(str).str.lower().isin(["true", "1"])
    else:
        epi_flag = df["background_annual_cases"].notna()

    epi_df = df[epi_flag].copy()

    agg = (
        df.assign(has_epi=epi_flag)
        .groupby("query_sample")
        .agg(
            query_collection_date=("query_collection_date", "first"),
            query_country=("query_country", "first"),
            query_admin1=("query_admin1", "first"),
            n_neighbors=("background_sample", "count"),
            n_neighbors_with_epi=("has_epi", "sum"),
            potential_cases_weighted=("weighted_cases", "sum"),
            potential_deaths_weighted=("weighted_deaths", "sum"),
            potential_cfr_weighted=("weighted_cfr", "sum"),
            mean_div_diff=("background_div_diff", "mean"),
            min_geo_km=("geo_distance_km", "min"),
        )
        .reset_index()
    )

    # --- likely source countries (top weighted epi-linked origins) ------------
    if not epi_df.empty:
        src = (
            epi_df.groupby(["query_sample", "background_country"])["weight_norm"]
            .sum()
            .reset_index()
            .sort_values(["query_sample", "weight_norm"], ascending=[True, False])
            .groupby("query_sample")["background_country"]
            .apply(lambda s: ";".join(s.head(3).astype(str)))
        )
        agg = agg.merge(
            src.rename("likely_source_countries"), on="query_sample", how="left"
        )
    else:
        agg["likely_source_countries"] = np.nan

    # --- predicted behaviour from closest epi-linked strain -------------------
    strain_col = (
        "background_strain" if "background_strain" in df.columns
        else "background_clade" if "background_clade" in df.columns
        else None
    )
    if (
        strain_profile is not None
        and not strain_profile.empty
        and strain_col is not None
        and not epi_df.empty
    ):
        label_map = dict(
            zip(strain_profile["strain"], strain_profile["behavior_label"])
        )
        epi_df["bg_behavior"] = epi_df[strain_col].map(label_map)
        beh = (
            epi_df.dropna(subset=["bg_behavior"])
            .groupby(["query_sample", "bg_behavior"])["weight_norm"]
            .sum()
            .reset_index()
            .sort_values(["query_sample", "weight_norm"], ascending=[True, False])
            .groupby("query_sample")["bg_behavior"]
            .first()
        )
        agg = agg.merge(
            beh.rename("predicted_behavior"), on="query_sample", how="left"
        )
    else:
        agg["predicted_behavior"] = np.nan

    # --- scores ---------------------------------------------------------------
    pc = agg["potential_cases_weighted"]
    if pc.max() > pc.min():
        agg["transmission_potential_score"] = (pc - pc.min()) / (pc.max() - pc.min())
    else:
        agg["transmission_potential_score"] = 0.0

    agg["severity_score"] = agg["potential_cfr_weighted"].clip(0, 1)

    # Queries with no epi-linked neighbours cannot be scored
    no_epi = agg["n_neighbors_with_epi"] == 0
    agg.loc[no_epi, "transmission_potential_score"] = np.nan
    agg["risk_label"] = _risk_label(agg["transmission_potential_score"])
    agg.loc[no_epi, "risk_label"] = "no_epi_link"

    return agg[SCORE_COLUMNS]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Model transmission potential (spatiotemporal kernel)."
    )
    ap.add_argument("--input", required=True, help="Path to transmission_potential.tsv")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--alpha", type=float, default=50.0, help="Genetic decay (per unit div)")
    ap.add_argument("--beta", type=float, default=0.001, help="Geographic decay (per km)")
    ap.add_argument("--gamma", type=float, default=0.01, help="Temporal decay (per day)")
    ap.add_argument(
        "--strain-profile",
        default=None,
        help="Path to strain_transmission_profile.tsv (default: same dir as --input)",
    )
    args = ap.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Optional strain profile for predicted behaviour
    sp_path = (
        Path(args.strain_profile)
        if args.strain_profile
        else Path(args.input).parent / "strain_transmission_profile.tsv"
    )
    strain_profile = pd.read_csv(sp_path, sep="\t") if sp_path.exists() else None

    out = model_potential(df, args.alpha, args.beta, args.gamma, strain_profile)
    out_path = outdir / "transmission_potential_score.tsv"
    out.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(out)} query scores to {out_path}")


if __name__ == "__main__":
    main()
