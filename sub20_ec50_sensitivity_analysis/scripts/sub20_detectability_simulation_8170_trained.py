from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import f


RNG = np.random.default_rng(20260726)

LOWER_BOUNDARY_PEC50_20UM = -np.log10(20e-6)

RAW_8170 = Path("E:/gradthesis/Raw_sq_8170_TMTMosaic.csv")
CURVE_8170 = Path("E:/gradthesis/curvecurator/8170/row/8170_output_curves_row.tsv")

OUT_DIR = Path("C:/Users/owner/Documents/Codex/2026-07-07/ni-ha/outputs/sub20_detectability_simulation_8170_trained")

CURRENT_DOSES = np.array([20, 40, 80, 160, 320, 640], dtype=float)
EXTENDED_DOSES = np.array([1, 2.5, 5, 10, 20, 40, 80, 160, 320, 640], dtype=float)

TRUE_EC50_UM = [1, 5, 10, 15, 20, 40]
TRUE_LOG2FC = [1, 2, 3]
SIGNAL_GROUPS = ["low", "medium", "high"]
N_PER_CELL = 20

TMT_CHANNELS = [
    "default~126_sn_sum",
    "default~127n_sn_sum",
    "default~127c_sn_sum",
    "default~128n_sn_sum",
    "default~128c_sn_sum",
    "default~129n_sn_sum",
    "default~129c_sn_sum",
    "default~130n_sn_sum",
    "default~130c_sn_sum",
    "default~131n_sn_sum",
    "default~131c_sn_sum",
    "default~132n_sn_sum",
    "default~132c_sn_sum",
    "default~133n_sn_sum",
    "default~133c_sn_sum",
    "default~134n_sn_sum",
    "default~134c_sn_sum",
    "default~135_sn_sum",
]


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    valid = ~np.isnan(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    tmp = np.empty(n)
    tmp[order] = adj
    out[valid] = tmp
    return out


def response_model(dose_um, bottom, log2fc, log10_ec50_um, hill):
    # Increasing 4PL-like response. top/bottom ratio is 2^log2fc.
    ec50 = 10 ** log10_ec50_um
    top = bottom * (2 ** log2fc)
    xh = np.power(dose_um, hill)
    eh = np.power(ec50, hill)
    return bottom + (top - bottom) * xh / (eh + xh)


def fit_curve(doses, y_mean):
    y = np.asarray(y_mean, dtype=float)
    if np.any(~np.isfinite(y)) or np.nanmin(y) <= 0:
        return None
    bottom0 = max(np.nanmin(y), 1e-6)
    observed_log2fc0 = np.log2(max(np.nanmax(y) / bottom0, 1.01))
    p0 = [bottom0, min(max(observed_log2fc0, 0.05), 5), np.log10(np.median(doses)), 1.0]
    bounds = ([1e-8, 0.0, np.log10(0.05), 0.2], [np.inf, 6.0, np.log10(2000), 8.0])
    try:
        pars, _ = curve_fit(
            response_model,
            doses,
            y,
            p0=p0,
            bounds=bounds,
            maxfev=3000,
        )
    except Exception:
        return None
    pred = response_model(doses, *pars)
    rss_curve = float(np.sum((y - pred) ** 2))
    rss_null = float(np.sum((y - np.mean(y)) ** 2))
    if rss_null <= 0:
        r2 = np.nan
    else:
        r2 = 1 - rss_curve / rss_null
    n = len(y)
    p_curve = 4
    df1 = 3
    df2 = max(n - p_curve, 1)
    if rss_curve <= 0:
        p_value = 0.0
    else:
        f_value = ((rss_null - rss_curve) / df1) / (rss_curve / df2)
        p_value = float(f.sf(max(f_value, 0), df1, df2))
    bottom, log2fc, log10_ec50_um, hill = pars
    ec50_um = 10 ** log10_ec50_um
    pec50 = -np.log10(ec50_um * 1e-6)
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    observed_range_log2fc = float(np.log2(np.nanmax(y) / np.nanmin(y)))
    return {
        "fit_success": True,
        "fitted_bottom": bottom,
        "fitted_log2FC_total": log2fc,
        "fitted_EC50_uM": ec50_um,
        "fitted_pEC50": pec50,
        "fitted_hill": hill,
        "curve_R2": r2,
        "curve_RMSE": rmse,
        "p_value": p_value,
        "observed_range_log2FC": observed_range_log2fc,
    }


def load_training_distributions():
    raw = pd.read_csv(RAW_8170)
    values = raw[TMT_CHANNELS].apply(pd.to_numeric, errors="coerce")
    raw["overall_signal"] = values.mean(axis=1, skipna=True)

    dose_groups = [TMT_CHANNELS[i * 3 : (i + 1) * 3] for i in range(6)]
    cvs = []
    for cols in dose_groups:
        dose_values = raw[cols].apply(pd.to_numeric, errors="coerce")
        cv = dose_values.std(axis=1, skipna=True) / dose_values.mean(axis=1, skipna=True)
        cvs.append(cv)
    raw["median_cv"] = pd.concat(cvs, axis=1).median(axis=1, skipna=True)
    raw["signal_group"] = pd.qcut(raw["overall_signal"], 3, labels=SIGNAL_GROUPS, duplicates="drop")

    signal_dist = {}
    cv_dist = {}
    for group in SIGNAL_GROUPS:
        sub = raw[raw["signal_group"].astype(str) == group]
        signal_dist[group] = sub["overall_signal"].dropna().clip(lower=1).to_numpy()
        cv_dist[group] = sub["median_cv"].dropna().clip(lower=0.01, upper=1.0).to_numpy()

    curve = pd.read_csv(CURVE_8170, sep="\t")
    slopes = curve.loc[
        (curve["pEC50"].between(-np.log10(640e-6), LOWER_BOUNDARY_PEC50_20UM))
        & (curve["Curve R2"] >= 0.8),
        "Curve Slope",
    ].dropna()
    slopes = slopes.abs().clip(lower=0.3, upper=6).to_numpy()
    if len(slopes) == 0:
        slopes = np.array([1.0])
    return signal_dist, cv_dist, slopes


def simulate_one(design_name, doses, true_ec50, true_log2fc, signal_group, signal_dist, cv_dist, slope_dist, sim_idx):
    bottom = float(RNG.choice(signal_dist[signal_group]))
    cv = float(RNG.choice(cv_dist[signal_group]))
    hill = float(RNG.choice(slope_dist))
    true_mean = response_model(doses, bottom, true_log2fc, np.log10(true_ec50), hill)
    # Log-normal noise with approximate coefficient of variation from 8170.
    sigma = np.sqrt(np.log(cv**2 + 1))
    reps = RNG.lognormal(mean=np.log(true_mean) - 0.5 * sigma**2, sigma=sigma, size=(3, len(doses)))
    y_mean = reps.mean(axis=0)
    fit = fit_curve(doses, y_mean)
    row = {
        "simulation_id": f"{design_name}_EC{true_ec50}_FC{true_log2fc}_{signal_group}_{sim_idx}",
        "design": design_name,
        "true_EC50_uM": true_ec50,
        "true_pEC50": -np.log10(true_ec50 * 1e-6),
        "true_log2FC_total": true_log2fc,
        "true_hill": hill,
        "signal_group": signal_group,
        "baseline_signal": bottom,
        "sampled_cv": cv,
    }
    if fit is None:
        row.update({"fit_success": False})
    else:
        row.update(fit)
    return row


def run_simulation():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    signal_dist, cv_dist, slope_dist = load_training_distributions()

    rows = []
    designs = {
        "current_20_640uM": CURRENT_DOSES,
        "extended_1_640uM": EXTENDED_DOSES,
    }
    for design_name, doses in designs.items():
        for true_ec50 in TRUE_EC50_UM:
            for true_log2fc in TRUE_LOG2FC:
                for signal_group in SIGNAL_GROUPS:
                    for sim_idx in range(N_PER_CELL):
                        rows.append(
                            simulate_one(
                                design_name,
                                doses,
                                true_ec50,
                                true_log2fc,
                                signal_group,
                                signal_dist,
                                cv_dist,
                                slope_dist,
                                sim_idx,
                            )
                        )

    sim = pd.DataFrame(rows)
    sim["p_adj_BH_within_design"] = np.nan
    for design_name, idx in sim.groupby("design").groups.items():
        sim.loc[idx, "p_adj_BH_within_design"] = bh_adjust(sim.loc[idx, "p_value"].to_numpy())

    sim["recovered_sub20_fit"] = sim["fitted_pEC50"] > LOWER_BOUNDARY_PEC50_20UM
    sim["detected_significant"] = (sim["p_adj_BH_within_design"] <= 0.05) & (sim["curve_R2"] >= 0.8)
    sim["strict_detected"] = (
        (sim["p_adj_BH_within_design"] <= 0.05)
        & (sim["curve_R2"] >= 0.99)
        & (sim["observed_range_log2FC"] >= 1)
    )
    sim["sub20_strict_recovered"] = sim["strict_detected"] & sim["recovered_sub20_fit"]

    group_cols = ["design", "true_EC50_uM", "true_log2FC_total", "signal_group"]
    summary = (
        sim.groupby(group_cols)
        .agg(
            n=("simulation_id", "size"),
            fit_success_rate=("fit_success", "mean"),
            sub20_fit_recovery_rate=("recovered_sub20_fit", "mean"),
            significant_detection_rate=("detected_significant", "mean"),
            strict_detection_rate=("strict_detected", "mean"),
            sub20_strict_recovery_rate=("sub20_strict_recovered", "mean"),
            median_fitted_pEC50=("fitted_pEC50", "median"),
            median_observed_range_log2FC=("observed_range_log2FC", "median"),
            median_R2=("curve_R2", "median"),
        )
        .reset_index()
    )

    ec50_summary = (
        sim.groupby(["design", "true_EC50_uM"])
        .agg(
            n=("simulation_id", "size"),
            fit_success_rate=("fit_success", "mean"),
            sub20_fit_recovery_rate=("recovered_sub20_fit", "mean"),
            significant_detection_rate=("detected_significant", "mean"),
            strict_detection_rate=("strict_detected", "mean"),
            sub20_strict_recovery_rate=("sub20_strict_recovered", "mean"),
            median_observed_range_log2FC=("observed_range_log2FC", "median"),
        )
        .reset_index()
    )

    sim.to_csv(OUT_DIR / "simulation_curve_level_results.csv", index=False)
    summary.to_csv(OUT_DIR / "simulation_detection_summary_by_cell.csv", index=False)
    ec50_summary.to_csv(OUT_DIR / "simulation_detection_summary_by_true_EC50.csv", index=False)

    md = [
        "# Sub-20 uM Detectability Simulation",
        "",
        "Training dataset: 8170 raw TMTMosaic data.",
        "",
        f"Current design: {list(CURRENT_DOSES)} uM.",
        f"Extended design: {list(EXTENDED_DOSES)} uM.",
        "",
        f"Simulations per EC50/log2FC/signal/design cell: {N_PER_CELL}.",
        "",
        "Detected sub-20 uM fitted EC50 was defined as fitted pEC50 > 4.699.",
        "",
        "## Summary by true EC50",
        "",
        "```text",
        ec50_summary.to_string(index=False),
        "```",
        "",
        "Interpretation: this simulation estimates detectability under the experimental design. It cannot identify specific real sites with missed sub-20 uM EC50 values.",
    ]
    (OUT_DIR / "simulation_detection_summary.md").write_text("\n".join(md), encoding="utf-8")

    print(ec50_summary.to_string(index=False))
    print(f"\nWrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    run_simulation()
