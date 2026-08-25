"""Reproduce the manuscript statistics from minimized derived data.

The input files contain de-identified derivatives of archived experimental
observations. This script does not simulate participant observations, generate
replacement data, or modify the released input files.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from patsy import build_design_matrices
from scipy.stats import (
    chi2,
    mannwhitneyu,
    norm,
    t,
    ttest_ind_from_stats,
    ttest_rel,
    wilcoxon,
)
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SCENES = ["TL", "PTL", "LL", "PLL", "S", "PLR", "LR", "PTR", "TR"]
TRANSITIONS = [f"{a}->{b}" for a in "LCR" for b in "LCR"]


def bh_fdr(values: pd.Series | np.ndarray) -> np.ndarray:
    return multipletests(np.asarray(values, dtype=float), method="fdr_bh")[1]


def auc_analyses(out: Path) -> pd.DataFrame:
    data = pd.read_csv(DATA / "auc_trial_anonymized.csv")

    maneuver_rows = []
    for scene in SCENES:
        raw = data.query("scene_code == @scene and condition == 'RAW'")["auc"].to_numpy()
        hud = data.query("scene_code == @scene and condition == 'HUD'")["auc"].to_numpy()
        u, p_value = mannwhitneyu(raw, hud, alternative="two-sided")
        maneuver_rows.append(
            {
                "scenario": scene,
                "n_raw": len(raw),
                "n_hud": len(hud),
                "raw_mean": raw.mean(),
                "hud_mean": hud.mean(),
                "U": u,
                "p": p_value,
                "r_hud_vs_raw": 1 - 2 * u / (len(raw) * len(hud)),
            }
        )
    maneuver = pd.DataFrame(maneuver_rows)
    maneuver["p_fdr_9"] = bh_fdr(maneuver["p"])
    maneuver.to_csv(out / "table3_maneuver_auc.csv", index=False)

    participant_mean = (
        data.groupby(["participant_id", "condition", "scene_code"], as_index=False)["auc"]
        .mean()
    )
    participant_rows = []
    for scene in SCENES:
        paired = (
            participant_mean.query("scene_code == @scene")
            .pivot(index="participant_id", columns="condition", values="auc")
            .dropna()
        )
        w, p_value = wilcoxon(paired["HUD"], paired["RAW"], alternative="two-sided")
        participant_rows.append(
            {"scenario": scene, "n_pairs": len(paired), "W": w, "p": p_value}
        )
    participant = pd.DataFrame(participant_rows)
    participant["p_fdr_9"] = bh_fdr(participant["p"])
    participant.to_csv(out / "participant_mean_auc_sensitivity.csv", index=False)
    return data


def transition_analyses(out: Path) -> None:
    maneuver_data = pd.read_csv(DATA / "transition_maneuver_anonymized.csv")
    rows = []
    for scene in SCENES:
        for transition in TRANSITIONS:
            raw = maneuver_data.query(
                "scene_code == @scene and condition == 'RAW' and transition == @transition"
            )["probability"].to_numpy()
            hud = maneuver_data.query(
                "scene_code == @scene and condition == 'HUD' and transition == @transition"
            )["probability"].to_numpy()
            if np.allclose(raw, raw[0]) and np.allclose(hud, hud[0]) and np.isclose(raw[0], hud[0]):
                u, p_value = len(raw) * len(hud) / 2, 1.0
            else:
                u, p_value = mannwhitneyu(raw, hud, alternative="two-sided")
            rows.append(
                {
                    "scenario": scene,
                    "transition": transition,
                    "raw_mean": raw.mean(),
                    "hud_mean": hud.mean(),
                    "U": u,
                    "p": p_value,
                    "r_hud_vs_raw": 1 - 2 * u / (len(raw) * len(hud)),
                }
            )
    maneuver = pd.DataFrame(rows)
    maneuver["p_fdr_81"] = bh_fdr(maneuver["p"])
    maneuver.to_csv(out / "all_81_maneuver_transition_tests.csv", index=False)
    maneuver.query("p_fdr_81 < 0.05").to_csv(
        out / "table2_significant_transitions.csv", index=False
    )

    # Round-trip parsing preserves the exact participant-mean values written by
    # the minimization step. This matters for Wilcoxon tie handling.
    participant_data = pd.read_csv(
        DATA / "transition_participant_aggregate.csv", float_precision="round_trip"
    )
    rows = []
    for scene in SCENES:
        for transition in TRANSITIONS:
            paired = (
                participant_data.query(
                    "scene_code == @scene and transition == @transition"
                )
                .pivot(index="participant_id", columns="condition", values="probability")
                .dropna()
            )
            if np.allclose(paired["HUD"], paired["RAW"]):
                w, p_value = 0.0, 1.0
            else:
                w, p_value = wilcoxon(
                    paired["HUD"], paired["RAW"], alternative="two-sided"
                )
            rows.append(
                {
                    "scenario": scene,
                    "transition": transition,
                    "n_pairs": len(paired),
                    "W": w,
                    "p": p_value,
                }
            )
    participant = pd.DataFrame(rows)
    participant["p_fdr_81"] = bh_fdr(participant["p"])
    participant.to_csv(out / "participant_mean_transition_sensitivity.csv", index=False)


def distance_analyses(out: Path) -> None:
    frame = pd.read_csv(DATA / "distance_frame_summary.csv")
    frame_rows = []
    for scene in SCENES:
        cell = frame.query("scenario_code == @scene").set_index("condition")
        hud, raw = cell.loc["HUD"], cell.loc["RAW"]
        equal_var = bool(hud["equal_variance_used"])
        test = ttest_ind_from_stats(
            mean1=hud["mean_distance_px"],
            std1=hud["sd_distance_px"],
            nobs1=int(hud["n_frames"]),
            mean2=raw["mean_distance_px"],
            std2=raw["sd_distance_px"],
            nobs2=int(raw["n_frames"]),
            equal_var=equal_var,
        )
        pooled = np.sqrt(
            (hud["sd_distance_px"] ** 2 + raw["sd_distance_px"] ** 2) / 2
        )
        frame_rows.append(
            {
                "scenario": scene,
                "hud_mean": hud["mean_distance_px"],
                "raw_mean": raw["mean_distance_px"],
                "cohens_d": (hud["mean_distance_px"] - raw["mean_distance_px"]) / pooled,
                "equal_variance_used": equal_var,
                "p": test.pvalue,
            }
        )
    pd.DataFrame(frame_rows).to_csv(out / "table4_frame_distance.csv", index=False)

    clip = pd.read_csv(DATA / "distance_clip_summary.csv")
    rows = []
    for scene in SCENES:
        paired = (
            clip.query("scenario_code == @scene")
            .pivot(index="source_clip_id", columns="condition", values="mean_distance_px")
            .dropna()
        )
        differences = paired["HUD"] - paired["RAW"]
        test = ttest_rel(paired["HUD"], paired["RAW"])
        standard_error = differences.std(ddof=1) / np.sqrt(len(differences))
        critical = t.ppf(0.975, len(differences) - 1)
        rows.append(
            {
                "scenario": scene,
                "n_clips": len(differences),
                "hud_minus_raw": differences.mean(),
                "ci_low": differences.mean() - critical * standard_error,
                "ci_high": differences.mean() + critical * standard_error,
                "p": test.pvalue,
            }
        )
    clip_result = pd.DataFrame(rows)
    clip_result["p_fdr_9"] = bh_fdr(clip_result["p"])
    clip_result.to_csv(out / "clip_paired_distance_sensitivity.csv", index=False)


def crossed_random_intercept_model(data: pd.DataFrame, out: Path) -> None:
    data = data.copy()
    data["condition"] = pd.Categorical(data["condition"], categories=["RAW", "HUD"])
    data["scene_code"] = pd.Categorical(data["scene_code"], categories=SCENES)
    data["one_group"] = "all"
    model = smf.mixedlm(
        "auc ~ C(condition) * C(scene_code)",
        data,
        groups=data["one_group"],
        re_formula="0",
        vc_formula={
            "participant": "0+C(participant_id)",
            "source_clip": "0+C(source_clip_id)",
        },
    )
    result = model.fit(reml=True, method="lbfgs", maxiter=2000, disp=False)
    fixed_names = list(result.fe_params.index)
    beta = result.fe_params.to_numpy()
    covariance = result.cov_params().loc[fixed_names, fixed_names].to_numpy()

    cells = pd.DataFrame(
        [(condition, scene) for scene in SCENES for condition in ("RAW", "HUD")],
        columns=["condition", "scene_code"],
    )
    cells["condition"] = pd.Categorical(cells["condition"], categories=["RAW", "HUD"])
    cells["scene_code"] = pd.Categorical(cells["scene_code"], categories=SCENES)
    design = np.asarray(build_design_matrices([model.data.design_info], cells)[0])
    contrast = design[cells["condition"] == "HUD"].mean(0) - design[
        cells["condition"] == "RAW"
    ].mean(0)
    estimate = float(contrast @ beta)
    standard_error = float(np.sqrt(contrast @ covariance @ contrast))

    interaction_indices = [
        index for index, name in enumerate(fixed_names) if ":" in name
    ]
    interaction_beta = beta[interaction_indices]
    interaction_covariance = covariance[np.ix_(interaction_indices, interaction_indices)]
    wald = float(
        interaction_beta
        @ np.linalg.pinv(interaction_covariance)
        @ interaction_beta
    )
    variance = dict(zip(model.exog_vc.names, result.vcomp))
    total_variance = result.scale + sum(result.vcomp)
    pd.DataFrame(
        [
            {
                "n": result.nobs,
                "converged": result.converged,
                "participant_icc": variance["participant"] / total_variance,
                "source_clip_icc": variance["source_clip"] / total_variance,
                "equal_weight_hud_minus_raw": estimate,
                "ci_low": estimate - 1.96 * standard_error,
                "ci_high": estimate + 1.96 * standard_error,
                "p": 2 * norm.sf(abs(estimate / standard_error)),
                "interaction_chi2": wald,
                "interaction_df": len(interaction_indices),
                "interaction_p": chi2.sf(wald, len(interaction_indices)),
            }
        ]
    ).to_csv(out / "crossed_random_intercept_auc_model.csv", index=False)


def write_summary(out: Path) -> None:
    auc = pd.read_csv(out / "participant_mean_auc_sensitivity.csv")
    transitions = pd.read_csv(out / "all_81_maneuver_transition_tests.csv")
    participant_transitions = pd.read_csv(
        out / "participant_mean_transition_sensitivity.csv"
    )
    model = pd.read_csv(out / "crossed_random_intercept_auc_model.csv").iloc[0]
    summary = pd.DataFrame(
        [
            ["maneuver_transition_FDR_significant_count", int((transitions.p_fdr_81 < 0.05).sum())],
            ["participant_transition_FDR_significant_count", int((participant_transitions.p_fdr_81 < 0.05).sum())],
            ["participant_transition_min_adjusted_p", participant_transitions.p_fdr_81.min()],
            ["participant_mean_AUC_FDR_significant_scenarios", ",".join(auc.loc[auc.p_fdr_9 < 0.05, "scenario"])],
            ["mixed_model_equal_weight_HUD_minus_RAW", model.equal_weight_hud_minus_raw],
            ["mixed_model_equal_weight_p", model.p],
            ["mixed_model_interaction_chi2", model.interaction_chi2],
        ],
        columns=["check", "value"],
    )
    summary.to_csv(out / "reproduction_summary.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Reproduce the reported statistical analyses from minimal derived data."
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "outputs",
        help="Directory for generated result tables (default: package/outputs).",
    )
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    auc = auc_analyses(out)
    transition_analyses(out)
    distance_analyses(out)
    crossed_random_intercept_model(auc, out)
    write_summary(out)
    print(f"Reproduction completed successfully: {out}")


if __name__ == "__main__":
    main()
