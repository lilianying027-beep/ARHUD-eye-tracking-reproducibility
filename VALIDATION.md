# Validation Record

The minimized package was originally tested on 20 August 2026 and re-tested on 25 August 2026 in an isolated working directory without reading from or referencing the raw-data or Tobii-processing directories. The re-test used only the packaged data and code together with the listed software dependencies. This is a released-files reproducibility check, not an external independent validation of the study, its data collection, or its ecological validity.

The five released data files are byte-for-byte unchanged from the package assembled on 20 August 2026. The 25 August update clarifies documentation and adds a provenance record and data checksums; it does not alter the statistical inputs or analytical calculations.

## Headline checks

- Maneuver-level transition results retained after global BH-FDR correction: **14**.
- Participant-mean transition results retained after global BH-FDR correction: **0**.
- Minimum participant-mean transition adjusted p-value: **0.0931243940**.
- Participant-mean AUC scenarios retained after nine-scenario BH-FDR correction: **PTL, PLL, and PLR**.
- Equal-weight mixed-model HUD-minus-RAW contrast: **-6.378231751**, p = **0.126631015**.
- Condition-by-scenario interaction: chi-square(8) = **169.137912653**.

## Comparison with the full internal reproduction workflow

The largest absolute numerical discrepancies between the outputs generated from this minimized package and the outputs from the full internal reproduction workflow were:

| Analysis output | Maximum absolute difference |
|---|---:|
| 81 maneuver-level transition tests | 1.11e-16 |
| Clip-level distance sensitivity | 2.84e-14 |
| Crossed-random-intercept AUC model | 7.82e-14 |
| Participant-mean AUC sensitivity | 0 |
| Participant-mean transition sensitivity | 0 |
| Table 2 retained transitions | 1.11e-16 |
| Table 3 maneuver-level AUC | 2.27e-13 |
| Table 4 frame-summary distance analysis | 5.68e-14 |

These differences are floating-point rounding at machine precision and do not affect any reported value or conclusion.
