# Minimal Statistical Reproducibility Package

This package contains the minimum de-identified, derived data needed to reproduce the statistical analyses reported in the AR-HUD eye-tracking manuscript. It is not a release of the original eye-tracking recordings.

## Disclosure scope

The package excludes:

- participant names and original participant numbers;
- recording dates and session timestamps;
- original 120-Hz gaze streams and fixation coordinates;
- Tobii Pro Lab project files and databases;
- road-video and experimental-stimulus files;
- calibration records, questionnaires, and demographic source forms;
- raw-data extraction, participant-mapping, and Tobii-processing code.

Participant and source-clip identifiers in the released files are newly generated pseudonyms. The participant labels do not preserve gaps in the original recording names.

## Files

### `data/auc_trial_anonymized.csv`

Minimum input for the maneuver-level AUC comparison, participant-mean AUC sensitivity analysis, and crossed-random-intercept AUC model.

| Field | Meaning |
|---|---|
| `participant_id` | Newly assigned participant pseudonym |
| `source_clip_id` | Source-clip pseudonym |
| `condition` | `HUD` or `RAW` |
| `scene_code` | Maneuver scenario code |
| `auc` | Simpson-integrated time-normalized horizontal-fixation AUC |

### `data/transition_maneuver_anonymized.csv`

Minimum input for the 81 exploratory maneuver-level Mann-Whitney comparisons. Participant, clip, trial, coordinate, and timing identifiers are deliberately omitted.

### `data/transition_participant_aggregate.csv`

Participant-mean transition probabilities required for the paired Wilcoxon sensitivity analysis. These are aggregated values, not gaze events.

### `data/distance_clip_summary.csv`

Mean fixation-to-arrow-reference distance for each pseudonymous source clip, condition, and scenario. It reproduces the clip-level paired analysis without releasing frame coordinates.

### `data/distance_frame_summary.csv`

Only the sample size, mean, standard deviation, and variance assumption used for each scenario-condition cell. These summaries reproduce Table 4's exploratory t-tests and Cohen's d values without releasing 60,710 frame observations.

### `manifest.csv`

Row counts and disclosure level of each released data file.

## Reproduction

Install Python 3.10 or later and the listed dependencies:

```bash
python -m pip install -r requirements.txt
python reproduce_statistics.py
```

To place results outside the package directory:

```bash
python reproduce_statistics.py --output-dir reproduced_outputs
```

The script produces the manuscript statistical tables and sensitivity results, including:

- the 81 transition tests and the 14 maneuver-level results retained after global BH-FDR correction;
- maneuver-level and participant-mean AUC analyses;
- the participant-mean transition sensitivity analysis;
- frame-summary and clip-level distance analyses; and
- the crossed participant/source-clip random-intercept AUC model.

## Important analytical definition

The transition-region definition used for the released results is:

- Left: `x < 700` pixels;
- Center: `700 <= x < 1220` pixels;
- Right: `x >= 1220` pixels.

The released transition probabilities were calculated before data minimization. Raw fixation coordinates are not included in this package.

## Data governance

The original 120-Hz eye-tracking streams, timestamped recording metadata, Tobii project files, video materials, and participant source documents remain under controlled access. Any access statement used in the manuscript should be consistent with the participant-consent documents and the responsible institution's ethics and data-governance requirements.

## Software environment used for validation

The package was independently validated with Python, pandas, NumPy, SciPy, statsmodels, and patsy. Small last-digit differences may occur across library versions, particularly during mixed-model optimization.
