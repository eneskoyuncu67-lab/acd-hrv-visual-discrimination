# Anticipatory Cardiac Deceleration (ACD) Analysis

MATLAB analysis and figures (including 3D visualizations) of anticipatory
cardiac deceleration, post-stimulus heart-rate response, and heart-rate
variability (HRV), from ECG-derived heartbeat timing around stimulus onset
in a reaction-time task. Tests whether heart rate decelerates before a
stimulus and how that relates to reaction time, accuracy, and experimental
condition, using linear mixed-effects models with a random intercept per
subject.

## How ACD and ARA are indexed

- **ACD (Anticipatory Cardiac Deceleration)**: `HR(stim-3) - HR(stim)`,
  in bpm: the change in heart rate over the 3 heartbeats immediately
  *before* stimulus onset. Positive ACD = heart rate slowing
  (deceleration) as the stimulus approaches. Indexed following:

  > Sroufe, L. A. (1971). Age changes in cardiac deceleration within a
  > fixed foreperiod reaction-time task: An index of attention.
  > *Developmental Psychology*, 5(2), 338-343.
  > https://doi.org/10.1037/h0031418

- **ARA**: `HR(stim+3) - HR(stim+1)`, in bpm: the analogous change in
  heart rate over the 3 heartbeats *after* stimulus onset.

Both measures are computed from R-peak-derived inter-beat intervals (IBI,
converted to bpm via `60000/IBI_ms`), with a 7-criterion exclusion
pipeline applied to both per trial: a physiological IBI range check
(444-1333 ms), an absolute per-subject outlier threshold, and a two-pass
z ± 3 outlier exclusion per beat position and per condition.

## How HRV is computed

Block-level heart-rate variability is computed independently of ACD/ARA,
from the continuous IBI series: time-domain metrics (mean IBI, mean HR,
SDNN, RMSSD, SDSD, pNN50) and frequency-domain metrics (LF, HF, LF/HF
ratio, total power) via Welch power spectral density on a cubic-spline-
interpolated IBI series. A companion 3D visualization plots the
power-spectral-density surface (frequency x time-window x power) per
condition.

## Models

**Predicting ACD** (does pre-stimulus deceleration depend on RT,
accuracy, or condition?)
- `ACD ~ RT + (1|SubjectID)`, fit separately per condition
- `ACD ~ RT + ConditionVal + (1|SubjectID)`, pooled across conditions
- `ACD ~ Acc * Condition + (1|SubjectID)`
- `ACD ~ ConditionVal + (1|SubjectID)`
- `RMSSD ~ ConditionVal + (1|SubjectID)`: HRV as an outcome, does
  condition predict it

**Predicting ARA** (does the post-stimulus response depend on condition,
RT, or accuracy?)
- `ARA ~ ConditionVal + (1|SubjectID)`
- `ARA ~ RT + ConditionVal + (1|SubjectID)`
- `ARA ~ Acc * ConditionVal + (1|SubjectID)`

**Does ACD predict ARA, and is that relationship moderated?**
- `ARA ~ ACD + Acc + ConditionVal + ACD:Acc + Acc:ConditionVal + (1|SubjectID)`
- `ARA ~ ACD * Acc + ConditionVal + (1|SubjectID)`
- `ARA ~ ACD * RT + ConditionVal + (1|SubjectID)` (fit twice, with two
  different visualization styles: a median-split Fast/Slow violin plot,
  and a continuous-RT scatter with a single fitted line)

**Descriptive only** (not a mixed model): Pearson correlation between ACD
and ARA, computed separately for correct vs. incorrect trials within each
condition.

## AI-assisted development

This repository's code was developed with the assistance of **Claude
Code** (Anthropic), used to organize, curate, and modify the original
analysis scripts. 
