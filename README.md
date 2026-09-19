# Anticipatory Cardiac Deceleration (ACD) Analysis

MATLAB analysis and figures (including 3D visualizations) of anticipatory
cardiac deceleration, post-stimulus heart-rate response, and heart-rate
variability (HRV), from ECG-derived heartbeat timing around stimulus onset
in a reaction-time task. Tests whether heart rate decelerates before a
stimulus and how that relates to reaction time, accuracy, and experimental
condition, using linear mixed-effects models with a random intercept per
subject.

## How ACD is indexed

- **ACD (Anticipatory Cardiac Deceleration)**: `HR(stim-3) - HR(stim)`,
  in bpm: the change in heart rate over the 3 heartbeats immediately
  *before* stimulus onset. Positive ACD = heart rate slowing
  (deceleration) as the stimulus approaches. Indexed following:

  > Sroufe, L. A. (1971). Age changes in cardiac deceleration within a
  > fixed foreperiod reaction-time task: An index of attention.
  > *Developmental Psychology*, 5(2), 338–343.
  > https://doi.org/10.1037/h0031418

- **ARA**: `HR(stim+3) - HR(stim+1)`, in bpm: the analogous change in
  heart rate over the 3 heartbeats *after* stimulus onset.

Both measures are computed from R-peak-derived inter-beat intervals, with
a 7-criterion exclusion pipeline applied to both per trial.

## AI-assisted development

This repository's code was developed with the assistance of **Claude
Code** (Anthropic), used to organize, curate, and modify the original
analysis scripts. All analytical decisions were made and reviewed by the
author.
