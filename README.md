# Autonomic regulation on heart rate and blood pressure

A unified closed-loop cardiovascular model that combines:

1. **ANS-mediated heart-rate modulation** — `Neuro_HR` from `ans_modulation.f90`
2. **ANS-mediated blood-pressure modulation** — `Neuro_BP` from `bp_modulation_reflected.f90`
3. **A simplified vascular-resistance disturbance** — analogous to `Disturbance` in `cvs.f90`
4. **Cardiac inflow with a phenomenological reflected-wave mechanism**
5. **A minimal closed-loop artery–vascular bed–vein circulation**

The integrated code is a **reduced conceptual model**. It is intended for
demonstration of ANS-mediated cardiovascular regulation, for sensitivity
analysis, and for sharing alongside the existing simplified component
models. It does **not** restore the full multi-vessel network or 1D wave
propagation of the original `cvs.f90`.

---

## 1. Files

| File | Role |
|---|---|
| `ans_cvs_integrated.f90` | Integrated Fortran 90/95 model. Independent of `cvs.f90`; the original file is not modified. |
| `plot_ans_cvs.py`        | Python visualisation (`numpy` + `matplotlib`). Reads the three `.dat` outputs. |
| `ANS_CVS_waveform.dat`   | Per-step waveform output (created on each run). |
| `ANS_activity.dat`       | Per-step ANS-activity output (created on each run). |
| `ANS_CVS_summary.dat`    | Per-cardiac-cycle summary output (created on each run). |

---

## 2. Required runtime environment

- A Fortran 90/95-capable compiler (`gfortran` recommended). No external Fortran libraries.
- Python 3 with `numpy` and `matplotlib`. No `pandas` required.
- macOS / Linux: install `gfortran` via the system package manager
  (`brew install gcc`, `apt install gfortran`, etc.).
- Windows: use WSL with the same Linux commands, or install MinGW-w64 and
  use `gfortran` from the MinGW shell.

---

## 3. Compile and run

```bash
gfortran -O2 -Wall ans_cvs_integrated.f90 -o ans_cvs_integrated
./ans_cvs_integrated
python3 plot_ans_cvs.py
```

Total runtime on a modern CPU is a few seconds. The three `.dat` files
and two `.png` files are written to the current working directory.

---

## 4. What is preserved versus what was adapted

The user-facing requirement was: do not change the core ANS regulation
formulas. The following table summarises which symbols are kept exactly
and which differences are only unit / interface adaptations.

| Quantity | Symbolic form preserved? | Notes |
|---|---|---|
| `ns_hr = 1/(1+(hr/miu_hr)^niu_s_hr)` | ✓ verbatim | from `ans_modulation.f90` |
| `np_hr = 1/(1+(hr/miu_hr)^niu_p_hr)` | ✓ verbatim | from `ans_modulation.f90` |
| `EneuHR_hr = AneuHR_hr·ns_hr − BneuHR_hr·np_hr + CneuHR_hr` | ✓ verbatim | from `ans_modulation.f90` |
| HR RK4 update toward `EneuHR_hr` | ✓ verbatim | from `ans_modulation.f90` |
| `ns_bp = 1/(1+(P_aff/miu_p)^niu_s_p)` | ✓ verbatim | from `bp_modulation_reflected.f90` |
| `np_bp = 1/(1+(P_aff/miu_p)^niu_p_p)` | ✓ verbatim | from `bp_modulation_reflected.f90` |
| `R_eq = R0 + AneuR·(ns_p−0.5) − BneuR·(np_p−0.5)` | ✓ verbatim | from `bp_modulation_reflected.f90` |
| First-order Euler relaxation of `R_vb` toward `R_eq` | ✓ verbatim | from `bp_modulation_reflected.f90` |
| Half-sine inflow + delay-buffer reflected copy | ✓ verbatim | from `bp_modulation_reflected.f90` |
| Mass-preserving rescaling `Q_peak ∝ 1/(1+k_reflect)` | ✓ verbatim | from `bp_modulation_reflected.f90` |
| `TneuHR_hr` (HR time constant) | unit clarified | documented as seconds, matching `TneuR` and `tau_aff` units in the BP model; the symbolic update formula is unchanged |
| Afferent BP signal | adapted | uses the low-pass-filtered `P_art` (same as `bp_modulation_reflected.f90`) in place of the three mid-vessel pressures averaged by the original `cvs.f90` |
| `R_vb` recovery dynamics | adapted | uses the first-order Euler form from `bp_modulation_reflected.f90` in place of the anchored exp-decay form `SVB_R_new = EneuR + (2*R0 − EneuR)·exp(−(nstep − 400000)/TneuR)` from `cvs.f90` |
| `Disturbance` | adapted | operates on a single `R_vb` rather than the 11 vessel beds of the original |
| Hard-coded SV perturbation from `bp_modulation_reflected.f90` | removed | SV is fixed; disturbance acts on resistance only |

---

## 5. Output files

All three files are written to the current working directory and start
with a `#`-prefixed header line. Numeric data follow as
space-separated columns. Sampling rate of the per-step files is every
`nprint = 100` integration steps (every 5 ms with the default `dt`).

### 5.1 `ANS_CVS_waveform.dat`

| # | Column | Unit | Meaning |
|---|---|---|---|
| 1 | `time`        | s | simulation time |
| 2 | `HR`          | beats/s | current heart rate, updated each step by `Neuro_HR` |
| 3 | `P_art`       | mmHg | arterial pressure |
| 4 | `P_ven`       | mmHg | venous pressure |
| 5 | `Q_forward`   | mL/s | forward cardiac inflow (half-sine) |
| 6 | `Q_reflected` | mL/s | delayed reflected inflow |
| 7 | `Q_in_total`  | mL/s | total inflow = `Q_forward + Q_reflected` |
| 8 | `Q_out`       | mL/s | artery → vein flow through the vascular bed |
| 9 | `Q_return`    | mL/s | venous return to the reservoir |
| 10 | `R_vb`       | mmHg·s/mL | current vascular-bed resistance |

### 5.2 `ANS_activity.dat`

| # | Column | Unit | Meaning |
|---|---|---|---|
| 1 | `time`        | s | simulation time |
| 2 | `HR_afferent` | beats/s | HR signal seen by `Neuro_HR` (the current `hr`, no filtering — this is exactly what `ans_modulation.f90` uses) |
| 3 | `BP_afferent` | mmHg | low-pass-filtered `P_art` seen by `Neuro_BP` |
| 4 | `ns_hr`       | — | SNS activity for HR (0–1) |
| 5 | `np_hr`       | — | PSNS activity for HR (0–1) |
| 6 | `ns_bp`       | — | SNS activity for BP (0–1) |
| 7 | `np_bp`       | — | PSNS activity for BP (0–1) |
| 8 | `HR_eq`       | beats/s | equilibrium HR computed from `EneuHR_hr` |
| 9 | `R_eq`        | mmHg·s/mL | equilibrium resistance toward which `R_vb` relaxes |
| 10 | `HR`         | beats/s | current HR (same as column 2; printed again for convenience) |
| 11 | `R_vb`       | mmHg·s/mL | current resistance |

### 5.3 `ANS_CVS_summary.dat`

One line per completed cardiac cycle. A cycle ends when `tcr` wraps from
near `tduration` back to near zero.

| # | Column | Unit | Meaning |
|---|---|---|---|
| 1 | `cycle`           | — | integer cycle index |
| 2 | `time`            | s | time at the end of the cycle |
| 3 | `HR`              | beats/s | HR during the cycle |
| 4 | `SBP`             | mmHg | maximum `P_art` within the cycle |
| 5 | `DBP`             | mmHg | minimum `P_art` within the cycle |
| 6 | `MAP`             | mmHg | arithmetic mean of `P_art` within the cycle |
| 7 | `pulse_pressure`  | mmHg | `SBP − DBP` |
| 8 | `R_vb`            | mmHg·s/mL | resistance at the cycle boundary |

---

## 6. Customising parameters

All tunable parameters live in `subroutine parameterset()` (numerical
constants, initial conditions, simulation timing) and inside
`subroutine neuro_hr()` / `subroutine neuro_bp()` (ANS gains, setpoints,
sigmoid slopes — preserved literally from the source files). To make a
change: edit the `.f90` file, recompile, rerun.

### 6.1 Simulation control (`parameterset()`)

| Parameter | Default | Unit | Meaning |
|---|---|---|---|
| `dt` | 5 × 10⁻⁵ | s | time step |
| `nlast` | 1 800 000 | steps | total simulation length (≈90 s) |
| `nprint` | 100 | steps | output every `nprint` steps |
| `n_onset` | 400 000 | step | step at which both ANS modules engage (≈20 s) |
| `n_dist_onset` | 100 000 | step | disturbance window start (≈5 s) |
| `n_dist_end` | 300 000 | step | disturbance window end (≈15 s) |

### 6.2 Heart-rate parameters

| Parameter | Default | Unit | Source / meaning |
|---|---|---|---|
| initial `tduration` | 0.6 | s | initial cardiac cycle (100 bpm) |
| `miu_hr` | 1.25 | beats/s | HR target (75 bpm; from `ans_modulation.f90`) |
| `niu_s_hr` | 7.0 | — | SNS sigmoid slope for HR |
| `niu_p_hr` | −7.0 | — | PSNS sigmoid slope for HR |
| `AneuHR_hr` | 1.0 | beats/s | SNS gain on HR |
| `BneuHR_hr` | 1.0 | beats/s | PSNS gain on HR |
| `CneuHR_hr` | 1.333 | beats/s | baseline offset (≈80 bpm) |
| `TneuHR_hr` | 5.0 | s | HR relaxation time constant |

### 6.3 Cardiac inflow and reflected wave

| Parameter | Default | Unit | Meaning |
|---|---|---|---|
| `t_systole` | 0.30 | s | systolic ejection duration |
| `SV` | 70.0 | mL | stroke volume (FIXED; not modulated) |
| `k_reflect_base` | 0.30 | — | reflection-strength scaling |
| `tau_reflect` | 0.25 | s | round-trip reflection delay |
| `mass_preserve` | `.true.` | logical | if true, rescale `Q_peak` by `1/(1+k_reflect)` so that ∫Q_in dt = SV per cycle |
| `couple_kreflect_R` | `.false.` | logical | if true, scale `k_reflect_eff` by `R_vb / R0_base` |

### 6.4 Vascular parameters

| Parameter | Default | Unit | Meaning |
|---|---|---|---|
| `C_art` | 1.5 | mL/mmHg | arterial compliance |
| `C_ven` | 100.0 | mL/mmHg | venous compliance |
| `R0_base` | 1.0 | mmHg·s/mL | unperturbed baseline resistance |
| `R_ven` | 0.04 | mmHg·s/mL | venous-return resistance |
| `P_RA` | 2.0 | mmHg | reservoir / right-atrial pressure |
| initial `P_art` | 80.0 | mmHg | arterial-pressure initial condition |
| initial `P_ven` | 5.0 | mmHg | venous-pressure initial condition |
| `tau_aff` | 2.0 | s | low-pass time constant for the afferent BP signal |

### 6.5 BP neuromodulation parameters (`neuro_bp()`)

| Parameter | Default | Unit | Meaning |
|---|---|---|---|
| `miu_p` | 95.0 | mmHg | BP setpoint |
| `niu_s_p` | 7.0 | — | SNS sigmoid slope for BP |
| `niu_p_p` | −7.0 | — | PSNS sigmoid slope for BP |
| `AneuR` | 0.6 | mmHg·s/mL | SNS gain on R |
| `BneuR` | 0.0 | mmHg·s/mL | PSNS gain on R |
| `TneuR` | 10.0 | s | R relaxation time constant |

### 6.6 Disturbance parameters

| Parameter | Default | Unit | Meaning |
|---|---|---|---|
| `R_dist_factor` | 1.5 | — | multiplier applied to `R0_base` during the disturbance window. Set to 2.0 to mimic the original `cvs.f90` exactly; set to 1.0 to disable. |

---

## 7. Using `disturbance()` to set different initial BP levels

`disturbance()` operates only on the vascular-bed resistance. Within
the window `n_dist_onset ≤ nstep ≤ n_dist_end` it forces
`R_vb = R_dist_factor · R0_base`. Outside the window, resistance is
governed by the rest of the model (held at `R0_base` during the warm-up
before ANS onset; regulated by `Neuro_BP` afterwards).

Common settings:

- `R_dist_factor = 1.0` — disables the disturbance.
- `R_dist_factor = 1.5` (default) — moderate hypertensive-like
  perturbation; MAP during the window reaches ≈140 mmHg.
- `R_dist_factor = 2.0` — matches the original `cvs.f90` setting; MAP
  during the window can exceed 200 mmHg in this reduced single-bed
  model.
- Shifting `n_dist_onset` / `n_dist_end` controls when and how long the
  perturbation is applied.

Disturbance is intentionally limited to resistance because the user
requested that HR, SV, target HR, target BP, compliance, and venous
pressure remain fixed at this stage. Subsequent dynamic regulation is
left to `Neuro_HR` and `Neuro_BP`.

---

## 8. Turning the reflected wave on or off and tuning it

In `parameterset()`:

- **Disable reflection**: set `k_reflect_base = 0.0d0`. The waveform
  becomes a single primary peak followed by exponential runoff,
  identical to the simplified BP model without reflection.
- **Stronger reflection**: increase `k_reflect_base` toward 0.5. The
  secondary peak grows.
- **Later reflection**: increase `tau_reflect` toward 0.30 s. The
  secondary peak shifts later in the cycle, from late systole into
  early diastole.
- **MAP invariance** under `k_reflect_base` sweeps requires
  `mass_preserve = .true.` (default).
- **Couple reflection to vascular tone**: set
  `couple_kreflect_R = .true.` so that periods of elevated `R_vb` (high
  SNS) increase the apparent reflection strength.

When `k_reflect_base = 0`, every other mechanism in the model continues
to function normally — HR modulation, BP modulation, and the
disturbance all act on the simplified single-peak waveform.

---

## 9. Coupling between HR modulation and BP modulation

HR and BP are coupled implicitly through the cardiovascular system:

- `Neuro_HR` updates `hr` and `tduration` every step from the current HR
  itself.
- `cardiac_inflow()` reads the current `hr` and `tduration` to determine
  the cardiac cycle position `tcr` and the systolic ejection waveform.
  With `mass_preserve = .true.` the cycle-integrated inflow is `SV`
  regardless of HR; changing HR therefore changes cardiac output
  `CO = HR × SV` linearly.
- The arterial-pressure waveform driven by this dynamic inflow then
  shapes `P_afferent`, which `Neuro_BP` reads to set `R_eq` and to relax
  `R_vb`.
- `R_vb` together with `CO` sets the steady-state MAP through
  `MAP ≈ CO · R_vb + P_RA`.

The two ANS modules can therefore work against each other or with each
other. With the default parameters (HR target 75 bpm, BP target 95 mmHg)
the system settles to those exact targets after the warm-up,
disturbance, and ANS onset transients have decayed.

---

## 10. How to interpret the main results

A representative run with the default settings produces the following
sequence (also shown in `fig_ans_cvs_overview.png`):

1. **Warm-up (0–5 s).** HR pinned at 100 bpm, `R_vb` at 1.0. The
   Windkessel relaxes from initial conditions to a quasi-steady
   pulsatile state. SBP/DBP/MAP ≈ 100/79/93 mmHg.
2. **Disturbance window (5–15 s).** `R_vb` is forced to 1.5 by
   `disturbance()`. With CO unchanged (HR still at 100 bpm), MAP rises
   toward ≈170 mmHg.
3. **Disturbance ends, ANS still off (15–20 s).** `R_vb` stays at 1.5
   (neither `Neuro_BP` nor `disturbance()` is acting). BP continues to
   build up toward the new equilibrium for this resistance level.
4. **ANS engages (t = 20 s).** Both `Neuro_HR` and `Neuro_BP` start.
   `Neuro_HR` sees HR (100 bpm) > target (75 bpm) → `ns_hr` falls,
   `np_hr` rises → HR drops toward 75 bpm with time constant
   `TneuHR_hr` ≈ 5 s. `Neuro_BP` sees `P_afferent` ≈ 180 mmHg ≫ 95 mmHg
   → `ns_bp → 0`, `np_bp → 1` → `R_eq → R0_base − 0.5·AneuR ≈ 0.7` →
   `R_vb` relaxes downward with time constant `TneuR` ≈ 10 s.
5. **Recovery (20–40 s).** HR settles at ≈76 bpm. `R_vb` overshoots
   below 1.0 briefly and then comes back. MAP first undershoots ≈85 mmHg
   and then climbs back toward the setpoint.
6. **Steady state (after ≈40 s).** HR ≈ 76 bpm, MAP ≈ 95 mmHg (matches
   setpoint), `R_vb` ≈ 1.0. The closed-loop control hits the target.

Each cardiac cycle still shows the reflected-wave morphology — primary
peak followed by a small secondary peak / shoulder — at all times. See
`fig_ans_cvs_zoom.png` for a single-cycle close-up at steady state.

---

## 11. Simplifications relative to the original `cvs.f90`

Removed or replaced:

- The full 1D vascular tree (`solver_1D`, `junction`, ~50 named vessels).
- The 0D four-chamber heart with elastance functions and the
  Cardiac–1D characteristic-variable coupling.
- The 11 parallel vascular beds → replaced by a single lumped bed.
- The respiratory module `RF` and the intrathoracic-pressure terms.
- The body-weight scaling powers (`RAW**0.45` etc.).
- The anchored exp-decay form of `SVB_R_new` → replaced by a continuous
  first-order ODE.
- The constant-offset form `EneuR = A·ns_p − B·np_p + CneuR` → replaced
  by the deviation form `R_eq = R0 + A·(ns_p − 0.5) − B·(np_p − 0.5)`,
  which guarantees `R_eq = R0` at neutral ANS.

---

## 12. Limitations

- The reflected wave is a phenomenological approximation built from one
  delay buffer; it does **not** represent real wave propagation,
  multiple reflection sites, or the dicrotic notch.
- No pulmonary circulation; no respiratory modulation; no preload
  coupling (SV is prescribed).
- ANS modulation is a phenomenological representation rather than a
  full neurophysiological model.
- Parameter values are tuned for conceptual demonstration and
  sensitivity analysis, not for any specific subject.
- The model must not be used for clinical diagnosis or
  patient-specific hemodynamic prediction.
- During the disturbance window with the default `R_dist_factor = 1.5`,
  MAP transiently reaches ≈170 mmHg; this is a property of the
  single-bed lumped model and is more extreme than the multi-bed
  parallel-resistance equivalent in the original `cvs.f90`. Reduce
  `R_dist_factor` for milder perturbations.

---

## 13. References for the physiological design choices

1. *Westerhof N, Lankhaar J-W, Westerhof BE.* The arterial Windkessel.
   **Med Biol Eng Comput** 47(2):131–141, 2009.
   <https://doi.org/10.1007/s11517-008-0359-2>
2. *Ursino M.* Interaction between carotid baroregulation and the
   pulsating heart: a mathematical model. **Am J Physiol Heart Circ
   Physiol** 275(5):H1733–H1747, 1998.
   <https://doi.org/10.1152/ajpheart.1998.275.5.H1733>
3. *Heldt T, Shim EB, Kamm RD, Mark RG.* Computational modeling of
   cardiovascular response to orthostatic stress. **J Appl Physiol**
   92(3):1239–1254, 2002.
   <https://doi.org/10.1152/japplphysiol.00241.2001>
4. *Guyenet PG.* The sympathetic control of blood pressure.
   **Nat Rev Neurosci** 7(5):335–346, 2006.
   <https://doi.org/10.1038/nrn1902>
5. *Eckberg DL.* Sympathovagal balance: a critical appraisal.
   **Circulation** 96(9):3224–3232, 1997.
   <https://doi.org/10.1161/01.CIR.96.9.3224>
6. *Westerhof N, Sipkema P, van den Bos GC, Elzinga G.* Forward and
   backward waves in the arterial system. **Cardiovasc Res**
   6(6):648–656, 1972.
   <https://doi.org/10.1093/cvr/6.6.648>
7. *Burattini R, Campbell KB.* Modified asymmetric T-tube model to
   infer arterial wave reflection at the aortic root. **IEEE Trans
   Biomed Eng** 36(8):805–814, 1989.
   <https://doi.org/10.1109/10.30806>
