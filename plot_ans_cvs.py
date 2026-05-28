"""
plot_ans_cvs.py

Visualisation script for the output of ans_cvs_integrated.f90.

Reads:
    ANS_CVS_waveform.dat  time HR P_art P_ven Q_forward Q_reflected
                          Q_in_total Q_out Q_return R_vb
    ANS_activity.dat      time HR_aff BP_aff ns_hr np_hr ns_bp np_bp
                          HR_eq R_eq HR R_vb
    ANS_CVS_summary.dat   cycle time HR SBP DBP MAP PP R_vb

Produces:
    fig_ans_cvs_overview.png         (6-panel full timeline)
    fig_ans_cvs_zoom.png             (zoom on a single cycle, reflected wave)
"""

import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------
wave = np.loadtxt("ANS_CVS_waveform.dat", comments="#")
act  = np.loadtxt("ANS_activity.dat",     comments="#")
summ = np.loadtxt("ANS_CVS_summary.dat",  comments="#")

t   = wave[:, 0];  HR_w   = wave[:, 1];  P_art = wave[:, 2]
P_ven = wave[:, 3]; Q_fwd = wave[:, 4];  Q_ref = wave[:, 5]
Q_in  = wave[:, 6]; Q_out = wave[:, 7];  Q_ret = wave[:, 8]
R_vb_w = wave[:, 9]

ta = act[:, 0]
ns_hr = act[:, 3]; np_hr = act[:, 4]
ns_bp = act[:, 5]; np_bp = act[:, 6]
HR_eq = act[:, 7]; R_eq  = act[:, 8]

cyc_t = summ[:, 1]; HR_c = summ[:, 2]
SBP   = summ[:, 3]; DBP  = summ[:, 4]
MAP   = summ[:, 5]; PP   = summ[:, 6]; R_vb_c = summ[:, 7]

# Vertical-line landmarks (match parameterset() defaults)
t_dist_on  = 5.0     # nstep 100000 * dt(5e-5) = 5 s
t_dist_off = 15.0    # nstep 300000
t_onset    = 20.0    # nstep 400000


# ----------------------------------------------------------------------
# Figure 1 — overview: 6 panels
# ----------------------------------------------------------------------
fig, ax = plt.subplots(6, 1, figsize=(12, 16), sharex=True)

# (1) HR
ax[0].plot(t, HR_w, lw=0.5, color="steelblue", alpha=0.5,
           label="HR (per step)")
ax[0].plot(cyc_t, HR_c, "o-", ms=3, color="navy",
           label="HR (per cycle)")
ax[0].set_ylabel("HR (beats/s)")
ax[0].set_title("Integrated ANS–CVS model: HR is dynamically updated by Neuro_HR")
ax[0].legend(loc="lower left")

# (2) Arterial pressure waveform + cycle SBP/DBP/MAP
ax[1].plot(t, P_art, lw=0.4, color="steelblue", alpha=0.7,
           label="P_art (waveform)")
ax[1].plot(cyc_t, SBP, "o-", ms=2, color="firebrick",  label="SBP")
ax[1].plot(cyc_t, DBP, "o-", ms=2, color="darkorange", label="DBP")
ax[1].plot(cyc_t, MAP, "o-", ms=2, color="black",     label="MAP")
ax[1].set_ylabel("Pressure (mmHg)")
ax[1].legend(loc="upper right", fontsize=8, ncol=2)

# (3) Cardiac flow decomposition
ax[2].plot(t, Q_fwd, lw=0.4, color="seagreen", alpha=0.6,
           label="Q_forward")
ax[2].plot(t, Q_ref, lw=0.4, color="purple",   alpha=0.8,
           label="Q_reflected")
ax[2].set_ylabel("Flow (mL/s)")
ax[2].legend(loc="upper right", fontsize=8)

# (4) ANS HR activities
ax[3].plot(ta, ns_hr, color="firebrick", label="ns_hr (SNS, HR)")
ax[3].plot(ta, np_hr, color="royalblue", label="np_hr (PSNS, HR)")
ax[3].axhline(0.5, color="k", lw=0.5, alpha=0.5)
ax[3].set_ylabel("Activity")
ax[3].legend(loc="upper right", fontsize=8)

# (5) ANS BP activities
ax[4].plot(ta, ns_bp, color="firebrick", label="ns_bp (SNS, BP)")
ax[4].plot(ta, np_bp, color="royalblue", label="np_bp (PSNS, BP)")
ax[4].axhline(0.5, color="k", lw=0.5, alpha=0.5)
ax[4].set_ylabel("Activity")
ax[4].legend(loc="upper right", fontsize=8)

# (6) Vascular bed resistance
ax[5].plot(t, R_vb_w, color="black", lw=1.2, label="R_vb (current)")
ax[5].plot(ta, R_eq, color="orange", ls="--", lw=1.0,
           label="R_eq (target)")
ax[5].set_ylabel("R_vb (mmHg·s/mL)")
ax[5].set_xlabel("Time (s)")
ax[5].legend(loc="upper right", fontsize=8)

# Annotate landmarks on every panel
for a in ax:
    a.axvline(t_dist_on,  ls=":",  color="gray", alpha=0.7)
    a.axvline(t_dist_off, ls=":",  color="gray", alpha=0.7)
    a.axvline(t_onset,    ls="--", color="gray", alpha=0.7)
    a.grid(alpha=0.2)
ax[0].text(t_dist_on,  ax[0].get_ylim()[1]*0.97,
           " disturbance on",   fontsize=8, va="top", color="gray")
ax[0].text(t_dist_off, ax[0].get_ylim()[1]*0.97,
           " disturbance off",  fontsize=8, va="top", color="gray")
ax[0].text(t_onset,    ax[0].get_ylim()[1]*0.97,
           " ANS onset",        fontsize=8, va="top", color="gray")

plt.tight_layout()
plt.savefig("fig_ans_cvs_overview.png", dpi=130)
print("Saved fig_ans_cvs_overview.png")


# ----------------------------------------------------------------------
# Figure 2 — single-cycle zoom showing the reflected wave
# Pick a late-time window when HR has settled and disturbance has cleared.
# ----------------------------------------------------------------------
m_late = (t > 80.0) & (t < 81.6)
fig2, ax2 = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
ax2[0].plot(t[m_late], P_art[m_late], color="steelblue", lw=1.4,
            label="P_art")
ax2[0].set_ylabel("P_art (mmHg)")
ax2[0].set_title("Single cycles at steady state: primary peak + secondary peak")
ax2[0].legend(loc="upper right"); ax2[0].grid(alpha=0.3)

ax2[1].plot(t[m_late], Q_fwd[m_late], color="seagreen", lw=1.2,
            label="Q_forward")
ax2[1].plot(t[m_late], Q_ref[m_late], color="purple",   lw=1.2,
            ls="--", label="Q_reflected")
ax2[1].plot(t[m_late], Q_in[m_late],  color="black",    lw=0.7,
            label="Q_in_total")
ax2[1].set_ylabel("Flow (mL/s)"); ax2[1].set_xlabel("Time (s)")
ax2[1].legend(loc="upper right", ncol=3); ax2[1].grid(alpha=0.3)

plt.tight_layout()
plt.savefig("fig_ans_cvs_zoom.png", dpi=130)
print("Saved fig_ans_cvs_zoom.png")
