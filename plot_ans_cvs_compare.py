"""
plot_ans_cvs_compare.py

Before/after comparison for the ANS-onset BP-spike fix in
ans_cvs_integrated.f90.

Expects four data files produced by running the OLD (pre-fix) and NEW
(post-fix) binaries and renaming their outputs:

    ANS_CVS_waveform_before.dat   ANS_activity_before.dat
    ANS_CVS_waveform_after.dat    ANS_activity_after.dat

Waveform columns:
    time HR P_art P_ven Q_forward Q_reflected Q_in_total Q_out Q_return R_vb
Activity columns (after-fix has an extra ans_gain column at the end):
    time HR_aff BP_aff ns_hr np_hr ns_bp np_bp HR_eq R_eq HR R_vb [ans_gain]

Produces:
    fig_spike_fix_compare.png   (before vs after, focused on onset)
    fig_spike_fix_detail.png    (after-fix mechanism detail incl. ans_gain)
"""

import numpy as np
import matplotlib.pyplot as plt

WCOLS = ["time", "HR", "P_art", "P_ven", "Q_fwd", "Q_ref",
         "Q_in", "Q_out", "Q_ret", "R_vb"]


def load_wave(path):
    d = np.loadtxt(path, comments="#")
    return {c: d[:, i] for i, c in enumerate(WCOLS)}


def load_act(path):
    d = np.loadtxt(path, comments="#")
    out = dict(time=d[:, 0], HR=d[:, 9], R_eq=d[:, 8], R_vb=d[:, 10],
               ns_bp=d[:, 5], np_bp=d[:, 6])
    out["ans_gain"] = d[:, 11] if d.shape[1] > 11 else np.zeros(len(d[:, 0]))
    return out


wb = load_wave("ANS_CVS_waveform_before.dat")
wa = load_wave("ANS_CVS_waveform_after.dat")
ab = load_act("ANS_activity_before.dat")
aa = load_act("ANS_activity_after.dat")

T_ONSET = 20.0

# ----------------------------------------------------------------------
# Figure 1: before vs after, zoomed on the onset region
# ----------------------------------------------------------------------
fig, ax = plt.subplots(3, 2, figsize=(14, 10), sharex=True)
lo, hi = 17.0, 30.0


def win(d):
    m = (d["time"] >= lo) & (d["time"] <= hi)
    return m


for col, (wd, ad, label) in enumerate(
        [(wb, ab, "BEFORE fix"), (wa, aa, "AFTER fix")]):
    mw = win(wd)
    ma = win(ad)

    # Row 0: arterial pressure waveform
    ax[0, col].plot(wd["time"][mw], wd["P_art"][mw], lw=0.8,
                    color="steelblue")
    ax[0, col].axvline(T_ONSET, ls="--", color="red", alpha=0.7)
    ax[0, col].set_title(f"{label}: arterial pressure")
    ax[0, col].set_ylabel("P_art (mmHg)")
    ax[0, col].grid(alpha=0.3)

    # Row 1: HR and R_vb
    ax[1, col].plot(ad["time"][ma], ad["HR"][ma], color="navy", label="HR")
    ax[1, col].axvline(T_ONSET, ls="--", color="red", alpha=0.7)
    ax[1, col].set_ylabel("HR (beats/s)")
    ax1b = ax[1, col].twinx()
    ax1b.plot(ad["time"][ma], ad["R_vb"][ma], color="darkgreen",
              label="R_vb")
    ax1b.plot(ad["time"][ma], ad["R_eq"][ma], color="orange", ls="--",
              label="R_eq")
    ax1b.set_ylabel("R (mmHg·s/mL)")
    ax[1, col].grid(alpha=0.3)
    lines1, labs1 = ax[1, col].get_legend_handles_labels()
    lines2, labs2 = ax1b.get_legend_handles_labels()
    ax[1, col].legend(lines1 + lines2, labs1 + labs2,
                      loc="upper right", fontsize=8)

    # Row 2: total inflow (shows the spurious inflow pulse before the fix)
    ax[2, col].plot(wd["time"][mw], wd["Q_in"][mw], lw=0.7,
                    color="purple")
    ax[2, col].axvline(T_ONSET, ls="--", color="red", alpha=0.7)
    ax[2, col].set_ylabel("Q_in_total (mL/s)")
    ax[2, col].set_xlabel("Time (s)")
    ax[2, col].grid(alpha=0.3)

# Shared y-limits per row for honest comparison
for r in range(3):
    lims = [ax[r, c].get_ylim() for c in range(2)]
    lo_y = min(l[0] for l in lims); hi_y = max(l[1] for l in lims)
    for c in range(2):
        ax[r, c].set_ylim(lo_y, hi_y)

fig.suptitle("ANS-onset BP spike: before vs after the fix "
             "(red dashed = ANS onset at t = 20 s)", fontsize=13)
plt.tight_layout()
plt.savefig("fig_spike_fix_compare.png", dpi=130)
print("Saved fig_spike_fix_compare.png")


# ----------------------------------------------------------------------
# Figure 2: after-fix mechanism detail, full timeline
# ----------------------------------------------------------------------
fig2, ax2 = plt.subplots(5, 1, figsize=(12, 13), sharex=True)
t = wa["time"]

ax2[0].plot(t, wa["P_art"], lw=0.4, color="steelblue")
ax2[0].set_ylabel("P_art (mmHg)")
ax2[0].set_title("After fix: smooth ANS engagement, no onset spike")

ax2[1].plot(aa["time"], aa["HR"], color="navy")
ax2[1].set_ylabel("HR (beats/s)")

ax2[2].plot(aa["time"], aa["R_vb"], color="black", label="R_vb")
ax2[2].plot(aa["time"], aa["R_eq"], color="orange", ls="--", label="R_eq")
ax2[2].set_ylabel("R (mmHg·s/mL)")
ax2[2].legend(loc="upper right", fontsize=8)

ax2[3].plot(t, wa["Q_fwd"], lw=0.4, color="seagreen", label="Q_forward")
ax2[3].plot(t, wa["Q_ref"], lw=0.4, color="purple",   label="Q_reflected")
ax2[3].set_ylabel("Flow (mL/s)")
ax2[3].legend(loc="upper right", fontsize=8)

ax2[4].plot(aa["time"], aa["ans_gain"], color="crimson", lw=1.5,
            label="ans_gain")
ax2[4].plot(aa["time"], aa["ns_bp"], color="firebrick", ls=":",
            label="ns_bp")
ax2[4].plot(aa["time"], aa["np_bp"], color="royalblue", ls=":",
            label="np_bp")
ax2[4].set_ylabel("gain / activity")
ax2[4].set_xlabel("Time (s)")
ax2[4].legend(loc="upper right", fontsize=8)

for a in ax2:
    a.axvline(T_ONSET, ls="--", color="red", alpha=0.6)
    a.axvline(24.0, ls=":", color="gray", alpha=0.6)
    a.grid(alpha=0.2)
ax2[4].text(20.0, 0.05, " onset", color="red", fontsize=8)
ax2[4].text(24.0, 0.05, " ramp end", color="gray", fontsize=8)

plt.tight_layout()
plt.savefig("fig_spike_fix_detail.png", dpi=130)
print("Saved fig_spike_fix_detail.png")
