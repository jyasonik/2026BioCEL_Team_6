"""
Exercise 2.1(b) — Repressor characterization analysis

Reads CSV files output by repressor_output.gro (one file per dox value) and
generates three figures:
  Plot A: Time-course — repressor vs time (all dox values overlaid)
  Plot B: Time-course — rfp vs time (all dox values overlaid)
  Plot C: Transfer curve — rfp vs repressor at multiple time points

Usage:
  1. Run the Gro simulation for each dox value you want to compare.
     Edit line 40 of repressor_output.gro (dox := ...) and re-run each time.
     Suggested values: 1e-5, 1e-4, 1e-3, 0.01, 0.1, 1.0
  2. Place all output CSV files in the same folder and point OUTPUT_DIR below
     to that folder.
  3. Run:  python analyze_repressor.py
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

# ── Configuration ──────────────────────────────────────────────────────────────

# Directory where Gro wrote the CSV files (change to match your output_path)
OUTPUT_DIR = "/Users/jyasonik/Downloads/Gro_Mac/output/rep1/"

# Time points (in minutes) at which to sample the transfer curve (Plot C)
TRANSFER_CURVE_TIMES = [25, 50, 75, 100]

# ── Load data ──────────────────────────────────────────────────────────────────

csv_files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "rfp_repressor_dox=*.csv")))

if not csv_files:
    raise FileNotFoundError(
        f"No CSV files found in {OUTPUT_DIR!r}.\n"
        "Run the Gro simulation first and check that OUTPUT_DIR points to the right folder."
    )

datasets = {}  # dox_value -> DataFrame

# Gro's fprint ends each row with a trailing comma, producing an 11th empty
# field. Adding "_drop" accounts for it so the real columns align correctly.
column_names = [
    "cell_id", "time", "volume", "dox",
    "bfp", "repressor", "rfp",
    "plasmid_constitutive", "plasmid_repressor", "plasmid_output", "_drop"
]

for path in csv_files:
    df = pd.read_csv(path, header=0, names=column_names, skipinitialspace=True)
    df = df.drop(columns=["_drop"], errors="ignore")
    # Strip whitespace from any string columns
    df = df.apply(lambda col: col.str.strip() if col.dtype == object else col)
    df = df.astype(float)
    dox_val = df["dox"].iloc[0]
    datasets[dox_val] = df

dox_values = sorted(datasets.keys())
print(f"Loaded {len(dox_values)} simulation(s): dox = {dox_values}")

# Find the last time point where any dataset still has plasmids (meaningful data).
# After plasmids dilute to zero from repeated division, all proteins follow to zero.
max_meaningful_time = max(
    df.loc[df["plasmid_constitutive"] > 0, "time"].max()
    for df in datasets.values()
)

# Color map — one color per dox value
colors = cm.viridis(np.linspace(0.15, 0.9, len(dox_values)))

# ── Plot A: Time-course — repressor vs time ────────────────────────────────────

fig_a, ax_a = plt.subplots(figsize=(7, 4))

for color, dox_val in zip(colors, dox_values):
    df = datasets[dox_val]
    ax_a.plot(df["time"], df["repressor"], color=color, label=f"dox = {dox_val:.0e}")

ax_a.set_xlabel("Time (min)")
ax_a.set_ylabel("Repressor (molecules)")
ax_a.set_title("Time-course: Repressor vs Time")
ax_a.set_xlim(0, max_meaningful_time)
ax_a.legend(title="dox", fontsize=8, title_fontsize=8)
ax_a.grid(True, alpha=0.3)
fig_a.tight_layout()
fig_a.savefig(os.path.join(OUTPUT_DIR, "plot_A_repressor_timecourse.pdf"))
print("Saved Plot A: repressor time-course")

# ── Plot B: Time-course — rfp vs time ─────────────────────────────────────────

fig_b, ax_b = plt.subplots(figsize=(7, 4))

for color, dox_val in zip(colors, dox_values):
    df = datasets[dox_val]
    ax_b.plot(df["time"], df["rfp"], color=color, label=f"dox = {dox_val:.0e}")

ax_b.set_xlabel("Time (min)")
ax_b.set_ylabel("RFP (molecules)")
ax_b.set_title("Time-course: RFP vs Time")
ax_b.set_xlim(0, max_meaningful_time)
ax_b.legend(title="dox", fontsize=8, title_fontsize=8)
ax_b.grid(True, alpha=0.3)
fig_b.tight_layout()
fig_b.savefig(os.path.join(OUTPUT_DIR, "plot_B_rfp_timecourse.pdf"))
print("Saved Plot B: RFP time-course")

# ── Plot C: Transfer curve — rfp vs repressor at multiple time points ──────────
#
# For each time point, take the snapshot (row closest to that time) from every
# dox condition to get one (repressor, rfp) point. Each time point becomes a
# curve across dox conditions.

fig_c, ax_c = plt.subplots(figsize=(7, 4))
time_colors = cm.plasma(np.linspace(0.15, 0.9, len(TRANSFER_CURVE_TIMES)))

for t_color, t in zip(time_colors, TRANSFER_CURVE_TIMES):
    rep_vals = []
    rfp_vals = []
    for dox_val in dox_values:
        df = datasets[dox_val]
        # Find the row closest to this time point
        idx = (df["time"] - t).abs().idxmin()
        rep_vals.append(df.loc[idx, "repressor"])
        rfp_vals.append(df.loc[idx, "rfp"])
    ax_c.plot(rep_vals, rfp_vals, "o-", color=t_color, label=f"t = {t} min")

ax_c.set_xlabel("Repressor (molecules)")
ax_c.set_ylabel("RFP (molecules)")
ax_c.set_title("Transfer Curve: RFP vs Repressor")
ax_c.legend(title="Time point", fontsize=8, title_fontsize=8)
ax_c.grid(True, alpha=0.3)
fig_c.tight_layout()
fig_c.savefig(os.path.join(OUTPUT_DIR, "plot_C_transfer_curve.pdf"))
print("Saved Plot C: transfer curve")

plt.show()
print("\nDone. All plots saved to", OUTPUT_DIR)
