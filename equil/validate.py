"""Validation plots for the lipid-droplet equilibration snapshot.

Compares the minimized structure against a short equilibration frame:
 - z-density profile (PL phosphates / TRIO core / water)
 - a y-slab cross-section side-by-side
Usage:  python validate.py [before.gro] [after.gro] [out.png]
"""
import sys
import numpy as np
import MDAnalysis as mda
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

before = sys.argv[1] if len(sys.argv) > 1 else "min.gro"
after = sys.argv[2] if len(sys.argv) > 2 else "step6.1_demo.gro"
out = sys.argv[3] if len(sys.argv) > 3 else "equilibration_validation.png"

colors = {"TRIO": "#ff8c1a", "POPC": "#2e7fff", "DOPE": "#22b04a", "SAPI": "#b04ad6"}
fig, axes = plt.subplots(2, 2, figsize=(14, 10), facecolor="white")

for col, (lab, f) in enumerate([("Minimized (step6.0)", before),
                                 ("After 5 ps eq (step6.1)", after)]):
    u = mda.Universe(f)
    box = u.dimensions
    # --- top row: z-density profile ---
    ax = axes[0][col]
    bins = np.linspace(0, box[2], 100)
    for sel, c, lab2 in [("resname POPC DOPE SAPI and name P", "k", "PL phosphate"),
                         ("resname TRIO", "#ff8c1a", "TRIO core"),
                         ("resname TIP3 and name OH2", "#1e90ff", "water O")]:
        z = u.select_atoms(sel).positions[:, 2]
        if len(z):
            ax.hist(z, bins=bins, histtype="step", lw=2, density=True, color=c, label=lab2)
    ax.set_title(lab, fontweight="bold"); ax.set_xlabel("z (Å)"); ax.set_ylabel("density")
    ax.legend(fontsize=8)
    # --- bottom row: cross-section ---
    ax2 = axes[1][col]
    ymid = box[1] / 2
    for rn, c in colors.items():
        ag = u.select_atoms(f"resname {rn}")
        m = np.abs(ag.positions[:, 1] - ymid) < 12
        p = ag.positions[m]
        ax2.scatter(p[:, 0], p[:, 2], s=3, c=c, alpha=0.6, linewidths=0, label=rn)
    ax2.set_aspect("equal"); ax2.set_facecolor("#f7f7f7")
    ax2.set_xlabel("x (Å)"); ax2.set_ylabel("z (Å)")
    if col == 0:
        ax2.legend(markerscale=3, fontsize=8, loc="upper right")

fig.suptitle("Pure lipid droplet — minimization vs short equilibration",
             fontsize=14, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig(out, dpi=120, bbox_inches="tight")
print("saved", out)
