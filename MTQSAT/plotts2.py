"""Satellite log plotter with true (geometric) eclipse detection.

Usage: python plotts.py [data_dir] [--dt 1.0] [--period 5806.2] [--thr 0.0]
                        [--xmax 10] [--power]

Files read from data_dir: MTB, wbn, qbn, PosN, Illum, svb (.42), optional
Hwhl.42 and SunN.42.

Sun direction in N (needed for geometric eclipse):
  1. SunN.42 if present (best, write it from 42: the true sun unit vector in N).
  2. Otherwise it is estimated as mean(CBN^T * svb) over samples where a CSS sees
     light. That assumes the sun direction is ~constant over the run (it moves
     ~1 deg/day), and that qbn.42 is q_BN with the scalar part last (42 convention).
     The spread of the estimate is printed as a sanity check.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RE = 6.378135e6  # m, same value as the ADCS C code

# ----------------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------------
ap = argparse.ArgumentParser()
ap.add_argument("data_dir", nargs="?", default=Path(__file__).resolve().parent, type=Path)
ap.add_argument("--dt", type=float, default=1.0, help="sample period [s]")
ap.add_argument("--period", type=float, default=96.77 * 60, help="orbit period [s]")
ap.add_argument("--thr", type=float, default=0.0, help="CSS illum threshold for 'sees sun'")
ap.add_argument("--xmax", type=float, default=None, help="x limit in orbits")
ap.add_argument("--power", action="store_true", help="also plot power generation factor")
args = ap.parse_args()
data_dir = args.data_dir


def load(name):
    # Drop rows with a different column count (e.g. a partially flushed last line)
    rows = [ln.split() for ln in (data_dir / name).read_text().splitlines() if ln.strip()]
    ncol = len(rows[0])
    return np.array([[float(v) for v in r] for r in rows if len(r) == ncol])


MTB, WBN, QBN = load("MTB.42"), load("wbn.42"), load("qbn.42")
POSN, ILLUM, SVB = load("PosN.42"), load("Illum.42"), load("svb.42")
sunn_path = data_dir / "SunN.42"
SUNN = load("SunN.42") if sunn_path.exists() and sunn_path.stat().st_size > 0 else None

arrays = [MTB, WBN, QBN, POSN, ILLUM, SVB] + ([SUNN] if SUNN is not None else [])
N = min(len(a) for a in arrays)  # files can differ in length if the run was cut short
MTB, WBN, QBN, POSN, ILLUM, SVB = (a[:N] for a in (MTB, WBN, QBN, POSN, ILLUM, SVB))
if SUNN is not None:
    SUNN = SUNN[:N]

t = np.arange(N) * args.dt
orbit = t / args.period


# ----------------------------------------------------------------------------
# Eclipse computation
# ----------------------------------------------------------------------------
def unit(v):
    n = np.linalg.norm(v, axis=1, keepdims=True)
    return v / np.where(n == 0, 1, n)


def cbn_from_qbn(q):
    """42 convention: q = [x y z w] (scalar last), returns C_BN (N,3,3)."""
    q = unit(q[:, :4])
    q0, q1, q2, q3 = q.T
    C = np.empty((len(q), 3, 3))
    C[:, 0, 0] = 1 - 2 * (q1**2 + q2**2)
    C[:, 0, 1] = 2 * (q0 * q1 + q2 * q3)
    C[:, 0, 2] = 2 * (q0 * q2 - q1 * q3)
    C[:, 1, 0] = 2 * (q0 * q1 - q2 * q3)
    C[:, 1, 1] = 1 - 2 * (q2**2 + q0**2)
    C[:, 1, 2] = 2 * (q1 * q2 + q0 * q3)
    C[:, 2, 0] = 2 * (q0 * q2 + q1 * q3)
    C[:, 2, 1] = 2 * (q1 * q2 - q0 * q3)
    C[:, 2, 2] = 1 - 2 * (q0**2 + q1**2)
    return C


pos = POSN[:, :3].copy()
if np.median(np.linalg.norm(pos, axis=1)) < 1e5:  # looks like km
    pos *= 1e3

lit_css = ILLUM.max(axis=1) > args.thr  # what the CSS report

if SUNN is not None:
    sun_n = unit(SUNN[:, :3])
    sun_src = "SunN.42"
else:
    C = cbn_from_qbn(QBN)
    sun_n_all = unit(np.einsum("nji,nj->ni", C, SVB[:, :3]))  # C^T * svb
    use = lit_css & (np.linalg.norm(SVB[:, :3], axis=1) > 0.5)
    if not use.any():
        raise SystemExit("No illuminated samples to estimate the sun direction; provide SunN.42")
    s_est = unit(sun_n_all[use].mean(axis=0, keepdims=True))[0]
    spread = np.degrees(np.arccos(np.clip(sun_n_all[use] @ s_est, -1, 1)))
    print(f"Sun direction estimated from svb+qbn: {s_est}, "
          f"spread median {np.median(spread):.2f} deg, 95% {np.percentile(spread, 95):.2f} deg")
    if np.percentile(spread, 95) > 10:
        print("WARNING: large spread -> quaternion convention/frame mismatch or "
              "svb is not the true sun vector. Eclipse result is unreliable; use SunN.42.")
    sun_n = np.tile(s_est, (N, 1))
    sun_src = "estimated"

r_par = np.sum(pos * sun_n, axis=1)
r_perp = np.linalg.norm(pos - r_par[:, None] * sun_n, axis=1)
eclipse = (r_par < 0) & (r_perp < RE)  # cylindrical shadow
lit_geom = ~eclipse


def spans(mask, x):
    e = np.diff(mask.astype(int), prepend=0, append=0)
    for s, f in zip(np.where(e == 1)[0], np.where(e == -1)[0] - 1):
        yield x[s], x[f]


def shade(ax):
    for a, b in spans(eclipse, orbit):
        ax.axvspan(a, b, color="0.85", zorder=0)


def finish(axs):
    for ax in np.atleast_1d(axs):
        ax.grid(True)
        shade(ax)
        if args.xmax:
            ax.set_xlim(0, args.xmax)


# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
print(f"Samples: {N}, eclipse fraction (geometric): {eclipse.mean() * 100:.1f}%  [sun: {sun_src}]")
tp = np.sum(eclipse & ~lit_css)
fp_css = np.sum(~eclipse & ~lit_css)
fn_css = np.sum(eclipse & lit_css)
print("Geometric eclipse vs CSS (Illum.max <= thr):")
print(f"  agree in eclipse (both dark):    {tp}")
print(f"  eclipse but CSS still see light: {fn_css}  <- CSS/threshold/albedo problem")
print(f"  CSS dark but geometrically lit:  {fp_css}  <- no face sees sun, or wrong sun estimate")

# ----------------------------------------------------------------------------
# Figure 1: MTB, WBN, QBN, Illum (eclipse shaded grey)
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(4, 1, num=1, clear=True, sharex=True, figsize=(10, 9))
for i, c in enumerate("rgb"):
    ax[0].plot(orbit, MTB[:, i], c)
ax[0].set(title="MTB (grey = eclipse)", ylabel="Am^2", ylim=(-0.03, 0.03))
for i, c in enumerate("rgb"):
    ax[1].plot(orbit, WBN[:, i], c)
ax[1].set(title="WBN", ylabel="rad/s", ylim=(-0.003, 0.003))
for i, c in enumerate("rgbk"):
    ax[2].plot(orbit, QBN[:, i], c)
ax[2].set(title="QBN", ylim=(-1, 1))
for i, s in enumerate(["r", "r--", "g", "g--", "b", "b--"][: ILLUM.shape[1]]):
    ax[3].plot(orbit, ILLUM[:, i], s)
ax[3].set(title="CSS Illum", xlabel="Orbits")
finish(ax)
fig.tight_layout()

# ----------------------------------------------------------------------------
# Figure 2: eclipse diagnostic
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(4, 1, num=2, clear=True, sharex=True, figsize=(10, 8))
ax[0].step(orbit, eclipse.astype(int), "k", where="post", label="geometric eclipse")
ax[0].step(orbit, (~lit_css).astype(int) * 0.9, "r--", where="post", label="CSS dark (Illum.max <= thr)")
ax[0].set(title=f"Eclipse: geometry vs CSS (sun: {sun_src})", yticks=[0, 1])
ax[0].legend(loc="upper right")
ax[1].plot(orbit, r_perp / RE, "b")
ax[1].axhline(1.0, color="k", ls="--")
ax[1].set(ylabel="r_perp / RE (sun axis)")
ax[2].plot(orbit, ILLUM.max(axis=1), "g")
ax[2].axhline(args.thr, color="k", ls="--")
ax[2].set(ylabel="max Illum")
ax[3].plot(orbit, np.linalg.norm(SVB[:, :3], axis=1), "m")
ax[3].set(ylabel="|svb|", xlabel="Orbits")
finish(ax)
fig.tight_layout()

# ----------------------------------------------------------------------------
# Figure 3: SVB while truly sunlit; Figure 4: 3D orbit coloured by eclipse
# ----------------------------------------------------------------------------
plt.figure(3, figsize=(10, 4)).clear()
for i, c in enumerate("rgb"):
    plt.plot(orbit[lit_geom], SVB[lit_geom, i], c + ".", ms=2)
plt.grid(True)
plt.title("SVB (geometrically sunlit only)")
plt.xlabel("Orbits")

fig = plt.figure(4, figsize=(6, 6))
fig.clear()
a3 = fig.add_subplot(projection="3d")
a3.plot(*pos[lit_geom].T / 1e3, "y.", ms=1, label="sunlit")
a3.plot(*pos[eclipse].T / 1e3, "k.", ms=1, label="eclipse")
a3.quiver(0, 0, 0, *(sun_n[0] * RE * 1.5 / 1e3), color="orange")
u_, v_ = np.mgrid[0:2 * np.pi:30j, 0:np.pi:15j]
a3.plot_wireframe(RE / 1e3 * np.cos(u_) * np.sin(v_), RE / 1e3 * np.sin(u_) * np.sin(v_),
                  RE / 1e3 * np.cos(v_), color="c", alpha=0.2, lw=0.5)
a3.set(title="Orbit (km) and sun direction", box_aspect=(1, 1, 1))
a3.legend()

# ----------------------------------------------------------------------------
# Figure 5: wheel/MTB summary (if Hwhl exists)
# ----------------------------------------------------------------------------
hwhl_path = data_dir / "Hwhl.42"
if hwhl_path.exists() and hwhl_path.stat().st_size > 0:
    HWHL = load("Hwhl.42")[:N]
    fig, ax = plt.subplots(4, 1, num=5, clear=True, sharex=True, figsize=(10, 9))
    for i, c in enumerate("rgb"):
        ax[0].plot(orbit, MTB[:, i], c)
    ax[0].set(title="MTB", ylabel="Am^2")
    ax[1].plot(orbit, HWHL[:, 0] / 2.1128e-6 / (2 * np.pi) * 60, "b")
    ax[1].set(title="Wheel Speed", ylabel="RPM")
    for i, c in enumerate("rgb"):
        ax[2].plot(orbit, WBN[:, i], c)
    ax[2].set(title="WBN", ylabel="rad/s")
    for i, c in enumerate("rgbk"):
        ax[3].plot(orbit, QBN[:, i], c)
    ax[3].set(title="QBN", xlabel="Orbits")
    finish(ax)
    fig.tight_layout()

# ----------------------------------------------------------------------------
# Figure 6: moving average of |cos(sun angle)| (sunlit samples only)
# ----------------------------------------------------------------------------
WF = 6000
ctita = np.where(lit_geom, SVB[:, 1], 0.0)
cossunavg = np.convolve(np.abs(ctita), np.ones(WF), mode="full") / WF
cossunavg = cossunavg[WF - 1 : len(cossunavg) - WF]
if len(cossunavg):
    plt.figure(6, figsize=(10, 3)).clear()
    plt.plot(np.arange(1, len(cossunavg) + 1) * args.dt / args.period, cossunavg, "b--")
    plt.grid(True)
    plt.title("Cos(TitaSol) moving average")
    plt.xlabel("Orbits")
    if len(cossunavg) > 9999:
        print("mean cos(sun) after 10000 samples:", np.mean(cossunavg[9999:]))

# ----------------------------------------------------------------------------
# Optional: mean power generation factor vs RAAN
# ----------------------------------------------------------------------------
if args.power:
    RAAN = [90, 67.5, 45, 35, 30, 22.5, 10, 0]
    POWFACTOR = [0.6167, 0.6019, 0.6407, 0.7, 0.82614, 0.8235, 0.8234, 0.8]
    POWFACTORm = [0.46, 0.45, 0.474, 0.57, 0.65, 0.63, 0.65, 0.64]
    POWFACTORw = [0.6367, 0.661, 0.7452, 0.8395, 0.9938, 0.99864, 0.998936, 0.9983]
    POWFACTORwm = [0.614, 0.628, 0.71, 0.7990, 0.96, 0.996, 0.9965, 0.985]
    plt.figure(10).clear()
    plt.plot(RAAN, POWFACTOR, "ro-")
    plt.plot(RAAN, POWFACTORm, "r-.")
    plt.plot(RAAN, POWFACTORw, "bo-")
    plt.plot(RAAN, POWFACTORwm, "b-.")
    plt.plot(RAAN[-3:], POWFACTORw[-3:], "b-")
    plt.grid(True)
    plt.title("Mean Power Generation Factor")
    plt.xlabel("RAAN")

plt.show()