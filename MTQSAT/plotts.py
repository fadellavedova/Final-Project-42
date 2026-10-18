"""Python port of plotts.m. Usage: python plotts.py [data_dir]"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent


def load(name):
    # Drop rows with a different column count (e.g. a partially flushed last line)
    rows = [ln.split() for ln in (data_dir / name).read_text().splitlines() if ln.strip()]
    ncol = len(rows[0])
    return np.array([[float(v) for v in r] for r in rows if len(r) == ncol])


MTB = load("MTB.42")
WBN = load("wbn.42")
QBN = load("qbn.42")
POSN = load("PosN.42")
ILLUM = load("Illum.42")
SVB = load("svb.42")

# Output files can end up with different lengths if the run was cut short
N = min(len(a) for a in (MTB, WBN, QBN, POSN, ILLUM, SVB))
MTB, WBN, QBN, POSN, ILLUM, SVB = (a[:N] for a in (MTB, WBN, QBN, POSN, ILLUM, SVB))

ORBIT_PERIOD = 96.77 * 60
orbit = np.arange(1, len(POSN) + 1) / ORBIT_PERIOD

# Figure 2: POSN z
plt.figure(2)
plt.plot(orbit, POSN[:, 2])
plt.grid(True)
plt.title("POSN")
plt.xlabel("Orbits")

# Figure 1: MTB, WBN, QBN
fig, ax = plt.subplots(3, 1, num=1, clear=True)
for i, c in enumerate("rgb"):
    ax[0].plot(orbit, MTB[:, i], c)
ax[0].set(title="MTB", xlabel="Orbits", ylabel="Am^2", xlim=(0, 10), ylim=(-0.03, 0.03))
for i, c in enumerate("rgb"):
    ax[1].plot(orbit, WBN[:, i], c)
ax[1].set(title="WBN", xlabel="Orbits", ylabel="rad/sec", xlim=(0, 10), ylim=(-0.003, 0.003))
for i, c in enumerate("rgbk"):
    ax[2].plot(orbit, QBN[:, i], c)
ax[2].set(title="QBN", xlabel="Orbits", xlim=(0, 10), ylim=(-1, 1))
for a in ax:
    a.grid(True)
fig.tight_layout()

# Figure 3: QBN(0:2) [deg]
plt.figure(3)
for i, c in enumerate("rgb"):
    plt.plot(orbit, QBN[:, i] * 180.0 / np.pi, c)
plt.grid(True)
plt.title("QBN(0:2) [deg]")
plt.xlabel("Orbits")

# Figure 4: SVB while illuminated
iiv = np.where(ILLUM.max(axis=1) > 0.0)[0]
plt.figure(4)
for i, c in enumerate("rgb"):
    plt.plot(orbit[iiv], SVB[iiv, i], c + ".")
plt.grid(True)
plt.title("SVB")
plt.xlabel("Orbits")

# Figure 41: SVB in 3D
fig41 = plt.figure(41)
ax41 = fig41.add_subplot(projection="3d")
ax41.plot(SVB[iiv, 0], SVB[iiv, 1], SVB[iiv, 2], "k.")

# Figure 5: ILLUM
plt.figure(5)
styles = ["r", "r--", "g", "g--", "b", "b--"]
for i, s in enumerate(styles[: ILLUM.shape[1]]):
    plt.plot(orbit, ILLUM[:, i], s)

# Figure 6: MTB, wheel speed, WBN, QBN
hwhl_path = data_dir / "Hwhl.42"
if hwhl_path.exists() and hwhl_path.stat().st_size > 0:
    HWHL = load("Hwhl.42")[:N]
    fig, ax = plt.subplots(4, 1, num=6, clear=True)
    for i, c in enumerate("rgb"):
        ax[0].plot(orbit, MTB[:, i], c)
    ax[0].set(title="MTB", xlabel="Orbits", ylabel="Am^2")
    ax[1].plot(orbit, HWHL[:, 0] / 2.1128e-6 / (2 * np.pi) * 60, "b")
    ax[1].set(title="Wheel Speed", xlabel="Orbits", ylabel="RPM")
    for i, c in enumerate("rgb"):
        ax[2].plot(orbit, WBN[:, i], c)
    ax[2].set(title="WBN", xlabel="Orbits", ylabel="rad/sec")
    for i, c in enumerate("rgbk"):
        ax[3].plot(orbit, QBN[:, i], c)
    ax[3].set(title="QBN", xlabel="Orbits")
    for a in ax:
        a.grid(True)
    fig.tight_layout()

# Figure 7: moving average of |cos(sun angle)|
WF = 6000
ctita = np.zeros(len(SVB))
ctita[iiv] = SVB[iiv, 1]
cossunavg = np.convolve(np.abs(ctita), np.ones(WF), mode="full") / WF
cossunavg = cossunavg[WF - 1 : len(cossunavg) - WF]
plt.figure(7)
plt.plot(np.arange(1, len(cossunavg) + 1) / ORBIT_PERIOD, cossunavg, "b--")
plt.grid(True)
plt.title("Cos(TitaSol)")
plt.xlabel("Orbits")
if len(cossunavg) > 9999:
    print(np.mean(cossunavg[9999:]))

# Figure 10: mean power generation factor vs RAAN
RAAN = [90, 67.5, 45, 35, 30, 22.5, 10, 0]
POWFACTOR = [0.6167, 0.6019, 0.6407, 0.7, 0.82614, 0.8235, 0.8234, 0.8]
POWFACTORm = [0.46, 0.45, 0.474, 0.57, 0.65, 0.63, 0.65, 0.64]
RAANw = [90, 67.5, 45, 35, 30, 22.5, 10, 0]
POWFACTORw = [0.6367, 0.661, 0.7452, 0.8395, 0.9938, 0.99864, 0.998936, 0.9983]
POWFACTORwm = [0.614, 0.628, 0.71, 0.7990, 0.96, 0.996, 0.9965, 0.985]
plt.figure(10)
plt.plot(RAAN, POWFACTOR, "ro-")
plt.plot(RAAN, POWFACTORm, "r-.")
plt.plot(RAANw, POWFACTORw, "bo-")
plt.plot(RAANw, POWFACTORwm, "b-.")
plt.plot(RAANw[-3:], POWFACTORw[-3:], "b-")
plt.grid(True)
plt.title("Mean Power Generation Factor")
plt.xlabel("RAAN")

plt.show()
