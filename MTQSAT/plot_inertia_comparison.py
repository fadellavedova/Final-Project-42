"""Compare two simulation runs: FSW inertia model full vs diagonal.

Usage: python MTQSAT/plot_inertia_comparison.py
Reads FswComparison_Jcompleta.csv and FswComparison_Jdiagonal.csv and writes
figuras/comparacion_inercia.{png,pdf}.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

BASE = Path(__file__).parent
RUNS = (
    ("FswComparison_Jcompleta.csv", "FSW con $J$ completa", "C0"),
    ("FswComparison_Jdiagonal.csv", "FSW con $J$ diagonal", "C1"),
)


def load(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        rows = [row for row in reader if len(row) == len(header)]
    values = np.asarray(rows, dtype=float)
    return dict(zip(header, values.T))


def vector(data, prefix, count=3):
    names = [f"{prefix}_{component}" for component in range(count)]
    return np.column_stack([data[name] for name in names])


def main():
    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 6), layout="constrained")
    eclipse_shading = None
    for filename, label, color in RUNS:
        data = load(BASE / filename)
        time = data["time"] / 60
        valid = data["ekf_valid"].astype(bool)
        state = np.where(valid[:, None], vector(data, "ekf_state", 6), np.nan)

        rmm_true = vector(data, "truth_residual_dipole")
        relative = (np.linalg.norm(state[:, 3:] - rmm_true, axis=1)
                    / np.linalg.norm(rmm_true, axis=1) * 100)
        axes[0].plot(time, relative, color=color, linewidth=1, label=label)

        rate_error = np.linalg.norm(state[:, :3] - vector(data, "truth_wbn"), axis=1)
        axes[1].plot(time, rate_error, color=color, linewidth=1, label=label)
        if eclipse_shading is None:
            eclipse_shading = (time, data["truth_eclipse"].astype(bool))

    axes[0].axhline(10, color="k", linestyle=":", linewidth=1, label="Umbral del 10 %")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Error relativo del dipolo [%]")
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Error de velocidad angular [rad/s]")
    axes[1].set_xlabel("Tiempo de simulación [min]")
    for axis in axes:
        axis.grid(True, alpha=0.25)
        axis.fill_between(*eclipse_shading[:1], 0, 1, where=eclipse_shading[1], step="post",
                          transform=axis.get_xaxis_transform(), color="0.5", alpha=0.14)
        axis.legend(fontsize="small", loc="upper right")

    output = BASE / "figuras"
    output.mkdir(exist_ok=True)
    for extension in ("png", "pdf"):
        figure.savefig(output / f"comparacion_inercia.{extension}", dpi=200)
    print(f"Figura guardada en {output / 'comparacion_inercia.png'}")

    # Summary stats for the thesis text
    print("\nError relativo del dipolo, ultimas 3 orbitas (t > 39480 s):")
    for filename, label, _ in RUNS:
        data = load(BASE / filename)
        t = data["time"]
        valid = data["ekf_valid"].astype(bool)
        state = np.where(valid[:, None], vector(data, "ekf_state", 6), np.nan)
        rmm_true = vector(data, "truth_residual_dipole")
        relative = (np.linalg.norm(state[:, 3:] - rmm_true, axis=1)
                    / np.linalg.norm(rmm_true, axis=1) * 100)
        window = t >= t[-1] - 3 * 5640
        print(f"  {label}: media={np.nanmean(relative[window]):.1f}%  "
              f"max={np.nanmax(relative[window]):.1f}%")


if __name__ == "__main__":
    main()
