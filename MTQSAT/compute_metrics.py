"""Compute the thesis results table (tab:metricas) from FswComparison.csv.

Usage: python MTQSAT/compute_metrics.py [CSV_OR_DIRECTORY]
Prints each metric and the LaTeX table rows ready to paste.
"""

import csv
import sys
from pathlib import Path

import numpy as np

ORBIT_PERIOD = 5640.0          # s (~10 orbits in 56400 s)
NADIR_COMMAND_TIME = 28200.0   # s
W_LOW = 0.005                  # rad/s, detumble_w_thld_low
RMM_TRUE_FALLBACK = np.array([0.001, -0.001, 0.001])  # A m^2, SC_MTQSAT.txt


def load(path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        rows = [row for row in reader if len(row) == len(header)]
    values = np.asarray(rows, dtype=float)
    return dict(zip(header, values.T))


def vector(data, prefix, count=3):
    names = [f"{prefix}_{component}" for component in range(count)]
    return np.column_stack([data[name] for name in names]) if all(name in data for name in names) else None


def angle_deg(first, second, valid):
    denominator = np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)
    usable = valid & np.isfinite(denominator) & (denominator > 1e-12)
    result = np.full(len(first), np.nan)
    cosine = np.sum(first[usable] * second[usable], axis=1) / denominator[usable]
    result[usable] = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    return result


def stats(values, which):
    values = values[np.isfinite(values)]
    if values.size == 0:
        return "sin datos"
    parts = {"media": np.mean, "mediana": np.median, "p95": lambda v: np.percentile(v, 95), "max": np.max}
    return " / ".join(f"{parts[name](values):.2f}" for name in which)


def main():
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent
    path = base / "FswComparison.csv" if base.is_dir() else base
    data = load(path)

    t = data["time"]
    eclipse = data["truth_eclipse"].astype(bool)
    mode = data["fsw_mode"]
    end = t[-1]

    # Windows
    sun_window = (t >= NADIR_COMMAND_TIME - 2 * ORBIT_PERIOD) & (t < NADIR_COMMAND_TIME)
    nadir_window = t >= end - 3 * ORBIT_PERIOD

    # 1. Automatic transition to sun pointing
    sun_idx = np.argmax(mode == 2) if np.any(mode == 2) else None
    row1 = f"$t = {t[sun_idx]:.0f}$~s" if sun_idx is not None else "no ocurre"

    # 2. True ||w|| below w_low
    w_true = vector(data, "truth_wbn")
    w_norm = np.linalg.norm(w_true, axis=1)
    below = np.argmax(w_norm < W_LOW) if np.any(w_norm < W_LOW) else None
    row2 = f"$t = {t[below]:.0f}$~s" if below is not None else "no ocurre"

    # 3. Sun estimation error (mean / p95), sun valid & sunlit
    sun_err = angle_deg(vector(data, "truth_svb"), vector(data, "fsw_svb"),
                        ~eclipse & data["fsw_sun_valid"].astype(bool))
    row3 = stats(sun_err, ("media", "p95")) + "$^\\circ$"

    # 4. Sun pointing error, 2 orbits before nadir, sunlit (mean/median/p95)
    sun_target = vector(data, "sun_target")
    point_err = angle_deg(vector(data, "truth_svb"), sun_target, ~eclipse)
    row4 = stats(point_err[sun_window], ("media", "mediana", "p95")) + "$^\\circ$"

    # 5. Nadir knowledge error with valid horizon (mean / p95)
    horizon = data.get("fsw_horizon_valid", np.ones(len(t))).astype(bool)
    nadir_know = data["nadir_err_deg"].copy()
    nadir_know[~horizon] = np.nan
    row5 = stats(nadir_know, ("media", "p95")) + "$^\\circ$"

    # 6. Nadir pointing error, last 3 orbits (mean/median/p95)
    truth_nadir = vector(data, "truth_nadir_b")
    nadir_target = np.tile((0.0, 0.0, 1.0), (len(t), 1))
    nadir_point = angle_deg(truth_nadir, nadir_target, np.ones(len(t), dtype=bool))
    row6 = stats(nadir_point[nadir_window], ("media", "mediana", "p95")) + "$^\\circ$"

    # 7-9. RMM relative error
    rmm_true = vector(data, "truth_residual_dipole")
    if rmm_true is None:
        rmm_true = np.tile(RMM_TRUE_FALLBACK, (len(t), 1))
    state = vector(data, "ekf_state", 6)
    valid = data["ekf_valid"].astype(bool)
    rmm_err = np.where(valid,
                       np.linalg.norm(state[:, 3:] - rmm_true, axis=1)
                       / np.linalg.norm(rmm_true, axis=1) * 100, np.nan)
    crossing = np.argmax(np.nan_to_num(rmm_err, nan=np.inf) < 10)
    row7 = f"$t = {t[crossing]:.0f}$~s" if rmm_err[crossing] < 10 else "no baja del 10\\%"
    row8 = stats(rmm_err[sun_window], ("media", "max")) + "\\%"
    row9 = stats(rmm_err[nadir_window], ("media", "max")) + "\\%"

    # 10. Maximum applied dipole
    actual = [name for name in data if "mtb" in name and name.endswith("_actual")]
    dipole_max = max(np.abs(data[name]).max() for name in actual) if actual else np.nan
    row10 = f"${dipole_max:.3f}~\\mathrm{{A\\,m^2}}$"

    rows = (
        ("Transición automática a apuntamiento al Sol", row1),
        ("Norma verdadera de $\\omega$ menor que $\\omega_{low}$", row2),
        ("Error de estimación del Sol (media / percentil 95)", row3),
        ("Error de apuntamiento al Sol, 2 órbitas previas a nadir, fuera de eclipse (media / mediana / percentil 95)", row4),
        ("Error de conocimiento del nadir con horizonte válido (media / percentil 95)", row5),
        ("Error de apuntamiento al nadir, últimas 3 órbitas (media / mediana / percentil 95)", row6),
        ("Instante en que el error de $\\hat{\\vec{m}}_r$ baja del $10\\%$", row7),
        ("Error de $\\hat{\\vec{m}}_r$, 2 órbitas previas a nadir (media / máximo)", row8),
        ("Error de $\\hat{\\vec{m}}_r$, últimas 3 órbitas (media / máximo)", row9),
        ("Dipolo máximo aplicado por los magnetorquers", row10),
    )
    print("% Filas generadas por compute_metrics.py — pegar en tab:metricas")
    for name, value in rows:
        print(f"{name} & {value} \\\\ \\hline")


if __name__ == "__main__":
    main()
