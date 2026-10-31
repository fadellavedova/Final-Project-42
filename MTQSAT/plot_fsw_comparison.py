"""Compare simulator truth with FSW telemetry from FswComparison.csv.

Usage: python MTQSAT/plot_fsw_comparison.py [CSV_OR_DIRECTORY] [--save DIR] [--all] [--show]
By default only the figures for the thesis report are written (PNG and PDF) to
MTQSAT/figuras, in Spanish and without titles (the caption goes in LaTeX).
--all writes every figure, in the original English form, with titles.
Dependencies: numpy, matplotlib. Run the rebuilt simulator to generate the CSV.
"""

import argparse
import csv
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np

# (figure name, output file stem, report section), in the order they appear in the report.
REPORT_FIGURES = (
    ("Angular velocity", "velocidad_angular", "Detumbling y apuntamiento al Sol"),
    ("Attitude quaternion", "cuaternion_actitud", "Detumbling y apuntamiento al Sol"),
    ("Torque rod dipoles", "dipolos_magnetorquers", "Detumbling y apuntamiento al Sol"),
    ("Estimation and pointing errors", "errores_apuntamiento", "Detumbling y apuntamiento al Sol"),
    ("EKF angular rate", "ekf_velocidad_angular", "Estimación del momento magnético residual"),
    ("Residual magnetic dipole", "ekf_dipolo_residual", "Estimación del momento magnético residual"),
    ("EKF uncertainty", "ekf_incertidumbre", "Estimación del momento magnético residual"),
    ("Nadir error and horizon sensor validity", "error_nadir_horizonte", "Apuntamiento al nadir"),
    ("Sun vector in body frame", "vector_sol_cuerpo", "Validación de sensores"),
    ("Nadir vector in body frame", "vector_nadir_cuerpo", "Validación de sensores"),
    ("Magnetic field in body frame", "campo_magnetico_cuerpo", "Validación de sensores"),
    ("Eclipse and Sun availability", "eclipse_disponibilidad_sol", "Validación de sensores"),
)

SPANISH = {
    "Simulation time [min]": "Tiempo de simulación [min]",
    "Simulation truth": "Simulación (verdad)",
    "CSS-derived estimate": "Estimación con CSS",
    "Processed gyros": "Giróscopos procesados",
    "FSW attitude": "Actitud del FSW",
    "Processed magnetometers": "Magnetómetros procesados",
    "EHS-derived estimate": "Estimación con sensor de horizonte",
    "EKF estimate": "Estimación del EKF",
    "Unit vector": "Vector unitario",
    "Dimensionless": "Adimensional",
    "microtesla": "Campo [µT]",
    "A m^2": "Dipolo [A·m²]",
    "degrees": "Error [°]",
    "Applied dipole": "Dipolo aplicado",
    "FSW rod command": "Comando del FSW",
    "True eclipse": "Eclipse real",
    "FSW no valid Sun (not geometric eclipse)": "FSW sin Sol válido",
    "FSW SunValid": "FSW SunValid",
    "Sun estimation error": "Error de estimación del Sol",
    "True pointing error": "Error de apuntamiento real",
    "Quaternion separation (q and -q equivalent)": "Separación entre cuaterniones",
    "Nadir knowledge error (FSW vs truth)": "Error de estimación del nadir",
    "FSW horizon valid": "Horizonte válido (FSW)",
    "Angular-rate sigma [rad/s]": "σ velocidad angular [rad/s]",
    "Dipole sigma [A m^2]": "σ dipolo [A·m²]",
}
language = "en"


def tr(text):
    return SPANISH.get(text, text) if language == "es" else text


def slug(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def load_telemetry(path):
    """Read complete numeric records; do not silently truncate a running log."""
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.reader(stream)
        header = next(reader, [])
        if not header or len(set(header)) != len(header):
            raise ValueError("CSV header is empty or contains duplicate columns")
        records = []
        for line, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise ValueError(f"Incomplete CSV record on line {line}; finish the simulation first")
            try:
                records.append([float(value) for value in row])
            except ValueError as error:
                raise ValueError(f"Non-numeric CSV record on line {line}") from error
    if not records:
        raise ValueError("CSV contains no telemetry records")
    values = np.asarray(records, dtype=float)
    data = dict(zip(header, values.T))
    for name in ("time", "truth_eclipse", "fsw_sun_valid", "fsw_no_sun"):
        if name not in data:
            raise ValueError(f"Missing required column: {name}")
    if not np.all(np.isfinite(data["time"])) or np.any(np.diff(data["time"]) <= 0):
        raise ValueError("Telemetry time must be finite and strictly increasing")
    for name in ("truth_eclipse", "fsw_sun_valid", "fsw_no_sun"):
        if not np.all(np.isin(data[name], [0, 1])):
            raise ValueError(f"Invalid binary flag: {name}")
    return data


def vector(data, prefix, count=3):
    names = [f"{prefix}_{component}" for component in range(count)]
    return np.column_stack([data[name] for name in names]) if all(name in data for name in names) else None


def angle_degrees(first, second, valid):
    denominator = np.linalg.norm(first, axis=1) * np.linalg.norm(second, axis=1)
    usable = valid & np.isfinite(denominator) & (denominator > 1e-12)
    result = np.full(len(first), np.nan)
    cosine = np.sum(first[usable] * second[usable], axis=1) / denominator[usable]
    result[usable] = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    return result


def build_figures(data, titles=True):
    """Use identical scales and recorded time for every left/right comparison."""
    time = data["time"] / 60
    eclipse = data["truth_eclipse"].astype(bool)
    figures = []

    def finish(figure, axes, name):
        for axis in np.asarray(axes).flat:
            axis.grid(True, alpha=0.25)
            axis.set_xlabel(tr("Simulation time [min]"))
            axis.set_ylabel(tr(axis.get_ylabel()))
            axis.set_title(tr(axis.get_title()))
            axis.fill_between(time, 0, 1, where=eclipse, step="post",
                              transform=axis.get_xaxis_transform(), color="0.5", alpha=0.14)
            handles, labels = axis.get_legend_handles_labels()
            if handles:
                axis.legend(handles, [tr(label) for label in labels], fontsize="small", loc="best")
        if titles:
            figure.suptitle(f"{name} | gray shading: true eclipse")
        figures.append((name, figure))

    def pair(name, left, right, units, labels, right_title="FSW", scale=1):
        if left is None or right is None:
            return
        figure, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(13, 4), layout="constrained")
        for axis, values, title in zip(axes, (left, right), ("Simulation truth", right_title)):
            for component, label in enumerate(labels):
                axis.plot(time, values[:, component] * scale, label=label, linewidth=1)
            axis.set_title(title)
            axis.set_ylabel(units)
        finish(figure, axes, name)

    for prefix, name, units, count, scale, source in (
        ("svb", "Sun vector in body frame", "Unit vector", 3, 1, "CSS-derived estimate"),
        ("wbn", "Angular velocity", "rad/s", 3, 1, "Processed gyros"),
        ("qbn", "Attitude quaternion", "Dimensionless", 4, 1, "FSW attitude"),
        ("bvb", "Magnetic field in body frame", "microtesla", 3, 1e6, "Processed magnetometers"),
        ("posn", "Inertial position", "km", 3, 1e-3, "FSW ephemeris"),
        ("veln", "Inertial velocity", "km/s", 3, 1e-3, "FSW ephemeris"),
        ("nadir_b", "Nadir vector in body frame", "Unit vector", 3, 1, "EHS-derived estimate"),
    ):
        passthrough = {"wbn": "gyro", "qbn": "attitude", "posn": "ephem", "veln": "ephem"}.get(prefix)
        if passthrough and np.any(data.get(f"{passthrough}_truth_passthrough", np.zeros(len(time)))):
            source += " (truth passthrough flagged in log)"
        labels = ("qx", "qy", "qz", "qw") if count == 4 else ("X", "Y", "Z")
        pair(name, vector(data, f"truth_{prefix}", count), vector(data, f"fsw_{prefix}", count),
             units, labels, source, scale)

    figure, axes = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(13, 4), layout="constrained")
    axes[0].step(time, data["truth_eclipse"], where="post", label="True eclipse")
    axes[1].step(time, data["fsw_no_sun"], where="post", label="FSW no valid Sun (not geometric eclipse)")
    axes[1].step(time, data["fsw_sun_valid"], where="post", label="FSW SunValid", alpha=0.65)
    for axis in axes:
        axis.set_ylim(-0.1, 1.1)
        axis.set_yticks([0, 1])
    finish(figure, axes, "Eclipse and Sun availability")

    for pattern, name, left_suffix, right_suffix, units, left_title, right_title in (
        (r"css_(\d+)_illum", "CSS illumination", "ideal_direct", "fsw_illum", "Scaled illumination", "Ideal direct sunlight", "FSW sensor reading (may include albedo)"),
        (r"css_(\d+)_valid", "CSS validity", "valid", "fsw_valid", "Valid flag", "Sensor model", "FSW input"),
        (r"mtb_(\d+)_actual", "Torque rod dipoles", "actual", "command", "A m^2", "Applied dipole", "FSW rod command"),
        (r"wheel_(\d+)_truth_h", "Wheel angular momentum", "truth_h", "fsw_h", "N m s", "Simulation truth", "FSW measurement"),
        (r"wheel_(\d+)_truth_torque", "Wheel torque", "truth_torque", "tcmd", "N m", "Applied torque", "FSW command"),
    ):
        matches = sorted((int(match.group(1)), key.rsplit("_", 1)[0])
                         for key in data if (match := re.fullmatch(pattern, key)))
        if not matches:
            continue
        family = pattern.split("_")[0]
        figure, axes = plt.subplots(len(matches), 2, sharex=True, sharey="row",
                                   squeeze=False, figsize=(13, max(3, 2 * len(matches))), layout="constrained")
        for row, (index, _) in enumerate(matches):
            prefix = f"{family}_{index}"
            for side, suffix in enumerate((left_suffix, right_suffix)):
                if f"{prefix}_{suffix}" in data:
                    axes[row, side].plot(time, data[f"{prefix}_{suffix}"], label=f"{family} {index}", linewidth=1)
                axes[row, side].set_ylabel(units)
            if family == "css" and left_suffix == "ideal_direct":
                axes[row, 0].plot(time, data[f"{prefix}_illum"], label="Sensor model", alpha=0.65)
                axes[row, 0].plot(time, data[f"{prefix}_albedo"] * data[f"{prefix}_scale"], label="Albedo (scaled)", alpha=0.65)
                axes[row, 1].plot(time, data["css_threshold"] * data[f"{prefix}_scale"], label="FSW threshold", linestyle=":")
            axes[row, 0].set_title(left_title if row == 0 else "")
            axes[row, 1].set_title(right_title if row == 0 else "")
        finish(figure, axes, name)

    truth_sun, fsw_sun = vector(data, "truth_svb"), vector(data, "fsw_svb")
    if truth_sun is not None and fsw_sun is not None:
        figure, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 6), layout="constrained")
        axes[0].plot(time, angle_degrees(truth_sun, fsw_sun, ~eclipse & data["fsw_sun_valid"].astype(bool)), label="Sun estimation error")
        target = vector(data, "sun_target")
        if target is not None:
            axes[0].plot(time, angle_degrees(truth_sun, target, ~eclipse), label="True pointing error")
        axes[0].set_ylabel("degrees")
        truth_q, fsw_q = vector(data, "truth_qbn", 4), vector(data, "fsw_qbn", 4)
        if truth_q is not None and fsw_q is not None:
            denominator = np.linalg.norm(truth_q, axis=1) * np.linalg.norm(fsw_q, axis=1)
            valid = np.isfinite(denominator) & (denominator > 1e-12)
            error = np.full(len(time), np.nan)
            cosine = np.abs(np.sum(truth_q[valid] * fsw_q[valid], axis=1) / denominator[valid])
            error[valid] = 2 * np.degrees(np.arccos(np.clip(cosine, 0, 1)))
            axes[1].plot(time, error, label="Quaternion separation (q and -q equivalent)")
        axes[1].set_ylabel("degrees")
        finish(figure, axes, "Estimation and pointing errors")

    if "nadir_err_deg" in data:
        figure, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 6), layout="constrained")
        axes[0].plot(time, data["nadir_err_deg"], label="Nadir knowledge error (FSW vs truth)")
        axes[0].set_ylabel("degrees")
        if "fsw_horizon_valid" in data:
            axes[1].step(time, data["fsw_horizon_valid"], where="post", label="FSW horizon valid")
            axes[1].set_ylim(-0.1, 1.1)
            axes[1].set_yticks([0, 1])
        finish(figure, axes, "Nadir error and horizon sensor validity")

    state = vector(data, "ekf_state", 6)
    if state is not None:
        valid = data["ekf_valid"].astype(bool)
        state = np.where(valid[:, None], state, np.nan)
        pair("EKF angular rate", vector(data, "truth_wbn"), state[:, :3], "rad/s", ("X", "Y", "Z"), "EKF estimate")
        pair("Residual magnetic dipole", vector(data, "truth_residual_dipole"), state[:, 3:], "A m^2", ("X", "Y", "Z"), "EKF estimate")
        variance = vector(data, "ekf_variance", 6)
        if variance is not None:
            figure, axes = plt.subplots(2, 1, sharex=True, figsize=(12, 6), layout="constrained")
            sigma = np.sqrt(np.where(valid[:, None] & (variance >= 0), variance, np.nan))
            for component, label in enumerate(("X", "Y", "Z")):
                axes[0].plot(time, sigma[:, component], label=label)
                axes[1].plot(time, sigma[:, component + 3], label=label)
            axes[0].set_ylabel("Angular-rate sigma [rad/s]")
            axes[1].set_ylabel("Dipole sigma [A m^2]")
            finish(figure, axes, "EKF uncertainty")

    for family, suffix, truth_prefix, units in (
        ("gyro", "rate", "truth_wbn", "rad/s"),
        ("mag", "field", "truth_bvb", "tesla"),
    ):
        truth = vector(data, truth_prefix)
        indices = sorted(int(match.group(1)) for key in data
                         if (match := re.fullmatch(rf"{family}_(\d+)_{suffix}", key)))
        if truth is None or not indices:
            continue
        figure, axes = plt.subplots(len(indices), 2, sharex=True, sharey="row", squeeze=False,
                                   figsize=(13, max(3, len(indices) * 2)), layout="constrained")
        for row, index in enumerate(indices):
            sensor_axis = vector(data, f"{family}_{index}_axis")
            if sensor_axis is not None:
                axes[row, 0].plot(time, np.sum(truth * sensor_axis, axis=1), label=f"Ideal axis projection {index}")
            axes[row, 1].plot(time, data[f"{family}_{index}_{suffix}"], label=f"FSW sensor input {index}")
            for axis in axes[row]:
                axis.set_ylabel(units)
        finish(figure, axes, f"Raw {family} measurements")

    flags = [name for name in data if name.endswith("_truth_passthrough") or name in
             ("fsw_mode", "fsw_mag_valid", "fsw_eph_valid", "fsw_st_valid", "fsw_horizon_valid", "ekf_valid")]
    if flags:
        figure, axes = plt.subplots(len(flags), 1, sharex=True, squeeze=False,
                                   figsize=(12, 1.6 * len(flags)), layout="constrained")
        for axis, name in zip(axes.flat, flags):
            axis.step(time, data[name], where="post", label=name)
        finish(figure, axes, "FSW mode validity and truth passthrough")
    return figures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=Path(__file__).parent)
    parser.add_argument("--save", type=Path, default=Path(__file__).parent / "figuras",
                        help="Output directory (default: MTQSAT/figuras)")
    parser.add_argument("--all", action="store_true", help="Write every figure instead of only the report ones")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf"], choices=["png", "pdf", "svg"],
                        help="Output formats (default: png pdf)")
    parser.add_argument("--show", action="store_true", help="Also open interactive windows")
    parser.add_argument("--no-show", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    global language
    language = "en" if args.all else "es"
    path = args.input / "FswComparison.csv" if args.input.is_dir() else args.input
    try:
        data = load_telemetry(path)
        figures = build_figures(data, titles=args.all)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Cannot plot {path}: {error}\nRebuild and run 42 to generate FswComparison.csv.\n")

    by_name = dict(figures)
    if args.all:
        selected = [(slug(name), figure, None) for name, figure in figures]
    else:
        selected = [(stem, by_name[name], section) for name, stem, section in REPORT_FIGURES if name in by_name]
        for name, _, _ in REPORT_FIGURES:
            if name not in by_name:
                print(f"Sin generar (faltan columnas en el CSV): {name}")

    args.save.mkdir(parents=True, exist_ok=True)
    for stem, figure, _ in selected:
        for extension in args.formats:
            figure.savefig(args.save / f"{stem}.{extension}", dpi=200)
    print(f"{len(selected)} figuras guardadas en {args.save.resolve()}")
    if not args.all:
        current = None
        for stem, _, section in selected:
            if section != current:
                current = section
                print(f"\n{section}")
            print(f"  {stem}.{args.formats[0]}")
    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()