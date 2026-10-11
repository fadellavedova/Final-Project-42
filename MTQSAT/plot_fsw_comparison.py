"""Compare simulator truth with FSW telemetry from FswComparison.csv.

Usage: python MTQSAT/plot_fsw_comparison.py [CSV_OR_DIRECTORY] [--save DIR]
Dependencies: numpy, matplotlib. Run the rebuilt simulator to generate the CSV.
"""

import argparse
import csv
from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np


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


def build_figures(data):
    """Use identical scales and recorded time for every left/right comparison."""
    time = data["time"] / 60
    eclipse = data["truth_eclipse"].astype(bool)
    figures = []

    def finish(figure, axes, name):
        for axis in np.asarray(axes).flat:
            axis.grid(True, alpha=0.25)
            axis.set_xlabel("Simulation time [min]")
            axis.fill_between(time, 0, 1, where=eclipse, step="post",
                              transform=axis.get_xaxis_transform(), color="0.5", alpha=0.14)
            handles, labels = axis.get_legend_handles_labels()
            if handles:
                axis.legend(fontsize="small", loc="best")
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
    parser.add_argument("--save", type=Path, help="Save PNG figures to this directory")
    parser.add_argument("--no-show", action="store_true", help="Do not open interactive windows")
    args = parser.parse_args()
    path = args.input / "FswComparison.csv" if args.input.is_dir() else args.input
    try:
        data = load_telemetry(path)
        figures = build_figures(data)
    except (OSError, ValueError, KeyError) as error:
        parser.exit(1, f"Cannot plot {path}: {error}\nRebuild and run 42 to generate FswComparison.csv.\n")
    if args.save:
        args.save.mkdir(parents=True, exist_ok=True)
        for name, figure in figures:
            filename = re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")
            figure.savefig(args.save / f"{filename}.png", dpi=150)
        print(f"Saved {len(figures)} figures to {args.save.resolve()}")
    if not args.no_show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()