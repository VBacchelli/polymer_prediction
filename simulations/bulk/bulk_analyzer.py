import os
import re
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Set backend to avoid display errors on servers
plt.switch_backend("Agg")

# ==========================================
#              HELPER FUNCTIONS
# ==========================================

def check_file_existence(filepath):
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"File not found: {filepath}")


def check_dir_existence(dirpath):
    if not os.path.exists(dirpath):
        os.makedirs(dirpath, exist_ok=True)


def xvg_2_numpy(xvg_file, col1=0, col2=1):
    """Parses an XVG file and returns selected columns as numpy arrays."""
    check_file_existence(xvg_file)
    data_x = []
    data_y = []

    with open(xvg_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("@", "#")):
                continue
            parts = line.split()
            if len(parts) > max(col1, col2):
                try:
                    data_x.append(float(parts[col1]))
                    data_y.append(float(parts[col2]))
                except ValueError:
                    continue

    return np.array([data_x, data_y])


def read_xvg_matrix(xvg_file):
    """
    Reads an XVG file and returns:
      - legends: dict mapping series index -> legend string
      - data: 2D numpy array of numeric rows

    Column 0 is usually time; series s0 corresponds to column 1.
    """
    check_file_existence(xvg_file)

    legend_map = {}
    rows = []
    legend_re = re.compile(r'@\s*s(\d+)\s+legend\s+"(.+?)"')

    with open(xvg_file, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("@"):
                match = legend_re.search(line)
                if match:
                    legend_map[int(match.group(1))] = match.group(2).strip()
                continue
            if line.startswith("#"):
                continue

            parts = line.split()
            if not parts:
                continue

            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue

    if not rows:
        return legend_map, np.empty((0, 0))

    n_cols = len(rows[0])
    rows = [r for r in rows if len(r) == n_cols]
    if not rows:
        return legend_map, np.empty((0, 0))

    return legend_map, np.array(rows, dtype=float)


def run_gmx_command(cmd, stdin_text=None, cwd=None):
    """Runs a GROMACS command, capturing stdout and stderr."""
    return subprocess.run(
        cmd,
        input=stdin_text,
        text=True,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def extract_floats_from_text(text):
    """Extracts normal and scientific-notation floats from text."""
    float_re = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")
    return [float(x) for x in float_re.findall(text)]


def first_meaningful_float_from_line(line):
    """
    Heuristic parser for a GROMACS text line.
    If the line says average/mean, use the first float; otherwise use the last.
    """
    vals = extract_floats_from_text(line)
    if not vals:
        return None

    lower = line.lower()
    if "average" in lower or "avg" in lower or "mean" in lower:
        return vals[0]
    return vals[-1]


def _linear_fit_r2(x, y, params):
    """Returns R² for y = params[0] * x + params[1]."""
    if len(x) < 2:
        return 0.0
    y_pred = params[0] * x + params[1]
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    if ss_tot <= 0.0:
        return 0.0
    return 1.0 - ss_res / ss_tot


def _range_or_none(arr):
    if arr is None or len(arr) == 0:
        return None
    return [float(np.min(arr)), float(np.max(arr))]


# ==========================================
#            POLYMER CHAIN PARSERS
# ==========================================

def parse_polystat_text(text):
    """
    Parses gmx polystat stdout/stderr text for:
      - radius_of_gyration
      - end_to_end_dist

    Returned values are in Angstrom.
    """
    measures = {}
    units = {}

    for line in text.splitlines():
        lower = line.lower()

        if "gyration" in lower or "radius of gyration" in lower or re.search(r"\brg\b", lower):
            val_nm = first_meaningful_float_from_line(line)
            if val_nm is not None and "radius_of_gyration" not in measures:
                measures["radius_of_gyration"] = f"{val_nm * 10.0:.2f}"
                units["radius_of_gyration"] = "Angstrom"

        if "end to end" in lower or "end-to-end" in lower or "end_to_end" in lower:
            val_nm = first_meaningful_float_from_line(line)
            if val_nm is not None and "end_to_end_dist" not in measures:
                measures["end_to_end_dist"] = f"{val_nm * 10.0:.2f}"
                units["end_to_end_dist"] = "Angstrom"

    return measures, units


def parse_polystat_xvg(xvg_file):
    """
    Fallback parser for gmx polystat XVG.
    Looks for legends containing gyration / end-to-end and averages the series.

    Returned values are in Angstrom.
    """
    measures = {}
    units = {}

    if not os.path.exists(xvg_file):
        return measures, units

    legends, data = read_xvg_matrix(xvg_file)
    if data.size == 0 or data.shape[1] < 2:
        return measures, units

    for s_idx, legend in legends.items():
        lower = legend.lower()
        col_idx = s_idx + 1
        if col_idx >= data.shape[1]:
            continue

        series = data[:, col_idx]
        if len(series) == 0:
            continue

        if "gyration" in lower or "radius of gyration" in lower or re.search(r"\brg\b", lower):
            if "radius_of_gyration" not in measures:
#  EDIT -----------------------------------------
#                measures["radius_of_gyration"] = f"{np.mean(series) * 10.0:.2f}"
#                units["radius_of_gyration"] = "Angstrom"

                mean = float(np.mean(series) * 10.0)
                std = float(np.std(series) * 10.0)

                measures["radius_of_gyration"] = f"{mean:.2f}"
                measures["radius_of_gyration_std"] = f"{std:.2f}"

                units["radius_of_gyration"] = "Angstrom"
                units["radius_of_gyration_std"] = "Angstrom"
# -----------------------------------------
        if "end to end" in lower or "end-to-end" in lower or "end_to_end" in lower:
            if "end_to_end_dist" not in measures:
#   EDIT ------------------------------
#                measures["end_to_end_dist"] = f"{np.mean(series) * 10.0:.2f}"
#                units["end_to_end_dist"] = "Angstrom"

                 mean = float(np.mean(series) * 10.0)
                 std = float(np.std(series) * 10.0)

                 measures["end_to_end_dist"] = f"{mean:.2f}"
                 measures["end_to_end_dist_std"] = f"{std:.2f}"

                 units["end_to_end_dist"] = "Angstrom"
                 units["end_to_end_dist_std"] = "Angstrom"
# -------------------------------------------
    return measures, units


def gmx_polystat_2_chain_info(stdout: str):
    """Backward-compatible wrapper for older code."""
    return parse_polystat_text(stdout)


# ==========================================
#              PLOTTING HELPERS
# ==========================================

def plot_block_avg(block_avg_data, P_threshold, img_path):
    """Plots block-averaged stress-strain data."""
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.plot(block_avg_data[:, 0] * 100, block_avg_data[:, 1], "bo", markersize=8, label="Block Avg")
    ax.plot(block_avg_data[:, 0] * 100, block_avg_data[:, 1], color="b", linestyle="--", linewidth=1.5)
    ax.axhline(y=P_threshold, color="orange", label=f"Threshold ({P_threshold} bar)")
    ax.set_xlabel("Block Averaged Strain$_{zz}$ (%)", fontsize=14)
    ax.set_ylabel("Block Averaged Stress$_{zz}$ (bar)", fontsize=14)
    ax.legend()
    plt.savefig(img_path)
    plt.close()


def plot_fit_data(full_x, full_y, params, x_fit, y_fit_data, xlabel, ylabel, img_path):
    """Plots data and a linear fit."""
    x_min, x_max = np.min(full_x), np.max(full_x)
    y_min, y_max = np.min(full_y), np.max(full_y)

    x_axis = np.linspace(x_min, x_max, 100)
    y_line = params[0] * x_axis + params[1]

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    ax.scatter(full_x, full_y, alpha=0.3, label="All Data")
    ax.scatter(x_fit, y_fit_data, color="navy", label="Fitted Region")
    ax.plot(x_axis, y_line, color="red", linewidth=2, label="Fit")

    ax.set_xlim((x_min, x_max))
    ax.set_ylim((y_min, y_max))
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.legend()
    plt.savefig(img_path)
    plt.close()
    return img_path


# ==========================================
#            PHYSICS ANALYZERS
# ==========================================

def evaluate_tg(vol_xvg, temp_xvg, plot_path, return_diagnostics=False):
    """
    Bilinear Tg estimate from Volume vs Temperature.

    Default return is backward-compatible: Tg float.
    If return_diagnostics=True, returns (Tg, diagnostics_dict).
    """
    diagnostics = {
        "method": "bilinear fit of Volume vs Temperature",
        "volume_xvg": vol_xvg,
        "temperature_xvg": temp_xvg,
        "plot": plot_path,
        "warnings": [],
    }

    vol_data = xvg_2_numpy(vol_xvg)
    temp_data = xvg_2_numpy(temp_xvg)
    if vol_data.size == 0 or temp_data.size == 0:
        diagnostics["warnings"].append("Missing or empty volume/temperature data.")
        return (0.0, diagnostics) if return_diagnostics else 0.0

    min_len = min(vol_data.shape[1], temp_data.shape[1])
    V = vol_data[1, :min_len]
    T = temp_data[1, :min_len]

    if len(T) < 10:
        diagnostics["warnings"].append("Too few points for Tg bilinear fit.")
        return (0.0, diagnostics) if return_diagnostics else 0.0

    sorted_indices = np.argsort(T)
    T = T[sorted_indices]
    V = V[sorted_indices]

    n_points = len(T)
    split_low = max(int(n_points * 0.30), 3)
    split_high = min(int(n_points * 0.70), n_points - 3)

    low_T, low_V = T[:split_low], V[:split_low]
    high_T, high_V = T[split_high:], V[split_high:]

    slope_g, c_g = np.polyfit(low_T, low_V, 1)
    slope_r, c_r = np.polyfit(high_T, high_V, 1)

    denominator = slope_g - slope_r
    if abs(denominator) < 1e-12:
        diagnostics["warnings"].append("Bilinear Tg fit failed: nearly parallel lines.")
        Tg = 0.0
    else:
        Tg = float((c_r - c_g) / denominator)

    r2_low = _linear_fit_r2(low_T, low_V, np.array([slope_g, c_g]))
    r2_high = _linear_fit_r2(high_T, high_V, np.array([slope_r, c_r]))

    diagnostics.update({
        "n_points": int(n_points),
        "temperature_range_K": _range_or_none(T),
        "low_fit_temperature_range_K": _range_or_none(low_T),
        "high_fit_temperature_range_K": _range_or_none(high_T),
        "low_fit_r2": float(r2_low),
        "high_fit_r2": float(r2_high),
        "tg_inside_temperature_range": bool(np.min(T) <= Tg <= np.max(T)) if Tg > 0 else False,
    })

    plt.figure(figsize=(8, 6))
    plt.scatter(T, V, s=5, color="gray", alpha=0.5)
    x_range = np.linspace(min(T), max(T), 100)
    plt.plot(x_range, slope_g * x_range + c_g, "b--", label="Low-T fit")
    plt.plot(x_range, slope_r * x_range + c_r, "r--", label="High-T fit")
    if Tg > 0:
        plt.axvline(x=Tg, color="k", linestyle=":", label=f"Tg = {Tg:.1f} K")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Volume")
    plt.title(f"Tg = {Tg:.1f} K" if Tg > 0 else "Tg fit failed")
    plt.legend()
    plt.savefig(plot_path)
    plt.close()

    return (Tg, diagnostics) if return_diagnostics else Tg


def evaluate_thermal_expansion(
    vol_xvg,
    temp_xvg,
    plot_path,
    return_diagnostics=False,
    temperature_bin_width=5.0,
    min_points_per_bin=3,
):
    """
    Estimates thermal expansion from an NPT temperature-ramp simulation.

    Volumetric coefficient:
        alpha_V = (1 / V) * dV/dT

    Approximate linear coefficient for isotropic bulk:
        alpha_L = alpha_V / 3

    The analysis bins the trajectory by temperature before fitting V(T). This
    avoids over-weighting any constant-temperature hold region in an annealing
    protocol.

    Returns (measures, units) by default, or (measures, units, diagnostics)
    when return_diagnostics=True.
    """
    measures = {}
    units = {}
    diagnostics = {
        "method": "temperature-binned linear fit of Volume vs Temperature; alpha_V = (1/V) dV/dT; alpha_L = alpha_V/3",
        "volume_xvg": vol_xvg,
        "temperature_xvg": temp_xvg,
        "plot": plot_path,
        "temperature_bin_width_K": float(temperature_bin_width),
        "min_points_per_bin": int(min_points_per_bin),
        "warnings": [],
    }

    vol_data = xvg_2_numpy(vol_xvg)
    temp_data = xvg_2_numpy(temp_xvg)
    if vol_data.size == 0 or temp_data.size == 0:
        diagnostics["warnings"].append("Missing or empty volume/temperature data.")
        return (measures, units, diagnostics) if return_diagnostics else (measures, units)

    min_len = min(vol_data.shape[1], temp_data.shape[1])
    V_raw = vol_data[1, :min_len]
    T_raw = temp_data[1, :min_len]

    mask = np.isfinite(V_raw) & np.isfinite(T_raw)
    V_raw = V_raw[mask]
    T_raw = T_raw[mask]

    if len(T_raw) < 10:
        diagnostics["warnings"].append("Too few points for thermal expansion fit.")
        return (measures, units, diagnostics) if return_diagnostics else (measures, units)

    # Temperature binning prevents a thermostat hold segment from dominating
    # the linear fit. One averaged point is used per populated temperature bin.
    t_min = float(np.min(T_raw))
    t_max = float(np.max(T_raw))
    if temperature_bin_width <= 0:
        temperature_bin_width = 5.0

    edges = np.arange(
        np.floor(t_min / temperature_bin_width) * temperature_bin_width,
        np.ceil(t_max / temperature_bin_width) * temperature_bin_width + temperature_bin_width,
        temperature_bin_width,
    )

    T_fit = []
    V_fit = []
    bin_counts = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = (T_raw >= lo) & (T_raw < hi)
        if np.sum(in_bin) >= min_points_per_bin:
            T_fit.append(float(np.mean(T_raw[in_bin])))
            V_fit.append(float(np.mean(V_raw[in_bin])))
            bin_counts.append(int(np.sum(in_bin)))

    T_fit = np.asarray(T_fit, dtype=float)
    V_fit = np.asarray(V_fit, dtype=float)

    if len(T_fit) < 4:
        diagnostics["warnings"].append("Too few populated temperature bins for thermal expansion fit.")
        return (measures, units, diagnostics) if return_diagnostics else (measures, units)

    order = np.argsort(T_fit)
    T_fit = T_fit[order]
    V_fit = V_fit[order]

    params = np.polyfit(T_fit, V_fit, 1)
    dVdT = float(params[0])
    V_mean = float(np.mean(V_fit))
    alpha_v = dVdT / V_mean if V_mean != 0.0 else 0.0
    alpha_l = alpha_v / 3.0
    r2 = _linear_fit_r2(T_fit, V_fit, params)

    measures["volumetric_thermal_expansion_coefficient"] = f"{alpha_v:.6e}"
    units["volumetric_thermal_expansion_coefficient"] = "1/K"
    measures["linear_thermal_expansion_coefficient_approx"] = f"{alpha_l:.6e}"
    units["linear_thermal_expansion_coefficient_approx"] = "1/K"

    diagnostics.update({
        "n_raw_points": int(len(T_raw)),
        "n_temperature_bins_used": int(len(T_fit)),
        "temperature_range_K": _range_or_none(T_fit),
        "volume_range": _range_or_none(V_fit),
        "fit_r2": float(r2),
        "dVdT": dVdT,
        "mean_volume": V_mean,
        "alpha_v_1_per_K": float(alpha_v),
        "alpha_l_approx_1_per_K": float(alpha_l),
        "assumption_for_alpha_l": "isotropic bulk: alpha_L = alpha_V / 3",
        "bin_counts": bin_counts,
    })

    x_fit = np.linspace(np.min(T_fit), np.max(T_fit), 100)
    y_fit = params[0] * x_fit + params[1]

    plt.figure(figsize=(8, 6))
    plt.scatter(T_raw, V_raw, s=3, color="lightgray", alpha=0.35, label="Raw data")
    plt.scatter(T_fit, V_fit, s=20, color="black", label="Temperature-bin averages")
    plt.plot(x_fit, y_fit, "r--", label=f"Linear fit, R² = {r2:.3f}")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Volume")
    plt.title(f"alpha_V = {alpha_v:.3e} 1/K; alpha_L ≈ {alpha_l:.3e} 1/K")
    plt.legend()
    plt.savefig(plot_path)
    plt.close()

    return (measures, units, diagnostics) if return_diagnostics else (measures, units)

def measure_young_poisson(work_dir: str, edr_file: str, start=100, P_threshold=200, n_avg=15, return_diagnostics=False):
    """
    Measures Young's modulus and Poisson's ratio from uniaxial deformation.

    Backward-compatible return:
        measures, units, young_img, poisson_img

    With return_diagnostics=True:
        measures, units, young_img, poisson_img, diagnostics
    """
    measures = {}
    units = {}
    diagnostics = {
        "method": "linear fits from uniaxial deformation: stress-strain and transverse-strain/axial-strain",
        "edr_file": edr_file,
        "start_ps": start,
        "stress_threshold_bar": P_threshold,
        "n_avg_blocks": n_avg,
        "warnings": [],
    }
    original_cwd = os.getcwd()

    young_img = ""
    poisson_img = ""

    try:
        check_dir_existence(work_dir)
        os.chdir(work_dir)
        check_file_existence(edr_file)

        data_filename = "data_stretch.xvg"
        cmd = f"echo 'Box-X\\nBox-Y\\nBox-Z\\nPres-ZZ' | gmx energy -f {edr_file} -b {start} -o {data_filename}"
        subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        raw_x = xvg_2_numpy(data_filename, col1=0, col2=1)
        raw_z = xvg_2_numpy(data_filename, col1=0, col2=3)
        raw_p = xvg_2_numpy(data_filename, col1=0, col2=4)

        if raw_x.size == 0:
            diagnostics["warnings"].append("No box/stress data extracted from EDR.")
            out = (measures, units, young_img, poisson_img, diagnostics) if return_diagnostics else (measures, units, young_img, poisson_img)
            return out

        box_x = raw_x[1]
        box_z = raw_z[1]
        pres_zz = raw_p[1]

        trans_strain_data = -(box_x - box_x[0]) / box_x[0]
        axial_strain_data = (box_z - box_z[0]) / box_z[0]
        pressure_data = -1.0 * pres_zz

        n_points_per_block = max(int(len(axial_strain_data) / n_avg), 1)
        block_avg_data = []
        for i in range(n_avg):
            s, e = i * n_points_per_block, (i + 1) * n_points_per_block
            if e > len(axial_strain_data):
                break
            block_avg_data.append([np.mean(axial_strain_data[s:e]), np.mean(pressure_data[s:e])])
        block_avg_data = np.array(block_avg_data)

        if len(block_avg_data) > 0:
            plot_block_avg(block_avg_data, P_threshold, "block_avg_stress_strain.png")
            valid_mask = block_avg_data[:, 1] < P_threshold
            if np.sum(valid_mask) == 0:
                valid_mask[0] = True
            last_valid = np.where(valid_mask)[0][-1]
            cutoff = min((last_valid + 1) * n_points_per_block, len(axial_strain_data))
        else:
            cutoff = len(axial_strain_data)
            diagnostics["warnings"].append("Block averaging produced no blocks; fitting full data.")

        diagnostics.update({
            "n_data_points": int(len(axial_strain_data)),
            "n_fit_points": int(cutoff),
            "axial_strain_range_full": _range_or_none(axial_strain_data),
            "stress_range_full_bar": _range_or_none(pressure_data),
        })

        # Young's modulus fit
        x_y = axial_strain_data[:cutoff]
        y_y = pressure_data[:cutoff]
        if len(x_y) > 5:
            params_y = np.polyfit(x_y, y_y, 1)
            r2_y = _linear_fit_r2(x_y, y_y, params_y)
            measures["young_modulus"] = f"{params_y[0] * 1e-4:.2f}"
            units["young_modulus"] = "GPa"
            diagnostics["young_modulus_fit"] = {
                "slope_bar": float(params_y[0]),
                "intercept_bar": float(params_y[1]),
                "fit_r2": float(r2_y),
                "strain_range": _range_or_none(x_y),
                "stress_range_bar": _range_or_none(y_y),
            }
            young_img = plot_fit_data(
                axial_strain_data,
                pressure_data,
                params_y,
                x_y,
                y_y,
                "Strain ZZ",
                "Stress ZZ (bar)",
                "young_modulus.png",
            )
        else:
            diagnostics["warnings"].append("Too few points for Young's modulus fit.")

        # Poisson ratio fit
        x_p = axial_strain_data[:cutoff]
        y_p = trans_strain_data[:cutoff]
        if len(x_p) > 5:
            params_p = np.polyfit(x_p, y_p, 1)
            r2_p = _linear_fit_r2(x_p, y_p, params_p)
            measures["poisson_ratio"] = f"{params_p[0]:.2f}"
            units["poisson_ratio"] = "dimensionless"
            diagnostics["poisson_ratio_fit"] = {
                "slope": float(params_p[0]),
                "intercept": float(params_p[1]),
                "fit_r2": float(r2_p),
                "axial_strain_range": _range_or_none(x_p),
                "transverse_strain_range": _range_or_none(y_p),
            }
            poisson_img = plot_fit_data(
                axial_strain_data,
                trans_strain_data,
                params_p,
                x_p,
                y_p,
                "Strain ZZ",
                "Strain Trans",
                "poisson_ratio.png",
            )
        else:
            diagnostics["warnings"].append("Too few points for Poisson ratio fit.")

    except Exception as e:
        diagnostics["warnings"].append(f"Young/Poisson failed: {e}")
        print(f"   [Error] Young/Poisson failed: {e}")
        out = (measures, units, young_img, poisson_img, diagnostics) if return_diagnostics else (measures, units, young_img, poisson_img)
        return out
    finally:
        os.chdir(original_cwd)

    return (measures, units, young_img, poisson_img, diagnostics) if return_diagnostics else (measures, units, young_img, poisson_img)


def measure_bulk_modulus(work_dir: str, edr_file: str, start=500, n_avg=20, return_diagnostics=False):
    """
    Calculates bulk modulus from a pressure-volume fit during isotropic deformation.

    Isothermal compressibility is intentionally not reported. For coating
    screening, bulk modulus is the more direct and less redundant output.

    Backward-compatible return:
        measures, units, bulk_img

    With return_diagnostics=True:
        measures, units, bulk_img, diagnostics
    """
    measures = {}
    units = {}
    diagnostics = {
        "method": "bulk modulus from linear pressure-volume fit; K = -V dP/dV",
        "edr_file": edr_file,
        "start_ps": start,
        "n_avg": n_avg,
        "warnings": [],
    }
    original_cwd = os.getcwd()
    img = os.path.join(work_dir, "bulk_modulus.png")

    try:
        check_dir_existence(work_dir)
        os.chdir(work_dir)
        check_file_existence(edr_file)

        subprocess.run(
            f"echo Volume | gmx energy -f {edr_file} -b {start} -o volume.xvg",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        subprocess.run(
            f"echo Pressure | gmx energy -f {edr_file} -b {start} -o pressure.xvg",
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )

        vol = xvg_2_numpy("volume.xvg")[1]
        pres = xvg_2_numpy("pressure.xvg")[1]

        min_len = min(len(vol), len(pres))
        if min_len < 10:
            diagnostics["warnings"].append("Too few volume/pressure points for bulk modulus fit.")
            out = (measures, units, img, diagnostics) if return_diagnostics else (measures, units, img)
            return out

        vol = vol[:min_len]
        pres = pres[:min_len]

        avg_data = []
        n_blocks = max(int(min_len / n_avg), 1)
        for i in range(n_blocks):
            s, e = i * n_avg, (i + 1) * n_avg
            if e > min_len:
                break
            avg_data.append([np.mean(vol[s:e]), np.mean(pres[s:e])])
        avg_data = np.array(avg_data)

        if len(avg_data) < 3:
            diagnostics["warnings"].append("Too few averaged points for pressure-volume fit.")
            out = (measures, units, img, diagnostics) if return_diagnostics else (measures, units, img)
            return out

        params = np.polyfit(avg_data[:, 0], avg_data[:, 1], 1)
        r2 = _linear_fit_r2(avg_data[:, 0], avg_data[:, 1], params)

        x_fit = np.linspace(np.min(avg_data[:, 0]), np.max(avg_data[:, 0]), 100)
        plt.figure(figsize=(8, 8))
        plt.scatter(avg_data[:, 0], avg_data[:, 1], label="Data")
        plt.plot(x_fit, params[0] * x_fit + params[1], "r", label="Fit")
        plt.xlabel("Volume (nm$^3$)")
        plt.ylabel("Pressure (bar)")
        plt.title(f"Bulk modulus fit, R² = {r2:.3f}")
        plt.legend()
        plt.savefig("bulk_modulus.png")
        plt.close()

        V0 = float(np.mean(vol))
        K_bar = float(-V0 * params[0])
        K_GPa = K_bar * 1e-4

        measures["bulk_modulus"] = f"{K_GPa:.2f}"
        units["bulk_modulus"] = "GPa"

        diagnostics.update({
            "plot": img,
            "n_raw_points": int(min_len),
            "n_fit_points": int(len(avg_data)),
            "fit_r2": float(r2),
            "volume_range_nm3": _range_or_none(avg_data[:, 0]),
            "pressure_range_bar": _range_or_none(avg_data[:, 1]),
            "slope_dP_dV_bar_per_nm3": float(params[0]),
            "intercept_bar": float(params[1]),
            "mean_volume_nm3": V0,
            "bulk_modulus_bar": K_bar,
            "removed_outputs": ["isothermal_compressibility"],
        })

    except Exception as e:
        diagnostics["warnings"].append(f"Bulk Modulus failed: {e}")
        print(f"   [Error] Bulk Modulus failed: {e}")
        out = (measures, units, img, diagnostics) if return_diagnostics else (measures, units, img)
        return out
    finally:
        os.chdir(original_cwd)

    return (measures, units, img, diagnostics) if return_diagnostics else (measures, units, img)

def measure_bulk_properties(work_dir: str, edr_file: str, xtc_file: str, tpr_file: str, start=10, return_diagnostics=False):
    """
    Calculates density, radius of gyration, and end-to-end distance.

    Improvements:
      - density and polystat are handled independently
      - gmx command return codes are checked
      - polystat is parsed from stdout+stderr, with XVG fallback
      - warnings are collected in diagnostics

    Backward-compatible return:
        measures, units

    With return_diagnostics=True:
        measures, units, diagnostics
    """
    measures = {}
    units = {}
    diagnostics = {
        "method": "density from EDR; chain conformation from gmx polystat",
        "edr_file": edr_file,
        "xtc_file": xtc_file,
        "tpr_file": tpr_file,
        "start_ps": start,
        "warnings": [],
    }
    original_cwd = os.getcwd()

    try:
        check_dir_existence(work_dir)
        os.chdir(work_dir)

        # -----------------------------
        # 1. Density
        # -----------------------------
        if os.path.exists(edr_file):
            density_cmd = ["gmx", "energy", "-f", edr_file, "-b", str(start), "-o", "density.xvg"]
            density_res = run_gmx_command(density_cmd, stdin_text="Density\n")

            if density_res.returncode != 0:
                msg = "gmx energy failed for Density."
                diagnostics["warnings"].append(msg)
                print(f"   [Warning] {msg}")
                if density_res.stderr.strip():
                    print(f"   [GROMACS stderr] {density_res.stderr.strip()}")
            elif os.path.exists("density.xvg"):
                try:
                    den = xvg_2_numpy("density.xvg")
                    if den.shape[1] > 0:
                        density_mean = float(np.mean(den[1]))
                        density_std = float(np.std(den[1]))

                 # EDIT ------------------
                        #measures["density"] = f"{density_mean:.2f}"
                        #units["density"] = "kg/m^3"
                        measures["density"] = f"{density_mean:.2f}"
                        measures["density_std"] = f"{density_std:.2f}"

                        units["density"] = "kg/m^3"
                        units["density_std"] = "kg/m^3"
                 # ------------------------------------
                        diagnostics["density"] = {
                            "mean": density_mean,
                            "std": density_std,
                            "n_points": int(den.shape[1]),
                            "source": os.path.join(work_dir, "density.xvg"),
                        }
                    else:
                        diagnostics["warnings"].append("density.xvg exists but contains no numeric data.")
                except Exception as e:
                    diagnostics["warnings"].append(f"Failed to parse density.xvg: {e}")
            else:
                diagnostics["warnings"].append("gmx energy completed but density.xvg was not created.")
        else:
            diagnostics["warnings"].append(f"EDR file missing for density measurement: {edr_file}")

        # -----------------------------
        # 2. Polymer statistics
        # -----------------------------
        if os.path.exists(xtc_file) and os.path.exists(tpr_file):
            polystat_cmd = [
                "gmx", "polystat",
                "-f", xtc_file,
                "-s", tpr_file,
                "-o", "poly_stats.xvg",
                "-b", str(start),
                "-v",
            ]
            polystat_res = run_gmx_command(polystat_cmd, stdin_text="System\n")

            combined_output = ""
            if polystat_res.stdout:
                combined_output += polystat_res.stdout + "\n"
            if polystat_res.stderr:
                combined_output += polystat_res.stderr + "\n"

            if polystat_res.returncode != 0:
                msg = "gmx polystat failed."
                diagnostics["warnings"].append(msg)
                print(f"   [Warning] {msg}")
                if polystat_res.stderr.strip():
                    print(f"   [GROMACS stderr] {polystat_res.stderr.strip()}")

            info_txt, units_txt = parse_polystat_text(combined_output)
            measures.update(info_txt)
            units.update(units_txt)

            need_rg = "radius_of_gyration" not in measures
            need_ete = "end_to_end_dist" not in measures

            if need_rg or need_ete:
                info_xvg, units_xvg = parse_polystat_xvg("poly_stats.xvg")
                for key, val in info_xvg.items():
                    if key not in measures:
                        measures[key] = val
                        units[key] = units_xvg.get(key, "")

            diagnostics["chain_conformation"] = {
                "source": os.path.join(work_dir, "poly_stats.xvg"),
                "radius_of_gyration_found": "radius_of_gyration" in measures,
                "end_to_end_dist_found": "end_to_end_dist" in measures,
            }

            if "radius_of_gyration" not in measures:
                diagnostics["warnings"].append("Radius of gyration could not be extracted from gmx polystat output.")
                print("   [Warning] Radius of gyration could not be extracted from gmx polystat output.")
            if "end_to_end_dist" not in measures:
                diagnostics["warnings"].append("End-to-end distance could not be extracted from gmx polystat output.")
                print("   [Warning] End-to-end distance could not be extracted from gmx polystat output.")
        else:
            if not os.path.exists(xtc_file):
                diagnostics["warnings"].append(f"XTC file missing for polystat measurement: {xtc_file}")
            if not os.path.exists(tpr_file):
                diagnostics["warnings"].append(f"TPR file missing for polystat measurement: {tpr_file}")

    except Exception as e:
        diagnostics["warnings"].append(f"General Properties failed: {e}")
        print(f"   [Error] General Properties failed: {e}")
    finally:
        os.chdir(original_cwd)

    return (measures, units, diagnostics) if return_diagnostics else (measures, units)