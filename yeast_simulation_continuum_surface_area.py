import os
import shutil
from dataclasses import dataclass, field

import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from PIL import Image, ImageOps


# =========================================================
# Data model
# =========================================================
@dataclass
class Cell:
    id: int
    parent_id: int | None
    center: np.ndarray
    radius: float
    mature_radius: float
    mature: bool = True
    has_active_bud: bool = False
    bud_id: int | None = None
    bud_axis: np.ndarray | None = None
    internal_pressure: float = 1.0
    birth_step: int | None = None  # for buds only
    used_bud_axes: list = field(default_factory=list)  # per-mother used directions

    # Gamma waiting-time scheduler (per-mother heterogeneity)
    bud_mean_interval_steps: float | None = None  # μ_i (steps)
    next_bud_step: int | None = None  # absolute step when next bud may initiate


# =========================================================
# Math helpers
# =========================================================
def unit_vector(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-12:
        return np.array([1.0, 0.0, 0.0], dtype=float)
    return v / n


def random_unit_vector(rng: np.random.Generator) -> np.ndarray:
    phi = rng.uniform(0, 2 * np.pi)
    cos_theta = rng.uniform(-1, 1)
    sin_theta = np.sqrt(1 - cos_theta ** 2)
    return np.array([sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta], dtype=float)


def volume_from_radius(r: float) -> float:
    return (4.0 / 3.0) * np.pi * r ** 3


def radius_from_volume(v: float) -> float:
    return ((3.0 * v) / (4.0 * np.pi)) ** (1.0 / 3.0)


def smooth_yield_activation(force: float, yield_force: float, sharpness: float = 6.0) -> float:
    denom = max(float(yield_force), 1e-12)
    x = (float(force) - float(yield_force)) / denom
    z = -sharpness * x
    z = float(np.clip(z, -60.0, 60.0))  # avoid overflow
    return 1.0 / (1.0 + np.exp(z))


def lognormal_params_from_mean_cv(mean: float, cv: float) -> tuple[float, float]:
    # mean = exp(mu + 0.5*sigma^2), cv^2 = exp(sigma^2) - 1
    sigma = np.sqrt(np.log(cv ** 2 + 1.0))
    mu = np.log(mean) - 0.5 * sigma ** 2
    return mu, sigma


# =========================================================
# Simulation
# =========================================================
@dataclass
class BuddingYeastPressureSimulation:
    # --- Physical / biological parameters ---
    yield_stress: float
    base_internal_pressure: float
    base_vol_growth_rate: float
    dt: float = 0.1
    target_cells: int = 50

    # --- Gamma waiting time parameters ---
    gamma_k: float = 1.1  # shared across all mothers
    mean_wait_steps: float = 20.0  # population mean (steps)
    mother_mean_cv: float = 0.30  # heterogeneity in μ_i

    # --- Growth-rate per-cell noise ---
    growth_noise_sigma: float = 0.10  # per-step multiplicative noise (std)

    # --- Geometry ---
    mature_radius: float = 1.0
    initial_bud_radius: float = 0.25
    contact_tol: float = 1e-4

    # --- Budding logic ---
    max_new_buds_per_step: int = 50
    max_direction_trials: int = 30
    bud_min_angle_deg: float = 25.0
    axis_noise: float = 0.05

    # --- Pressure-force model ---
    use_projected_area: bool = False  # False => full surface area (4πr²)

    # --- Rendering/frames ---
    out_dir: str = "output"
    frame_every: int = 15
    fig_w_px: int = 640
    fig_h_px: int = 640
    dpi: int = 100
    elev: float = 25
    azim: float = 35
    sphere_u: int = 14
    sphere_v: int = 10

    # --- Constant zoom / fixed world box ---
    dynamic_bounds: bool = False
    fixed_world_half_extent: float = 20.0
    warn_if_out_of_bounds: bool = True

    # --- Runtime ---
    rng_seed: int = 7
    step_idx: int = 0
    frame_idx: int = 0
    next_cell_id: int = 1
    cells: dict[int, Cell] = field(default_factory=dict)

    # --- Metrics tracking ---
    bud_start_times: dict[int, list] = field(default_factory=dict)
    bud_volume_records: list[tuple[int, float]] = field(default_factory=list)
    time_series: list[float] = field(default_factory=list)
    cell_counts: list[int] = field(default_factory=list)
    contact_counts: list[int] = field(default_factory=list)
    colony_volumes: list[float] = field(default_factory=list)

    def __post_init__(self):
        self.rng = np.random.default_rng(self.rng_seed)
        self._lognorm_mu, self._lognorm_sigma = lognormal_params_from_mean_cv(
            self.mean_wait_steps, self.mother_mean_cv
        )

    def initialize(self):
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)
        os.makedirs(self.out_dir, exist_ok=True)

        mother = Cell(
            id=0,
            parent_id=None,
            center=np.array([0.0, 0.0, 0.0], dtype=float),
            radius=self.mature_radius,
            mature_radius=self.mature_radius,
            mature=True,
            internal_pressure=self.base_internal_pressure
        )
        self._init_mother_waiting(mother, base_step=0)

        self.cells = {0: mother}
        self.step_idx = 0
        self.frame_idx = 0
        self.next_cell_id = 1

        self.bud_start_times = {}
        self.bud_volume_records = []
        self.time_series = []
        self.cell_counts = []
        self.contact_counts = []
        self.colony_volumes = []

    def pressure_to_force(self, pressure: float, radius: float) -> float:
        area = np.pi * radius ** 2 if self.use_projected_area else 4.0 * np.pi * radius ** 2
        return pressure * area

    def _draw_mother_mean_interval(self) -> float:
        # Per-mother heterogeneity: μ_i ~ LogNormal(mean=mean_wait_steps, CV=mother_mean_cv)
        return float(self.rng.lognormal(self._lognorm_mu, self._lognorm_sigma))

    def _draw_wait_steps(self, mean_steps: float) -> int:
        # Gamma waiting time per cycle, mean = k * theta -> theta = mean / k
        theta = max(mean_steps / max(self.gamma_k, 1e-12), 1e-6)
        wait = float(self.rng.gamma(shape=self.gamma_k, scale=theta))
        return max(1, int(np.ceil(wait)))

    def _init_mother_waiting(self, mother: Cell, base_step: int):
        mother.bud_mean_interval_steps = self._draw_mother_mean_interval()
        wait_steps = self._draw_wait_steps(mother.bud_mean_interval_steps)
        mother.next_bud_step = base_step + wait_steps

    def _schedule_next_bud(self, mother: Cell, base_step: int):
        if mother.bud_mean_interval_steps is None:
            mother.bud_mean_interval_steps = self._draw_mother_mean_interval()
        wait_steps = self._draw_wait_steps(mother.bud_mean_interval_steps)
        mother.next_bud_step = base_step + wait_steps

    def _blocked_by_contact_caps(self, mother: Cell, axis: np.ndarray) -> bool:
        # Contact inhibition: block budding into any contact cap
        for other in self.cells.values():
            if other.id == mother.id:
                continue
            dvec = other.center - mother.center
            d = np.linalg.norm(dvec)
            if d <= mother.radius + other.radius + self.contact_tol:
                if d < 1e-12:
                    return True
                nhat = dvec / d
                cos_arg = (d * d + mother.radius * mother.radius - other.radius * other.radius) / (2.0 * d * mother.radius)
                cos_arg = float(np.clip(cos_arg, -1.0, 1.0))
                alpha = np.arccos(cos_arg)
                if float(np.dot(axis, nhat)) >= np.cos(alpha):
                    return True
        return False

    def _is_valid_bud_direction(self, mother: Cell, axis: np.ndarray, bud_r: float) -> bool:
        # Reject if direction too close to previous bud axes
        min_angle_rad = np.deg2rad(self.bud_min_angle_deg)
        cos_min = np.cos(min_angle_rad)
        for prev_axis in mother.used_bud_axes:
            if np.dot(prev_axis, axis) > cos_min:
                return False

        # Contact inhibition: reject if budding into contact cap
        if self._blocked_by_contact_caps(mother, axis):
            return False

        # Reject if bud sphere would overlap any existing cell (other than mother)
        bud_center = mother.center + axis * (mother.radius + bud_r + self.contact_tol)
        for other in self.cells.values():
            if other.id == mother.id:
                continue
            d = np.linalg.norm(bud_center - other.center)
            if d < (bud_r + other.radius + self.contact_tol):
                return False

        return True

    def maybe_create_buds(self):
        if len(self.cells) >= self.target_cells:
            return

        candidates = [c for c in self.cells.values() if c.mature and not c.has_active_bud]
        self.rng.shuffle(candidates)

        n_created = 0
        for mother in candidates:
            if len(self.cells) >= self.target_cells:
                break
            if n_created >= self.max_new_buds_per_step:
                break
            if mother.next_bud_step is None or self.step_idx < mother.next_bud_step:
                continue

            bud_r = self.initial_bud_radius
            axis = None
            for _ in range(self.max_direction_trials):
                cand = unit_vector(random_unit_vector(self.rng) + self.axis_noise * random_unit_vector(self.rng))
                if self._is_valid_bud_direction(mother, cand, bud_r):
                    axis = cand
                    break
            if axis is None:
                # retry soon if no valid direction
                mother.next_bud_step = self.step_idx + 1
                continue

            bud_center = mother.center + axis * (mother.radius + bud_r + self.contact_tol)

            bud = Cell(
                id=self.next_cell_id,
                parent_id=mother.id,
                center=bud_center,
                radius=bud_r,
                mature_radius=self.mature_radius,
                mature=False,
                internal_pressure=mother.internal_pressure,
                birth_step=self.step_idx
            )

            self.cells[bud.id] = bud
            mother.has_active_bud = True
            mother.bud_id = bud.id
            mother.bud_axis = axis
            mother.used_bud_axes.append(axis)
            self.next_cell_id += 1

            t = self.step_idx * self.dt
            self.bud_start_times.setdefault(mother.id, []).append(t)
            self.bud_volume_records.append((0, volume_from_radius(bud_r)))

            n_created += 1

    def grow_buds(self):
        for c in list(self.cells.values()):
            if c.mature:
                continue

            mother = self.cells.get(c.parent_id, None)
            if mother is None:
                continue

            shared_p = 0.5 * (mother.internal_pressure + c.internal_pressure)
            mother.internal_pressure = shared_p
            c.internal_pressure = shared_p

            # Size-scaled yield resistance (force and yield force both scale with area)
            area = np.pi * c.radius ** 2 if self.use_projected_area else 4.0 * np.pi * c.radius ** 2
            growth_force = shared_p * area
            yield_force = self.yield_stress * area

            gfac = smooth_yield_activation(growth_force, yield_force, sharpness=6.0)

            # Per-cell, per-step metabolic noise
            eps = float(self.rng.normal(0.0, self.growth_noise_sigma))
            eps = float(np.clip(eps, -0.3, 0.3))
            noise_factor = max(0.0, 1.0 + eps)

            dV = self.base_vol_growth_rate * noise_factor * gfac * self.dt

            V_new = volume_from_radius(c.radius) + dV
            r_new = radius_from_volume(V_new)

            if r_new >= c.mature_radius:
                r_new = c.mature_radius
                c.mature = True
                c.internal_pressure = self.base_internal_pressure

            c.radius = r_new

            axis = unit_vector(c.center - mother.center)
            if np.linalg.norm(c.center - mother.center) < 1e-12:
                axis = random_unit_vector(self.rng)
            c.center = mother.center + axis * (mother.radius + c.radius + self.contact_tol)

            if c.mature and mother.has_active_bud and mother.bud_id == c.id:
                mother.has_active_bud = False
                mother.bud_id = None
                mother.bud_axis = None
                # Schedule mother's next bud after completion
                self._schedule_next_bud(mother, base_step=self.step_idx)

                # Initialize this newly mature cell as a mother
                if c.bud_mean_interval_steps is None:
                    self._init_mother_waiting(c, base_step=self.step_idx)

    def relax_nonbonded_overlaps(self, n_iter: int = 2):
        ids = list(self.cells.keys())
        for _ in range(n_iter):
            for i in range(len(ids)):
                ci = self.cells[ids[i]]
                for j in range(i + 1, len(ids)):
                    cj = self.cells[ids[j]]

                    if ci.parent_id == cj.id or cj.parent_id == ci.id:
                        continue

                    dvec = cj.center - ci.center
                    dist = np.linalg.norm(dvec)
                    min_dist = ci.radius + cj.radius + self.contact_tol
                    overlap = min_dist - dist

                    if overlap > 0:
                        direction = random_unit_vector(self.rng) if dist < 1e-12 else dvec / dist
                        shift = 0.5 * overlap * direction
                        ci.center -= shift
                        cj.center += shift

    def step(self):
        self.maybe_create_buds()
        self.grow_buds()
        self.relax_nonbonded_overlaps(n_iter=2)
        self.step_idx += 1

    # ---------- Metrics ----------
    def record_metrics(self):
        t = self.step_idx * self.dt
        self.time_series.append(t)
        self.cell_counts.append(len(self.cells))

        mother = self.cells.get(0, None)
        if mother is None:
            self.contact_counts.append(0)
        else:
            count = 0
            for c in self.cells.values():
                if c.id == 0:
                    continue
                d = np.linalg.norm(c.center - mother.center)
                if d <= (mother.radius + c.radius + self.contact_tol):
                    count += 1
            self.contact_counts.append(count)

        total_vol = sum(volume_from_radius(c.radius) for c in self.cells.values())
        self.colony_volumes.append(total_vol)

        for c in self.cells.values():
            if c.mature:
                continue
            if c.birth_step is None:
                continue
            age_steps = self.step_idx - c.birth_step
            vol = volume_from_radius(c.radius)
            self.bud_volume_records.append((age_steps, vol))

    # ---------- Rendering ----------
    def _plot_sphere(self, ax, center, radius, color, alpha=0.92):
        u = np.linspace(0, 2 * np.pi, self.sphere_u)
        v = np.linspace(0, np.pi, self.sphere_v)
        x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radius * np.outer(np.ones_like(u), np.cos(v))
        ax.plot_surface(x, y, z, color=color, linewidth=0, antialiased=True, shade=True, alpha=alpha)

    def _bounds(self, margin=2.0):
        if not self.dynamic_bounds:
            h = float(self.fixed_world_half_extent)
            lo = np.array([-h, -h, -h], dtype=float)
            hi = np.array([h, h, h], dtype=float)
            return lo, hi

        centers = np.array([c.center for c in self.cells.values()])
        radii = np.array([c.radius for c in self.cells.values()])
        mn = centers.min(axis=0) - radii.max() - margin
        mx = centers.max(axis=0) + radii.max() + margin
        span = np.max(mx - mn)
        mid = 0.5 * (mn + mx)
        return mid - span / 2.0, mid + span / 2.0

    def _check_bounds_warning(self):
        if self.dynamic_bounds or not self.warn_if_out_of_bounds:
            return
        h = float(self.fixed_world_half_extent)
        for c in self.cells.values():
            if np.any(c.center - c.radius < -h) or np.any(c.center + c.radius > h):
                print(
                    f"[WARN] Cell(s) exceed fixed world bounds (+/-{h}). "
                    f"Increase fixed_world_half_extent to avoid clipping."
                )
                return

    def save_frame(self):
        self._check_bounds_warning()

        fig = plt.figure(figsize=(self.fig_w_px / self.dpi, self.fig_h_px / self.dpi), dpi=self.dpi)
        ax = fig.add_subplot(111, projection="3d")

        for c in self.cells.values():
            if not c.mature:
                color = "red"
            elif c.id == 0:
                color = "blue"
            else:
                color = "green"

            self._plot_sphere(ax, c.center, c.radius, color=color)

            if c.parent_id is not None:
                p = self.cells[c.parent_id]
                ax.plot(
                    [p.center[0], c.center[0]],
                    [p.center[1], c.center[1]],
                    [p.center[2], c.center[2]],
                    color="k",
                    linewidth=0.6,
                    alpha=0.45
                )

        lo, hi = self._bounds()
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])

        ax.set_box_aspect((1, 1, 1))
        ax.view_init(elev=self.elev, azim=self.azim)
        ax.grid(False)
        ax.set_title(
            f"Yield={self.yield_stress:.3f} | BaseP={self.base_internal_pressure:.3f} | "
            f"Step={self.step_idx} | Cells={len(self.cells)}"
        )
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

        frame_path = os.path.join(self.out_dir, f"frame_{self.frame_idx:05d}.png")
        fig.savefig(frame_path, dpi=self.dpi)
        plt.close(fig)
        self.frame_idx += 1

    def run(self, max_steps: int = 250000):
        self.initialize()
        self.save_frame()
        self.record_metrics()

        while len(self.cells) < self.target_cells and self.step_idx < max_steps:
            self.step()
            self.record_metrics()

            if self.step_idx % self.frame_every == 0:
                self.save_frame()

        if self.step_idx % self.frame_every != 0:
            self.save_frame()

        return {
            "yield_stress": self.yield_stress,
            "steps": self.step_idx,
            "time": self.step_idx * self.dt,
            "cells": len(self.cells),
            "output_dir": self.out_dir,
            "times": self.time_series,
            "cell_counts": self.cell_counts,
            "contact_counts": self.contact_counts,
            "colony_volumes": self.colony_volumes,
            "bud_volume_records": self.bud_volume_records,
            "bud_start_times": self.bud_start_times,
            "dt": self.dt
        }


# =========================================================
# PNG -> MP4 utilities
# =========================================================
def padded_size(w: int, h: int, block: int = 16) -> tuple[int, int]:
    w2 = ((w + block - 1) // block) * block
    h2 = ((h + block - 1) // block) * block
    return w2, h2


def normalize_frame_to_target(path: str, target_wh: tuple[int, int]) -> np.ndarray:
    img = Image.open(path).convert("RGB")
    tw, th = target_wh

    dw = max(0, tw - img.width)
    dh = max(0, th - img.height)
    left = dw // 2
    right = dw - left
    top = dh // 2
    bottom = dh - top
    if dw > 0 or dh > 0:
        img = ImageOps.expand(img, border=(left, top, right, bottom), fill=(255, 255, 255))

    if img.width > tw or img.height > th:
        img = ImageOps.fit(img, (tw, th), method=Image.Resampling.LANCZOS, centering=(0.5, 0.5))

    return np.array(img)


def frames_to_mp4(frame_dir: str, mp4_path: str, fps: int = 12, macro_block: int = 16):
    frames = sorted([f for f in os.listdir(frame_dir) if f.lower().endswith(".png")])
    if not frames:
        print(f"[WARN] No PNG files in {frame_dir}. Skipping MP4.")
        return

    max_w, max_h = 0, 0
    for fn in frames:
        p = os.path.join(frame_dir, fn)
        with Image.open(p) as im:
            max_w = max(max_w, im.width)
            max_h = max(max_h, im.height)

    target_w, target_h = padded_size(max_w, max_h, block=macro_block)

    os.makedirs(os.path.dirname(mp4_path), exist_ok=True)
    try:
        with imageio.get_writer(
            mp4_path,
            fps=fps,
            codec="libx264",
            quality=8,
            macro_block_size=macro_block
        ) as writer:
            for fn in frames:
                p = os.path.join(frame_dir, fn)
                arr = normalize_frame_to_target(p, (target_w, target_h))
                writer.append_data(arr)
    except Exception as e:
        raise RuntimeError(
            "Failed to write MP4. Ensure ffmpeg is available (pip install imageio-ffmpeg). "
            f"Underlying error: {e}"
        )

    print(f"[MP4] Saved: {mp4_path}")


# =========================================================
# Plot helpers
# =========================================================
def plot_cellcount_overlay(summaries: list, root_output: str):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    for s in summaries:
        ax.plot(s["times"], s["cell_counts"], linewidth=2, label=f"yield={s['yield_stress']:.3f}")
    ax.set_xlabel("Simulation time (model units)")
    ax.set_ylabel("Number of cells")
    ax.set_title("Cell count vs simulation time")
    ax.grid(True, alpha=0.3)
    ax.legend()
    out_png = os.path.join(root_output, "cell_count_vs_time_overlay.png")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_contactcount_overlay(summaries: list, root_output: str):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    for s in summaries:
        ax.plot(s["times"], s["contact_counts"], linewidth=2, label=f"yield={s['yield_stress']:.3f}")
    ax.set_xlabel("Simulation time (model units)")
    ax.set_ylabel("Cells in contact with original mother")
    ax.set_title("Contacts with original mother vs time")
    ax.grid(True, alpha=0.3)
    ax.legend()
    out_png = os.path.join(root_output, "contacts_vs_time_overlay.png")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_colonyvolume_overlay(summaries: list, root_output: str):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)
    for s in summaries:
        ax.plot(s["times"], s["colony_volumes"], linewidth=2, label=f"yield={s['yield_stress']:.3f}")
    ax.set_xlabel("Simulation time (model units)")
    ax.set_ylabel("Total colony volume")
    ax.set_title("Total colony volume vs time")
    ax.grid(True, alpha=0.3)
    ax.legend()
    out_png = os.path.join(root_output, "colony_volume_vs_time_overlay.png")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_bud_volume_vs_age(summaries: list, root_output: str):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)

    for s in summaries:
        dt = s["dt"]
        age_dict = {}
        for age_steps, vol in s["bud_volume_records"]:
            age_dict.setdefault(age_steps, []).append(vol)

        ages = sorted(age_dict.keys())
        if not ages:
            continue

        times = [a * dt for a in ages]
        means = [np.mean(age_dict[a]) for a in ages]
        stds = [np.std(age_dict[a], ddof=0) for a in ages]

        ax.plot(times, means, linewidth=2, label=f"yield={s['yield_stress']:.3f}")
        ax.fill_between(times, np.array(means) - np.array(stds), np.array(means) + np.array(stds), alpha=0.2)

    ax.set_xlabel("Bud age (model time units)")
    ax.set_ylabel("Bud volume")
    ax.set_title("Bud volume vs bud age (mean ± SD)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    out_png = os.path.join(root_output, "bud_volume_vs_age_overlay.png")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_interbudding_intervals(summaries: list, root_output: str):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=140)

    for s in summaries:
        ys = s["yield_stress"]
        intervals = s["interbudding_intervals"]
        if len(intervals) == 0:
            continue

        x = np.full(len(intervals), ys, dtype=float)
        jitter = (np.random.rand(len(intervals)) - 0.5) * 0.02
        ax.scatter(x + jitter, intervals, alpha=0.6, s=18)

        mean = float(np.mean(intervals))
        std = float(np.std(intervals, ddof=0))
        ax.errorbar([ys], [mean], yerr=[std], fmt="o", capsize=4, linewidth=2, color="black")

    ax.set_xlabel("Yield stress")
    ax.set_ylabel("Inter-budding time (model units)")
    ax.set_title("Inter-budding intervals (scatter + mean ± SD)")
    ax.grid(True, alpha=0.3)
    out_png = os.path.join(root_output, "interbudding_intervals_overlay.png")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


# =========================================================
# Batch run
# =========================================================
def run_yield_stress_conditions(
    yield_stress_list,
    root_output="simulation_outputs",
    base_internal_pressure=1.0,
    base_vol_growth_rate=0.08,
    dt=0.1,
    target_cells=50,
    frame_every=15,
    mp4_fps=12,
    macro_block=16,
    fixed_world_half_extent=20.0,
    gamma_k=1.1,
    mean_wait_steps=20.0,
    mother_mean_cv=0.30,
    growth_noise_sigma=0.10
):
    os.makedirs(root_output, exist_ok=True)
    summaries = []

    for ys in yield_stress_list:
        out_dir = os.path.join(root_output, f"yield_{ys:.3f}")
        rng_seed = int(np.random.SeedSequence().entropy)

        sim = BuddingYeastPressureSimulation(
            yield_stress=ys,
            base_internal_pressure=base_internal_pressure,
            base_vol_growth_rate=base_vol_growth_rate,
            dt=dt,
            target_cells=target_cells,
            out_dir=out_dir,
            frame_every=frame_every,
            rng_seed=rng_seed,
            dynamic_bounds=False,
            fixed_world_half_extent=fixed_world_half_extent,
            fig_w_px=640,
            fig_h_px=640,
            sphere_u=14,
            sphere_v=10,
            gamma_k=gamma_k,
            mean_wait_steps=mean_wait_steps,
            mother_mean_cv=mother_mean_cv,
            growth_noise_sigma=growth_noise_sigma
        )

        s = sim.run()

        # compute inter-budding intervals
        intervals = []
        for _, times in s["bud_start_times"].items():
            if len(times) < 2:
                continue
            times_sorted = sorted(times)
            diffs = np.diff(times_sorted)
            intervals.extend(diffs.tolist())

        s["interbudding_intervals"] = intervals
        s["rng_seed"] = rng_seed
        summaries.append(s)

        print(
            f"[DONE] yield={s['yield_stress']:.3f} | cells={s['cells']} | "
            f"steps={s['steps']} | time={s['time']:.3f} | seed={rng_seed}"
        )

    # ---------- CSV outputs (long format) ----------
    with open(os.path.join(root_output, "cell_count_vs_time.csv"), "w", encoding="utf-8") as f:
        f.write("yield_stress,time,value\n")
        for s in summaries:
            for t, v in zip(s["times"], s["cell_counts"]):
                f.write(f"{s['yield_stress']},{t},{v}\n")

    with open(os.path.join(root_output, "contacts_vs_time.csv"), "w", encoding="utf-8") as f:
        f.write("yield_stress,time,value\n")
        for s in summaries:
            for t, v in zip(s["times"], s["contact_counts"]):
                f.write(f"{s['yield_stress']},{t},{v}\n")

    with open(os.path.join(root_output, "colony_volume_vs_time.csv"), "w", encoding="utf-8") as f:
        f.write("yield_stress,time,value\n")
        for s in summaries:
            for t, v in zip(s["times"], s["colony_volumes"]):
                f.write(f"{s['yield_stress']},{t},{v}\n")

    with open(os.path.join(root_output, "bud_volume_vs_age.csv"), "w", encoding="utf-8") as f:
        f.write("yield_stress,time,value\n")
        for s in summaries:
            dt = s["dt"]
            for age_steps, vol in s["bud_volume_records"]:
                t = age_steps * dt
                f.write(f"{s['yield_stress']},{t},{vol}\n")

    with open(os.path.join(root_output, "interbudding_intervals.csv"), "w", encoding="utf-8") as f:
        f.write("yield_stress,interval\n")
        for s in summaries:
            for v in s["interbudding_intervals"]:
                f.write(f"{s['yield_stress']},{v}\n")

    with open(os.path.join(root_output, "interbudding_intervals_summary.csv"), "w", encoding="utf-8") as f:
        f.write("yield_stress,mean,std,n\n")
        for s in summaries:
            vals = s["interbudding_intervals"]
            if len(vals) == 0:
                f.write(f"{s['yield_stress']},,,0\n")
            else:
                mean = float(np.mean(vals))
                std = float(np.std(vals, ddof=0))
                f.write(f"{s['yield_stress']},{mean},{std},{len(vals)}\n")

    # ---------- Plots ----------
    plot_cellcount_overlay(summaries, root_output)
    plot_contactcount_overlay(summaries, root_output)
    plot_colonyvolume_overlay(summaries, root_output)
    plot_bud_volume_vs_age(summaries, root_output)
    plot_interbudding_intervals(summaries, root_output)

    # ---------- MP4 creation ----------
    for s in summaries:
        frame_dir = s["output_dir"]
        mp4_path = os.path.join(frame_dir, f"yield_{s['yield_stress']:.3f}.mp4")
        frames_to_mp4(frame_dir, mp4_path, fps=mp4_fps, macro_block=macro_block)

    return summaries


# =========================================================
# Main
# =========================================================
if __name__ == "__main__":
    CONVERT_EXISTING_ONLY = False

    root_output = "simulation_outputs"
    mp4_fps = 12
    macro_block = 16

    if CONVERT_EXISTING_ONLY:
        yield_dirs = [
            os.path.join(root_output, d)
            for d in os.listdir(root_output)
            if os.path.isdir(os.path.join(root_output, d)) and d.startswith("yield_")
        ]
        for d in sorted(yield_dirs):
            mp4_path = os.path.join(d, f"{os.path.basename(d)}.mp4")
            frames_to_mp4(d, mp4_path, fps=mp4_fps, macro_block=macro_block)
    else:
        yield_stress_conditions = [0.05, 0.25, 0.5, 1.0, 1.5]

        base_internal_pressure = 1.0
        base_vol_growth_rate = 0.08
        dt = 0.1
        target_cells = 100
        frame_every = 15
        fixed_world_half_extent = 20.0

        run_yield_stress_conditions(
            yield_stress_list=yield_stress_conditions,
            root_output=root_output,
            base_internal_pressure=base_internal_pressure,
            base_vol_growth_rate=base_vol_growth_rate,
            dt=dt,
            target_cells=target_cells,
            frame_every=frame_every,
            mp4_fps=mp4_fps,
            macro_block=macro_block,
            fixed_world_half_extent=fixed_world_half_extent,
            gamma_k=1.1,
            mean_wait_steps=20.0,
            mother_mean_cv=0.30,
            growth_noise_sigma=0.10
        )
