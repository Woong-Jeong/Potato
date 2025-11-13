#!/usr/bin/env python3
# ------------------------------
# Module: tool_potato_peeling_path_generator.py
# Purpose:
#   - Generate helical trajectory for rotating potato with ISO-scallop adaptive pitch
#   - Robot moves in continuous spiral while potato rotates on Z-axis
#   - Curvature-based pitch adaptation for uniform surface quality
#
# Key Features:
#   - Exact curvature calculation with gradient correction
#   - ISO-scallop adaptive pitch (h = κ·s²/8)
#   - Surface following with blade engagement control
# ------------------------------

from typing import Tuple, Optional
import numpy as np
import trimesh
import open3d as o3d
from scipy.interpolate import RectBivariateSpline, splprep, splev
from scipy.ndimage import gaussian_filter


# ========================================
# CONFIGURATION PARAMETERS
# ========================================

# -------- File I/O --------
MESH_PATH = "../mesh/mesh_1_processed.stl"
OUTPUT_PATH = "../output/peeling_trajectory_adaptive.csv"

# -------- Mesh Sampling Resolution --------
NUM_Z_SLICES = 240          # Z-direction slices for cylindrical profile
NUM_THETA_SAMPLES = 360     # Angular samples (1 degree resolution)

# -------- Tool Geometry --------
TOOL_DIAMETER = 0.003       # Tool diameter: 10mm (reference only, not used for PITCH_MAX)

# -------- ISO-scallop Adaptive Pitch Parameters --------
# These parameters control the adaptive pitch based on local surface curvature
H_SCALLOP = 0.0003          # Target scallop height: 0.3mm (surface quality)
PITCH_MIN = 0.001           # Minimum pitch: 1mm (used in high curvature regions)
PITCH_MAX = 0.020           # Maximum pitch: 20mm (used in low curvature regions)
# NOTE: PITCH_MAX is independent of tool diameter in ISO-scallop method
#       Set based on acceptable surface quality and machining efficiency
#       Typical range: 1.5-3× tool diameter for smooth surfaces

# -------- Fixed Pitch Parameters (Non-adaptive Mode) --------
NUM_ROTATIONS = 150.0       # Number of full rotations during descent
NUM_SAMPLES = 10000         # Trajectory resolution

# -------- Surface Following --------
SURFACE_OFFSET = 0.001      # TCP offset from surface: 1.0mm (blade clearance)
ENGAGEMENT_MARGIN = 0.001   # Blade engagement threshold: 1mm

# -------- Smoothing --------
SPLINE_SMOOTHING = 0.5      # B-spline smoothing (0=interpolation, >0=smoothing)
CURVATURE_SMOOTHING_SIGMA = 1.5  # Gaussian filter sigma for curvature maps

# -------- Mode Selection --------
USE_ADAPTIVE = True         # True: ISO-scallop adaptive, False: fixed pitch

# -------- Adaptive Generation --------
ANGULAR_STEP = 0.01         # Angular increment: ~0.57 degrees
MAX_STEPS = 100000          # Safety limit for adaptive generation

# -------- Curvature Computation --------
CURVATURE_WINDOW = 2        # Neighborhood size for quadric fitting (5×5 window)
MIN_POINTS_FOR_FIT = 6      # Minimum points required for valid fit

# -------- Visualization --------
VIZ_MAX_POINTS = 5000       # Maximum points to show in Open3D viewer

# -------- Debug Output --------
DEBUG_CURVATURE = True      # Show detailed curvature-pitch calculations
DEBUG_FIRST_N_STEPS = 10    # Number of initial steps to show in detail
DEBUG_STEP_INTERVAL = 1000  # Show debug info every N steps

# ========================================


# -------- Mesh Loading --------
def load_mesh(mesh_path: str) -> trimesh.Trimesh:
    """Load STL mesh and compute bounds"""
    mesh = trimesh.load(mesh_path)
    print(f"Mesh loaded: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print(f"Z bounds: [{mesh.bounds[0][2]:.4f}, {mesh.bounds[1][2]:.4f}]m")
    return mesh


# -------- Radial Profile Extraction --------
def compute_cylindrical_profile(
    mesh: trimesh.Trimesh,
    num_z_slices: int = NUM_Z_SLICES,
    num_theta_samples: int = NUM_THETA_SAMPLES
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute radial profile r(θ, z) in cylindrical coordinates

    Returns:
        z_array: Z heights (increasing)
        theta_array: Angles [0, 2π)
        r_matrix: Radial distances (N_z × N_theta)
    """
    z_min, z_max = mesh.bounds[0][2], mesh.bounds[1][2]
    z_array = np.linspace(z_min, z_max, num_z_slices)
    theta_array = np.linspace(0, 2*np.pi, num_theta_samples, endpoint=False)

    r_matrix = np.zeros((num_z_slices, num_theta_samples))

    print(f"Computing cylindrical profile: {num_z_slices} Z × {num_theta_samples} θ...")

    for i, z in enumerate(z_array):
        if i % 20 == 0:
            print(f"  Slice {i+1}/{num_z_slices} (Z={z:.4f}m)")

        slice_result = mesh.section(plane_origin=[0, 0, z],
                                   plane_normal=[0, 0, 1])

        if slice_result is None:
            r_matrix[i, :] = 0
            continue

        try:
            planar, _ = slice_result.to_2D()
            vertices_2d = planar.vertices

            if len(vertices_2d) < 3:
                r_matrix[i, :] = 0
                continue

            # Compute polar coordinates
            contour_r = np.sqrt(vertices_2d[:, 0]**2 + vertices_2d[:, 1]**2)
            contour_theta = np.arctan2(vertices_2d[:, 1], vertices_2d[:, 0])

            # Sort by angle
            sorted_idx = np.argsort(contour_theta)
            contour_theta_sorted = contour_theta[sorted_idx]
            contour_r_sorted = contour_r[sorted_idx]

            # Wrap for interpolation
            contour_theta_wrapped = np.concatenate([
                contour_theta_sorted - 2*np.pi,
                contour_theta_sorted,
                contour_theta_sorted + 2*np.pi
            ])
            contour_r_wrapped = np.tile(contour_r_sorted, 3)

            # Interpolate
            r_interp = np.interp(theta_array, contour_theta_wrapped, contour_r_wrapped)
            r_matrix[i, :] = r_interp

        except Exception as e:
            print(f"  Warning: Failed at Z={z:.4f}m: {e}")
            r_matrix[i, :] = 0

    return z_array, theta_array, r_matrix


# -------- Curvature Computation Helper Functions --------
def compute_local_curvature(
    x_grid: np.ndarray,
    y_grid: np.ndarray,
    z_grid: np.ndarray,
    r_matrix: np.ndarray,
    theta_array: np.ndarray,
    i: int,
    j: int,
    window: int = CURVATURE_WINDOW
) -> Tuple[float, float]:
    """
    Compute principal curvatures at a single grid point using local quadric fitting

    Uses exact curvature formula: κ = |∂²z/∂n²| / (1 + |∇z|²)^(3/2)
    where ∇z is the gradient and ∂²z/∂n² is the directional second derivative.

    Args:
        x_grid, y_grid, z_grid: Cartesian coordinate grids
        r_matrix: Radial distance matrix
        theta_array: Angular positions
        i, j: Grid indices
        window: Neighborhood size

    Returns:
        kappa_radial: Exact curvature in radial direction (1/m)
        kappa_tangential: Exact curvature in tangential direction (1/m)
    """
    N_z, N_theta = r_matrix.shape

    # Extract local neighborhood (with boundary handling)
    i_min = max(0, i - window)
    i_max = min(N_z, i + window + 1)
    j_min = j - window
    j_max = j + window + 1

    # Handle θ periodicity
    j_indices = np.arange(j_min, j_max) % N_theta

    # Local surface patch
    x_local = x_grid[i_min:i_max, :][:, j_indices].flatten()
    y_local = y_grid[i_min:i_max, :][:, j_indices].flatten()
    z_local = z_grid[i_min:i_max, :][:, j_indices].flatten()

    # Remove invalid points (r=0)
    valid = r_matrix[i_min:i_max, :][:, j_indices].flatten() > 0
    if np.sum(valid) < MIN_POINTS_FOR_FIT:
        return 0.0, 0.0

    x_local = x_local[valid]
    y_local = y_local[valid]
    z_local = z_local[valid]

    # Fit quadric surface: z = a·x² + b·y² + c·xy + d·x + e·y + f
    try:
        A_fit = np.column_stack([
            x_local**2,
            y_local**2,
            x_local * y_local,
            x_local,
            y_local,
            np.ones_like(x_local)
        ])

        coeffs, _, _, _ = np.linalg.lstsq(A_fit, z_local, rcond=None)
        a, b, c, d, e, f = coeffs

        # Compute curvatures at center point (x_ij, y_ij)
        x_ij = x_grid[i, j]
        y_ij = y_grid[i, j]
        theta_ij = theta_array[j]

        # First derivatives (gradient)
        z_x = d  # ∂z/∂x at center (linear term)
        z_y = e  # ∂z/∂y at center (linear term)

        # Second derivatives (Hessian)
        z_xx = 2*a  # ∂²z/∂x²
        z_yy = 2*b  # ∂²z/∂y²
        z_xy = c    # ∂²z/∂x∂y

        # Radial and tangential directions
        cos_t = np.cos(theta_ij)
        sin_t = np.sin(theta_ij)

        # Directional second derivatives (numerator of curvature)
        kappa_radial = z_xx * cos_t**2 + 2*z_xy * cos_t*sin_t + z_yy * sin_t**2
        kappa_tangential = z_xx * sin_t**2 - 2*z_xy * cos_t*sin_t + z_yy * cos_t**2

        # Apply exact curvature formula: κ = |∂²z/∂n²| / (1 + |∇z|²)^(3/2)
        grad_magnitude_sq = z_x**2 + z_y**2
        denominator = (1 + grad_magnitude_sq) ** 1.5

        # Avoid division by very small numbers (nearly vertical surface)
        if denominator < 1e-10:
            return 0.0, 0.0

        kappa_radial /= denominator
        kappa_tangential /= denominator

        return abs(kappa_radial), abs(kappa_tangential)

    except np.linalg.LinAlgError:
        return 0.0, 0.0


def compute_curvature_map(
    r_matrix: np.ndarray,
    z_array: np.ndarray,
    theta_array: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute principal curvatures at each r(θ, z) grid point

    Uses local quadric surface fitting on Cartesian coordinates (x, y, z)
    with exact curvature formula: κ = |∂²z/∂n²| / (1 + |∇z|²)^(3/2)

    This accounts for surface slope (1st derivatives) to avoid overestimating
    curvature on inclined regions.

    Args:
        r_matrix: (N_z × N_theta) radial distances
        z_array: Z heights
        theta_array: Angular positions

    Returns:
        kappa_z: (N_z × N_theta) exact curvature in z-direction (1/m, for pitch adaptation)
        kappa_theta: (N_z × N_theta) exact curvature in θ-direction (1/m, feed direction)
    """
    print(f"\nComputing curvature map on r(θ, z) grid...")

    N_z, N_theta = r_matrix.shape
    kappa_z = np.zeros_like(r_matrix)
    kappa_theta = np.zeros_like(r_matrix)

    # Convert r(θ, z) grid to Cartesian (x, y, z) surface points
    theta_grid, z_grid = np.meshgrid(theta_array, z_array)
    x_grid = r_matrix * np.cos(theta_grid)
    y_grid = r_matrix * np.sin(theta_grid)

    print(f"  Processing {N_z} × {N_theta} grid points...")

    for i in range(N_z):
        if i % 30 == 0:
            print(f"    Row {i+1}/{N_z}")

        for j in range(N_theta):
            kappa_z[i, j], kappa_theta[i, j] = compute_local_curvature(
                x_grid, y_grid, z_grid, r_matrix, theta_array, i, j
            )

    # Smooth curvature map to remove noise
    kappa_z = gaussian_filter(kappa_z, sigma=CURVATURE_SMOOTHING_SIGMA)
    kappa_theta = gaussian_filter(kappa_theta, sigma=CURVATURE_SMOOTHING_SIGMA)

    print(f"  Curvature statistics:")
    print(f"    κ_z (radial):  mean={np.mean(kappa_z):.3f}, max={np.max(kappa_z):.3f}, min={np.min(kappa_z):.3f}")
    print(f"    κ_θ (tangent): mean={np.mean(kappa_theta):.3f}, max={np.max(kappa_theta):.3f}, min={np.min(kappa_theta):.3f}")

    # Show what pitch would result from these curvatures (debug mode)
    if DEBUG_CURVATURE:
        print(f"\n  Expected pitch range (from curvature map):")
        pitch_from_min = np.sqrt(8 * H_SCALLOP / max(np.max(kappa_z), 1e-6))
        pitch_from_max = np.sqrt(8 * H_SCALLOP / max(np.min(kappa_z[kappa_z > 0]), 1e-6))
        pitch_from_mean = np.sqrt(8 * H_SCALLOP / max(np.mean(kappa_z), 1e-6))
        print(f"    Max κ ({np.max(kappa_z):.3f} 1/m) → pitch_theory = {pitch_from_min*1000:.2f}mm")
        print(f"    Mean κ ({np.mean(kappa_z):.3f} 1/m) → pitch_theory = {pitch_from_mean*1000:.2f}mm")
        print(f"    Min κ ({np.min(kappa_z[kappa_z > 0]) if np.any(kappa_z > 0) else 0:.3f} 1/m) → pitch_theory = {pitch_from_max*1000:.2f}mm")
        print(f"    PITCH_MAX limit: {PITCH_MAX*1000:.2f}mm")
        print(f"    → Clamping expected for κ < {(8*H_SCALLOP/PITCH_MAX**2):.3f} 1/m")

    return kappa_z, kappa_theta


def compute_adaptive_pitch(
    
    kappa: float,
    h_scallop: float = H_SCALLOP,
    pitch_min: float = PITCH_MIN,
    pitch_max: float = PITCH_MAX
) -> float:
    """
    Compute adaptive pitch based on ISO-scallop theory for flat tool

    For a flat tool on curved surface:
        h = κ · s² / 8

    Solving for step-over s (= pitch in our helical case):
        s = √(8h / κ)

    Args:
        kappa: Local curvature (1/m)
        h_scallop: Target scallop height (m)
        pitch_min: Minimum pitch (m)
        pitch_max: Maximum pitch (m)

    Returns:
        pitch: Adaptive pitch value (m)
    """
    # Avoid division by zero
    kappa_safe = max(kappa, 1e-6)

    # ISO-scallop formula for flat tool
    pitch = np.sqrt(8 * h_scallop / kappa_safe)

    # Clamp to bounds
    pitch = np.clip(pitch, pitch_min, pitch_max)

    return pitch


# -------- Surface Following Helper Functions --------
def compute_radial_normals(closest_points: np.ndarray) -> np.ndarray:
    """
    Compute radial outward normals for surface points

    Safe for convex-like shapes centered at origin.
    Uses direction from origin to surface point.

    Args:
        closest_points: Nx3 array of surface points

    Returns:
        normals: Nx3 array of unit normal vectors (radial outward)
    """
    normals = np.zeros_like(closest_points)

    for i in range(len(closest_points)):
        # Direction from origin to surface point
        radial_dir = closest_points[i] - np.array([0, 0, 0])
        radial_norm = np.linalg.norm(radial_dir)

        if radial_norm > 1e-6:
            normals[i] = radial_dir / radial_norm
        else:
            # Fallback to Z-up if point is at origin
            normals[i] = np.array([0, 0, 1])

    return normals


def compute_blade_engagement(
    closest_points: np.ndarray,
    r_max: float,
    margin: float = ENGAGEMENT_MARGIN
) -> np.ndarray:
    """
    Compute blade engagement states based on radial distance

    Args:
        closest_points: Nx3 array of surface points
        r_max: Maximum radius (unpeeled reference)
        margin: Engagement threshold (m)

    Returns:
        blade_states: N-length array of {0, 1} (0=retracted, 1=engaged)
    """
    r_surface = np.sqrt(closest_points[:, 0]**2 + closest_points[:, 1]**2)
    is_unpeeled = r_surface < (r_max - margin)
    return is_unpeeled.astype(int)


def apply_surface_offset(
    closest_points: np.ndarray,
    normals: np.ndarray,
    offset: float = SURFACE_OFFSET
) -> np.ndarray:
    """
    Apply surface offset along normal direction

    Args:
        closest_points: Nx3 array of surface points
        normals: Nx3 array of unit normals
        offset: Offset distance (m)

    Returns:
        tcp_points: Nx3 array of offset points (TCP positions)
    """
    return closest_points + normals * offset


def apply_surface_following(
    mesh: trimesh.Trimesh,
    positions: np.ndarray,
    r_max: float,
    surface_offset: float = SURFACE_OFFSET
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Apply surface-following logic: find nearest surface, compute offset and blade state

    Args:
        mesh: Trimesh object
        positions: Nx3 array of query points
        r_max: Maximum radius for blade engagement
        surface_offset: Offset distance from surface (m)

    Returns:
        tcp_points: Nx3 array of TCP positions (offset from surface)
        blade_states: N-length array of blade engagement {0, 1}
    """
    print(f"  Finding nearest surface points...")
    closest_points, _, _ = mesh.nearest.on_surface(positions)

    print(f"  Computing radial outward directions...")
    normals = compute_radial_normals(closest_points)

    print(f"  Computing blade engagement...")
    blade_states = compute_blade_engagement(closest_points, r_max)

    print(f"  Applying surface offset...")
    tcp_points = apply_surface_offset(closest_points, normals, surface_offset)

    engagement_ratio = np.sum(blade_states) / len(blade_states) * 100
    print(f"  Blade engagement: {engagement_ratio:.1f}%")

    return tcp_points, blade_states


# -------- Helical Trajectory Generation (Fixed Pitch) --------
def generate_helical_curve(
    mesh: trimesh.Trimesh,
    z_array: np.ndarray,
    theta_array: np.ndarray,
    r_matrix: np.ndarray,
    num_rotations: float = NUM_ROTATIONS,
    surface_offset: float = SURFACE_OFFSET,
    num_samples: int = NUM_SAMPLES
) -> Tuple[np.ndarray, float]:
    """
    Generate helical single-curve trajectory following mesh surface (fixed pitch)

    Args:
        mesh: trimesh object
        z_array: Z heights
        theta_array: Angular samples
        r_matrix: Radial profile (for initial guess)
        num_rotations: Number of full rotations during descent
        surface_offset: Distance from surface (m)
        num_samples: Number of trajectory points

    Returns:
        trajectory: Nx4 array [x, y, z, blade_engaged]
        r_max: Maximum radius
    """
    print(f"\nGenerating helical surface-following trajectory:")
    print(f"  Number of rotations: {num_rotations}")
    print(f"  Trajectory samples: {num_samples}")
    print(f"  Surface offset: {surface_offset*1000:.2f}mm")

    # Create 2D interpolator r(z, θ) for initial guess
    r_interp_func = RectBivariateSpline(z_array, theta_array, r_matrix, kx=1, ky=1)

    # Global max radius
    r_max = np.max(r_matrix)
    print(f"  Global max radius: {r_max:.4f}m")

    # Helix parameters
    z_start = z_array[-1]  # Top
    z_end = z_array[0]     # Bottom
    theta_start = 0.0
    theta_end = 2 * np.pi * num_rotations

    # Generate helix parameter
    t = np.linspace(0, 1, num_samples)

    # Z descends linearly
    z_traj = z_start - t * (z_start - z_end)

    # θ increases linearly (helix)
    theta_traj = theta_start + t * (theta_end - theta_start)

    # Wrap theta to [0, 2π) for r_matrix lookup
    theta_wrapped = theta_traj % (2 * np.pi)

    # Generate initial query points
    print(f"  Generating query points...")
    query_points = np.zeros((num_samples, 3))

    for i in range(num_samples):
        z = z_traj[i]
        theta = theta_wrapped[i]

        # Initial guess from cylindrical profile
        r_guess = float(r_interp_func(z, theta))

        # Query point at initial guess
        query_points[i] = [
            r_guess * np.cos(theta_traj[i]),
            r_guess * np.sin(theta_traj[i]),
            z
        ]

    # Apply surface following
    print(f"  Computing surface-following path...")
    tcp_points, blade_states = apply_surface_following(mesh, query_points, r_max, surface_offset)

    # Combine into trajectory
    trajectory = np.column_stack([tcp_points, blade_states])

    return trajectory, r_max


# -------- Adaptive Helical Trajectory Generation --------
def generate_adaptive_helical_curve(
    mesh: trimesh.Trimesh,
    z_array: np.ndarray,
    theta_array: np.ndarray,
    r_matrix: np.ndarray,
    kappa_z: np.ndarray,
    kappa_theta: np.ndarray,
    h_scallop: float = H_SCALLOP,
    pitch_min: float = PITCH_MIN,
    pitch_max: float = PITCH_MAX,
    surface_offset: float = SURFACE_OFFSET
) -> Tuple[np.ndarray, float]:
    """
    Generate helical trajectory with ISO-scallop adaptive pitch

    Pitch adapts based on local curvature:
    - High curvature → small pitch (dense sampling)
    - Low curvature → large pitch (efficient coverage)

    Args:
        mesh: trimesh object
        z_array: Z heights
        theta_array: Angular samples
        r_matrix: Radial profile
        kappa_z: Curvature map in z-direction (N_z × N_theta)
        kappa_theta: Curvature map in θ-direction
        h_scallop: Target scallop height (m)
        pitch_min: Minimum pitch (m)
        pitch_max: Maximum pitch (m)
        surface_offset: Distance from surface (m)

    Returns:
        trajectory: Nx4 array [x, y, z, blade_engaged]
        r_max: Maximum radius
    """
    print(f"\nGenerating adaptive helical trajectory (ISO-scallop):")
    print(f"  Target scallop height: {h_scallop*1000:.2f}mm")
    print(f"  Pitch range: [{pitch_min*1000:.1f}, {pitch_max*1000:.1f}]mm")
    print(f"  Surface offset: {surface_offset*1000:.2f}mm")

    # Create interpolators
    r_interp = RectBivariateSpline(z_array, theta_array, r_matrix, kx=1, ky=1)
    kappa_interp = RectBivariateSpline(z_array, theta_array, kappa_z, kx=1, ky=1)

    # Global max radius
    r_max = np.max(r_matrix)
    print(f"  Global max radius: {r_max:.4f}m")

    # Start from top
    z_start = z_array[-1]
    z_end = z_array[0]
    theta = 0.0

    # Generate trajectory adaptively
    trajectory_list = []
    z = z_start

    print(f"  Generating adaptive path from Z={z_start:.4f}m to Z={z_end:.4f}m...")

    step_count = 0
    pitches = []
    curvatures = []
    pitch_max_count = 0
    pitch_min_count = 0

    if DEBUG_CURVATURE:
        print(f"\n  DEBUG: First {DEBUG_FIRST_N_STEPS} steps (curvature → pitch calculation):")

    while z > z_end and step_count < MAX_STEPS:
        step_count += 1

        # Get current curvature
        theta_wrapped = theta % (2 * np.pi)
        kappa_current = float(kappa_interp(z, theta_wrapped)[0, 0])
        curvatures.append(kappa_current)

        # Compute adaptive pitch
        pitch = compute_adaptive_pitch(kappa_current, h_scallop, pitch_min, pitch_max)
        pitches.append(pitch)

        # Debug output for first N steps
        if DEBUG_CURVATURE and step_count <= DEBUG_FIRST_N_STEPS:
            kappa_safe = max(kappa_current, 1e-6)
            pitch_unclamped = np.sqrt(8 * h_scallop / kappa_safe)
            clamped_status = ""
            if pitch >= pitch_max:
                clamped_status = " [MAX]"
            elif pitch <= pitch_min:
                clamped_status = " [MIN]"
            print(f"    Step {step_count}: κ={kappa_current:.4f} 1/m → pitch_theory={pitch_unclamped*1000:.2f}mm → pitch={pitch*1000:.2f}mm{clamped_status}")

        # Count clamping events
        if pitch >= pitch_max:
            pitch_max_count += 1
        elif pitch <= pitch_min:
            pitch_min_count += 1

        # Get radius at current position
        r_current = float(r_interp(z, theta_wrapped)[0, 0])

        # Convert to Cartesian
        x = r_current * np.cos(theta)
        y = r_current * np.sin(theta)

        # Store point
        trajectory_list.append([x, y, z])

        # Advance: small angular increment, descend by pitch
        theta += ANGULAR_STEP
        z -= pitch * ANGULAR_STEP / (2 * np.pi)  # Pitch per revolution

        # Debug output every N steps
        if DEBUG_CURVATURE and step_count % DEBUG_STEP_INTERVAL == 0:
            print(f"    Step {step_count}: Z={z:.4f}m, θ={theta/(2*np.pi):.2f} rev, κ={kappa_current:.3f} 1/m, pitch={pitch*1000:.2f}mm")

    print(f"  Generated {step_count} points, {theta/(2*np.pi):.2f} rotations")
    print(f"  Pitch clamping: MAX={pitch_max_count} ({pitch_max_count/step_count*100:.1f}%), MIN={pitch_min_count} ({pitch_min_count/step_count*100:.1f}%)")

    # Convert to numpy array
    positions = np.array(trajectory_list)

    # Apply surface following
    print(f"  Computing surface-following path...")
    tcp_points, blade_states = apply_surface_following(mesh, positions, r_max, surface_offset)

    # Combine into trajectory
    trajectory = np.column_stack([tcp_points, blade_states])

    # Statistics on adaptive pitch and curvature
    pitches = np.array(pitches)
    curvatures = np.array(curvatures)

    print(f"\n  Adaptive pitch statistics:")
    print(f"    Mean: {np.mean(pitches)*1000:.2f}mm")
    print(f"    Min:  {np.min(pitches)*1000:.2f}mm")
    print(f"    Max:  {np.max(pitches)*1000:.2f}mm")
    print(f"    Std:  {np.std(pitches)*1000:.2f}mm")

    # Detailed statistics (debug mode)
    if DEBUG_CURVATURE:
        print(f"\n  Curvature statistics (sampled along path):")
        print(f"    Mean: {np.mean(curvatures):.3f} 1/m")
        print(f"    Min:  {np.min(curvatures):.3f} 1/m")
        print(f"    Max:  {np.max(curvatures):.3f} 1/m")
        print(f"    Std:  {np.std(curvatures):.3f} 1/m")

        print(f"\n  Curvature-Pitch relationship:")
        low_curv_mask = curvatures < np.percentile(curvatures, 25)
        high_curv_mask = curvatures > np.percentile(curvatures, 75)
        print(f"    Low curvature (25%): κ={np.mean(curvatures[low_curv_mask]):.3f} 1/m → pitch={np.mean(pitches[low_curv_mask])*1000:.2f}mm")
        print(f"    High curvature (75%): κ={np.mean(curvatures[high_curv_mask]):.3f} 1/m → pitch={np.mean(pitches[high_curv_mask])*1000:.2f}mm")

    return trajectory, r_max


# -------- Spline Smoothing --------
def smooth_curve_with_spline(
    trajectory: np.ndarray,
    smoothing_factor: float = SPLINE_SMOOTHING
) -> np.ndarray:
    """
    Smooth trajectory with B-spline for natural curve
    Preserve blade state as discrete

    Args:
        trajectory: Nx4 array [x, y, z, blade_engaged]
        smoothing_factor: Spline smoothing parameter (0=interpolation)

    Returns:
        smooth_traj: Nx4 array (smoothed)
    """
    print(f"\nSmoothing curve with B-spline...")

    path_points = trajectory[:, :3]
    blade_states = trajectory[:, 3]

    try:
        # Fit 3D spline
        tck, u = splprep([path_points[:, 0],
                          path_points[:, 1],
                          path_points[:, 2]],
                         s=smoothing_factor,
                         k=3)

        # Evaluate at same resolution
        u_eval = np.linspace(0, 1, len(trajectory))
        smooth_points = np.column_stack(splev(u_eval, tck))

        # Combine with blade states
        smooth_traj = np.column_stack([smooth_points, blade_states])

        print(f"  Smoothed {len(trajectory)} points")
        return smooth_traj

    except Exception as e:
        print(f"  Spline failed: {e}, using original")
        return trajectory


# -------- Visualization --------
def visualize_trajectory(
    mesh: trimesh.Trimesh,
    trajectory: np.ndarray,
    r_max: float
) -> None:
    """
    Visualize mesh and helical trajectory with Open3D

    Args:
        mesh: Trimesh object
        trajectory: Nx4 array [x, y, z, blade_engaged]
        r_max: Maximum radius
    """
    # Mesh (semi-transparent to see path better)
    mesh_o3d = o3d.geometry.TriangleMesh()
    mesh_o3d.vertices = o3d.utility.Vector3dVector(mesh.vertices)
    mesh_o3d.triangles = o3d.utility.Vector3iVector(mesh.faces)
    mesh_o3d.compute_vertex_normals()
    mesh_o3d.paint_uniform_color([0.7, 0.7, 0.7])

    # Extract path
    path_points = trajectory[:, :3]
    blade_states = trajectory[:, 3]

    # Downsample for visualization
    step = max(1, len(path_points) // VIZ_MAX_POINTS)
    viz_points = path_points[::step]
    viz_states = blade_states[::step]

    print(f"\nVisualization info:")
    print(f"  Total points: {len(path_points)}")
    print(f"  Showing: {len(viz_points)} points (every {step}th point)")

    # Create line set
    lines = []
    colors = []

    for i in range(len(viz_points) - 1):
        lines.append([i, i + 1])

        if viz_states[i] == 1:
            colors.append([0.0, 0.9, 0.2])  # Green: cutting
        else:
            colors.append([1.0, 0.5, 0.0])  # Orange: retracted

    line_set = o3d.geometry.LineSet()
    line_set.points = o3d.utility.Vector3dVector(viz_points)
    line_set.lines = o3d.utility.Vector2iVector(lines)
    line_set.colors = o3d.utility.Vector3dVector(colors)

    # Coordinate frame
    coord_frame = o3d.geometry.TriangleMesh.create_coordinate_frame(size=0.04)

    print("\n=== Open3D Visualization ===")
    print("Green curve: Blade engaged (unpeeled regions)")
    print("Orange curve: Already peeled regions (max radius)")
    print("\nNote: Path is continuous (no retraction jumps)")

    o3d.visualization.draw_geometries(
        [mesh_o3d, line_set, coord_frame],
        window_name="Helical Potato Peeling Trajectory",
        width=1600, height=1000
    )


# -------- Export --------
def export_trajectory(trajectory: np.ndarray, output_path: str) -> None:
    """
    Export trajectory to CSV file

    Args:
        trajectory: Nx4 array [x, y, z, blade_engaged]
        output_path: Output file path
    """
    np.savetxt(output_path, trajectory, delimiter=',',
               header='x,y,z,blade_engaged', comments='',
               fmt='%.6f,%.6f,%.6f,%d')
    print(f"\nTrajectory saved: {output_path}")
    print(f"  {len(trajectory)} points")


# -------- Main Pipeline --------
def main():
    """Main trajectory generation pipeline"""
    print("=" * 70)
    print("Helical Potato Peeling Trajectory Generator")
    print("Mode: " + ("ADAPTIVE (ISO-scallop)" if USE_ADAPTIVE else "FIXED PITCH"))
    print("=" * 70)

    # Print configuration
    print(f"\nConfiguration:")
    print(f"  Tool diameter: {TOOL_DIAMETER*1000:.1f} mm (reference)")
    print(f"  Target scallop height: {H_SCALLOP*1000:.2f} mm")
    print(f"  Pitch range: [{PITCH_MIN*1000:.1f}, {PITCH_MAX*1000:.1f}] mm")
    print(f"  PITCH_MAX / Tool diameter: {PITCH_MAX/TOOL_DIAMETER:.1f}x")

    # Step 1: Load mesh
    print("\n[1/5] Loading mesh...")
    mesh = load_mesh(MESH_PATH)

    # Step 2: Compute cylindrical profile
    print("\n[2/5] Computing cylindrical surface profile...")
    z_array, theta_array, r_matrix = compute_cylindrical_profile(
        mesh, NUM_Z_SLICES, NUM_THETA_SAMPLES
    )

    # Step 3: Compute curvature map (for adaptive mode)
    if USE_ADAPTIVE:
        print("\n[3/5] Computing curvature map...")
        kappa_z, kappa_theta = compute_curvature_map(r_matrix, z_array, theta_array)
    else:
        kappa_z, kappa_theta = None, None
        print("\n[3/5] Skipping curvature computation (fixed pitch mode)")

    # Step 4: Generate helical curve
    print("\n[4/5] Generating helical trajectory...")
    if USE_ADAPTIVE:
        trajectory, r_max = generate_adaptive_helical_curve(
            mesh, z_array, theta_array, r_matrix,
            kappa_z, kappa_theta,
            h_scallop=H_SCALLOP,
            pitch_min=PITCH_MIN,
            pitch_max=PITCH_MAX,
            surface_offset=SURFACE_OFFSET
        )
    else:
        trajectory, r_max = generate_helical_curve(
            mesh, z_array, theta_array, r_matrix,
            num_rotations=NUM_ROTATIONS,
            surface_offset=SURFACE_OFFSET,
            num_samples=NUM_SAMPLES
        )

    # Step 5: Smooth curve
    print("\n[5/5] Smoothing curve...")
    smooth_traj = smooth_curve_with_spline(trajectory, SPLINE_SMOOTHING)

    # Export
    export_trajectory(smooth_traj, OUTPUT_PATH)

    # Visualize
    print("\nLaunching visualization...")
    visualize_trajectory(mesh, smooth_traj, r_max)

    print("\n" + "=" * 70)
    print("Pipeline complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
