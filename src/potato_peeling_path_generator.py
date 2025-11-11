#!/usr/bin/env python3
# ------------------------------
# Module: potato_peeling_path_generator.py
# Purpose:
#   - Generate helical single-curve trajectory for rotating potato
#   - Robot moves in continuous spiral while potato rotates
#   - Blade follows surface contour, engages only on unpeeled regions
#   - Output: smooth 3D curve synchronized with rotation
# Notes:
#   - Potato rotates on Z-axis
#   - Robot traces helical path (θ increases with Z descent)
#   - Radial distance modulated to follow surface
# ------------------------------

import numpy as np
import trimesh
import open3d as o3d
from scipy.interpolate import RectBivariateSpline, splprep, splev
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt


# -------- Mesh Loading --------
def load_mesh(mesh_path):
    """Load STL mesh and compute bounds"""
    mesh = trimesh.load(mesh_path)
    print(f"Mesh loaded: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print(f"Z bounds: [{mesh.bounds[0][2]:.4f}, {mesh.bounds[1][2]:.4f}]m")
    return mesh


# -------- Radial Profile Extraction --------
def compute_cylindrical_profile(mesh, num_z_slices=120, num_theta_samples=360):
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


# -------- Helical Trajectory Generation with Surface Following --------
def generate_helical_curve(mesh, z_array, theta_array, r_matrix,
                          num_rotations=2.0,
                          surface_offset=0.0001,  # 0.1mm from surface
                          num_samples=5000):
    """
    Generate helical single-curve trajectory following mesh surface

    Uses actual mesh geometry to compute precise surface-following path:
    - Ray casting to find closest surface point
    - Normal-based offset for accurate surface tracking
    - 0.1mm offset from mesh surface

    Args:
        mesh: trimesh object
        z_array: Z heights
        theta_array: Angular samples
        r_matrix: Radial profile (for initial guess)
        num_rotations: Number of full rotations during descent
        surface_offset: Distance from surface (m, default 0.1mm)
        num_samples: Number of trajectory points

    Returns:
        trajectory: Nx4 array [x, y, z, blade_engaged]
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

    # Generate trajectory with surface following
    trajectory = np.zeros((num_samples, 4))

    print(f"  Computing surface-following path...")

    # Batch process for speed
    query_points = []
    for i in range(num_samples):
        z = z_traj[i]
        theta = theta_wrapped[i]

        # Initial guess from cylindrical profile
        r_guess = float(r_interp_func(z, theta))

        # Query point at initial guess
        x_guess = r_guess * np.cos(theta_traj[i])
        y_guess = r_guess * np.sin(theta_traj[i])
        query_points.append([x_guess, y_guess, z])

    query_points = np.array(query_points)

    # Batch query nearest surface points
    print(f"  Finding nearest surface points...")
    closest_points, distances, face_ids = mesh.nearest.on_surface(query_points)

    # Get surface normals
    normals = mesh.face_normals[face_ids]

    # Compute blade engagement
    engagement_margin = 0.001  # 1mm
    retract_distance = 0.01    # 10mm

    print(f"  Applying surface offset and blade logic...")
    for i in range(num_samples):
        if i % 1000 == 0:
            print(f"    Sample {i}/{num_samples}")

        # Check if this is an unpeeled region
        r_surface = np.sqrt(closest_points[i, 0]**2 + closest_points[i, 1]**2)
        is_unpeeled = r_surface < (r_max - engagement_margin)

        if is_unpeeled:
            # Blade engaged: follow surface with precise offset
            # Move along surface normal
            tcp_point = closest_points[i] + normals[i] * surface_offset
            blade_engaged = 1
        else:
            # Blade retracted: move away from max radius zone
            theta = theta_traj[i]
            r_retract = r_max + retract_distance
            tcp_point = np.array([
                r_retract * np.cos(theta),
                r_retract * np.sin(theta),
                z_traj[i]
            ])
            blade_engaged = 0

        trajectory[i] = [tcp_point[0], tcp_point[1], tcp_point[2], blade_engaged]

    engagement_ratio = np.sum(trajectory[:, 3]) / num_samples * 100
    print(f"  Blade engagement: {engagement_ratio:.1f}%")

    # Debug info
    engaged_points = trajectory[trajectory[:, 3] == 1, :3]
    if len(engaged_points) > 0:
        # Check actual distance from surface for engaged points
        actual_distances = np.linalg.norm(engaged_points - closest_points[trajectory[:, 3] == 1], axis=1)
        print(f"\n  Debug - Surface offset accuracy (engaged points):")
        print(f"    Target offset: {surface_offset*1000:.2f}mm")
        print(f"    Mean actual offset: {np.mean(actual_distances)*1000:.2f}mm")
        print(f"    Min offset: {np.min(actual_distances)*1000:.2f}mm")
        print(f"    Max offset: {np.max(actual_distances)*1000:.2f}mm")

    return trajectory, r_max


# -------- Spline Smoothing --------
def smooth_curve_with_spline(trajectory, smoothing_factor=0.0):
    """
    Smooth trajectory with B-spline for natural curve
    Preserve blade state as discrete
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
def visualize_trajectory(mesh, trajectory, r_max):
    """
    Visualize mesh and helical trajectory with Open3D
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

    # Downsample for visualization (show more points for detail)
    step = max(1, len(path_points) // 5000)
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

    # Start/End markers
    start_marker = o3d.geometry.TriangleMesh.create_sphere(radius=0.008)
    start_marker.translate(path_points[0])
    start_marker.paint_uniform_color([0.2, 0.4, 1.0])  # Blue

    end_marker = o3d.geometry.TriangleMesh.create_sphere(radius=0.008)
    end_marker.translate(path_points[-1])
    end_marker.paint_uniform_color([1.0, 0.2, 0.2])  # Red

    # Max radius cylinder (for reference)
    cylinder = o3d.geometry.TriangleMesh.create_cylinder(
        radius=r_max, height=mesh.bounds[1][2] - mesh.bounds[0][2]
    )
    cylinder.translate([0, 0, mesh.bounds[0][2]])
    cylinder.paint_uniform_color([0.3, 0.3, 0.8])
    cylinder_wireframe = o3d.geometry.LineSet.create_from_triangle_mesh(cylinder)
    cylinder_wireframe.paint_uniform_color([0.5, 0.5, 0.9])

    print("\n=== Open3D Visualization ===")
    print("Green curve: Blade engaged (cutting unpeeled regions)")
    print("Orange curve: Blade retracted (avoiding max radius)")
    print("Blue wireframe cylinder: Global max radius reference")
    print("Blue sphere: Start (top)")
    print("Red sphere: End (bottom)")

    o3d.visualization.draw_geometries(
        [mesh_o3d, line_set, coord_frame, start_marker, end_marker, cylinder_wireframe],
        window_name="Helical Potato Peeling Trajectory",
        width=1600, height=1000
    )


# -------- Export --------
def export_trajectory(trajectory, output_path):
    """Export trajectory: x, y, z, blade_engaged"""
    np.savetxt(output_path, trajectory, delimiter=',',
               header='x,y,z,blade_engaged', comments='',
               fmt='%.6f,%.6f,%.6f,%d')
    print(f"\nTrajectory saved: {output_path}")
    print(f"  {len(trajectory)} points")


# -------- Main Pipeline --------
def main():
    # -------- Configuration --------
    MESH_PATH = "mesh/mesh_1_processed.stl"
    OUTPUT_PATH = "output/peeling_trajectory.csv"

    # Mesh sampling resolution
    NUM_Z_SLICES = 120
    NUM_THETA_SAMPLES = 360

    # Helix parameters
    NUM_ROTATIONS = 15.0     # Number of full rotations during descent (more = tighter spacing)
    NUM_SAMPLES = 10000      # Trajectory resolution
    SURFACE_OFFSET = 0.0001  # 0.1mm offset from surface

    # Smoothing
    SPLINE_SMOOTHING = 0.0   # 0 = interpolation, >0 = smoothing

    print("=" * 70)
    print("Helical Potato Peeling Trajectory Generator")
    print("=" * 70)

    # Step 1: Load mesh
    print("\n[1/4] Loading mesh...")
    mesh = load_mesh(MESH_PATH)

    # Step 2: Compute cylindrical profile
    print("\n[2/4] Computing cylindrical surface profile...")
    z_array, theta_array, r_matrix = compute_cylindrical_profile(
        mesh, NUM_Z_SLICES, NUM_THETA_SAMPLES
    )

    # Step 3: Generate helical curve
    print("\n[3/4] Generating helical trajectory...")
    trajectory, r_max = generate_helical_curve(
        mesh, z_array, theta_array, r_matrix,
        num_rotations=NUM_ROTATIONS,
        surface_offset=SURFACE_OFFSET,
        num_samples=NUM_SAMPLES
    )

    # Step 4: Smooth curve
    print("\n[4/4] Smoothing curve...")
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
