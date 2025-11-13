#!/usr/bin/env python3
"""
곡률 감소 현상 분석 도구

Step이 진행될수록 곡률이 감소하는 원인을 파악
"""
import numpy as np
import trimesh
from scipy.interpolate import RectBivariateSpline

# 파라미터 (main 파일과 동일)
MESH_PATH = "mesh/mesh_1_processed.stl"
NUM_Z_SLICES = 240
NUM_THETA_SAMPLES = 360

print("=" * 70)
print("곡률 감소 현상 분석")
print("=" * 70)

# 1. Mesh 로드
print("\n[1] Loading mesh...")
mesh = trimesh.load(MESH_PATH)
z_min, z_max = mesh.bounds[0][2], mesh.bounds[1][2]
print(f"  Z bounds: [{z_min:.4f}, {z_max:.4f}]m")
print(f"  Z range: {z_max - z_min:.4f}m")

# 2. 간단한 cylindrical profile (샘플링)
print("\n[2] Sampling cylindrical profile...")
z_array = np.linspace(z_min, z_max, NUM_Z_SLICES)
theta_array = np.linspace(0, 2*np.pi, NUM_THETA_SAMPLES, endpoint=False)

# 몇 개 Z 레벨에서만 샘플링
sample_z_indices = [0, NUM_Z_SLICES//4, NUM_Z_SLICES//2, 3*NUM_Z_SLICES//4, NUM_Z_SLICES-1]
sample_z_values = [z_array[i] for i in sample_z_indices]
sample_labels = ["Bottom (z_min)", "25%", "Middle", "75%", "Top (z_max)"]

print(f"\n  Sampling Z levels:")
for label, z_val in zip(sample_labels, sample_z_values):
    print(f"    {label:20s}: Z={z_val:.4f}m")

# 3. 각 Z 레벨에서 반지름 분석
print(f"\n[3] Analyzing radius at each Z level...")

for label, z_val in zip(sample_labels, sample_z_values):
    # Slice mesh at this Z
    slice_result = mesh.section(plane_origin=[0, 0, z_val],
                               plane_normal=[0, 0, 1])

    if slice_result is None:
        print(f"    {label:20s}: No intersection")
        continue

    try:
        planar, _ = slice_result.to_2D()
        vertices_2d = planar.vertices

        if len(vertices_2d) < 3:
            print(f"    {label:20s}: Insufficient vertices")
            continue

        # Compute radius
        r_values = np.sqrt(vertices_2d[:, 0]**2 + vertices_2d[:, 1]**2)
        r_mean = np.mean(r_values)
        r_std = np.std(r_values)
        r_min = np.min(r_values)
        r_max = np.max(r_values)

        print(f"    {label:20s}: r_mean={r_mean*1000:.1f}mm, r_std={r_std*1000:.1f}mm, r_range=[{r_min*1000:.1f}, {r_max*1000:.1f}]mm")

    except Exception as e:
        print(f"    {label:20s}: Error - {e}")

# 4. 곡률 추정 (간단한 방법)
print(f"\n[4] Estimating curvature trend (simple method)...")
print(f"  NOTE: 실제 곡률 계산은 quadric fitting을 사용하지만, 여기서는 간단한 방법 사용")

# Z 방향 곡률 근사: 반지름의 2차 미분
sample_r_means = []
for z_val in sample_z_values:
    slice_result = mesh.section(plane_origin=[0, 0, z_val],
                               plane_normal=[0, 0, 1])
    if slice_result is not None:
        try:
            planar, _ = slice_result.to_2D()
            vertices_2d = planar.vertices
            r_values = np.sqrt(vertices_2d[:, 0]**2 + vertices_2d[:, 1]**2)
            sample_r_means.append(np.mean(r_values))
        except:
            sample_r_means.append(0)
    else:
        sample_r_means.append(0)

sample_r_means = np.array(sample_r_means)

# 2차 미분 근사 (중심 차분)
print(f"\n  Radius profile along Z:")
for i, (label, z_val, r_mean) in enumerate(zip(sample_labels, sample_z_values, sample_r_means)):
    print(f"    {label:20s}: Z={z_val:.4f}m, r={r_mean*1000:.1f}mm")

# 5. Path simulation
print(f"\n[5] Simulating path generation...")
print(f"  Starting from TOP (z_max) and descending to BOTTOM (z_min)")

z_start = z_array[-1]  # Top
z_end = z_array[0]     # Bottom
theta = 0.0
ANGULAR_STEP = 0.01
PITCH_TEST = 0.010  # Fixed pitch for test

print(f"\n  Simulating first 10 steps:")
print(f"  {'Step':>5} | {'Z (m)':>10} | {'Z (mm)':>10} | {'θ (rad)':>10} | {'θ (deg)':>10} | Direction")
print(f"  {'-'*80}")

z = z_start
for step in range(1, 11):
    theta_deg = theta * 180 / np.pi
    z_mm = z * 1000

    direction = "TOP" if step == 1 else ("BOTTOM" if z <= z_end else "→")

    print(f"  {step:>5} | {z:>10.4f} | {z_mm:>10.1f} | {theta:>10.4f} | {theta_deg:>10.1f} | {direction}")

    # Advance
    theta += ANGULAR_STEP
    z -= PITCH_TEST * ANGULAR_STEP / (2 * np.pi)

print(f"\n  → Z는 감소 (하강)")
print(f"  → θ는 증가 (회전)")

# 6. 가능한 원인 분석
print(f"\n[6] 가능한 원인 분석")
print(f"=" * 70)

print(f"\n원인 1: 감자 형태 - 위쪽이 뾰족, 아래쪽이 완만")
print(f"  증상: Top에서 시작 → 높은 곡률")
print(f"        Bottom으로 하강 → 낮은 곡률")
print(f"  확인: 위의 radius profile을 보세요")
print(f"        - Top의 r_std가 크면 (불규칙) → 높은 곡률")
print(f"        - Bottom의 r_std가 작으면 (원형) → 낮은 곡률")

print(f"\n원인 2: Z 방향 혼동")
print(f"  증상: Mesh의 위아래가 뒤집혀 있을 수 있음")
print(f"  확인: Z bounds가 예상과 맞는지 확인")
print(f"        - 일반적으로 감자는 Z축을 따라 세로로 놓임")

print(f"\n원인 3: θ 방향 영향")
print(f"  증상: 특정 θ 범위에서 곡률이 낮음")
print(f"  확인: 한 바퀴(2π) 회전하는 동안 곡률이 다시 올라가는지 확인")
print(f"        - 만약 계속 감소한다면 θ가 아닌 Z의 영향")

print(f"\n원인 4: 곡률 계산 버그")
print(f"  증상: Interpolation이나 인덱싱 문제")
print(f"  확인: compute_curvature_map() 출력에서")
print(f"        - κ_z statistics를 확인")
print(f"        - Expected pitch range를 확인")

print(f"\n[7] 권장 디버깅 단계")
print(f"=" * 70)

print(f"""
1. 실제 프로그램 실행:
   $ uv run python src/potato_peeling_path_generator.py

2. 디버그 출력 확인 (DEBUG_CURVATURE=True):
   - "DEBUG: First 10 steps" 섹션에서 κ 값 추이
   - Step 0, 1000, 2000, ... 에서 Z, θ, κ 값

3. Curvature map statistics 확인:
   - κ_z: mean, max, min
   - Expected pitch range
   - "Clamping expected for κ < X" 메시지

4. 만약 κ가 계속 감소한다면:
   a. θ를 2π 증가시켰을 때 κ가 다시 원래 값으로 돌아오는지 확인
   b. 돌아오면 → θ 방향 영향 (정상)
   c. 돌아오지 않으면 → Z 방향 영향 (감자 형태 or 버그)

5. Mesh 시각화:
   - Open3D visualization에서 mesh 확인
   - 위아래가 예상과 맞는지 확인
""")

print(f"\n실행 명령:")
print(f"  uv run python analyze_curvature_trend.py")
print(f"\n그 다음 실제 path generator를 실행하고 출력을 비교하세요:")
print(f"  uv run python src/potato_peeling_path_generator.py")

print(f"\n" + "=" * 70)
