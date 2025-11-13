#!/usr/bin/env python3
"""
곡률 급감 현상 상세 진단

Step 1-10: κ ≈ 45-51 1/m
Step 1000: κ = 0.553 1/m (100배 감소!)

원인 규명
"""
import numpy as np

print("=" * 70)
print("곡률 급감 현상 진단")
print("=" * 70)

# 관찰된 데이터
steps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1000, 2000]
kappas = [48.09, 49.92, 51.30, 51.52, 50.41, 48.09, 45.88, 43.79, 44.21, 45.32, 0.553, 0.889]
z_values = [None] * 10 + [0.0365, 0.0046]  # Step 1000, 2000의 Z 값
theta_revs = [None] * 10 + [1.59, 3.18]

print("\n[1] 관찰된 데이터")
print(f"{'Step':>6} | {'κ (1/m)':>10} | {'Z (m)':>10} | {'θ (rev)':>10} | {'변화율':>10}")
print("-" * 70)

for i, (step, kappa) in enumerate(zip(steps, kappas)):
    z_str = f"{z_values[i]:.4f}" if z_values[i] is not None else "?"
    theta_str = f"{theta_revs[i]:.2f}" if theta_revs[i] is not None else "?"

    if i > 0:
        ratio = kappa / kappas[i-1]
        change = f"{ratio:.3f}x"
    else:
        change = "-"

    print(f"{step:>6} | {kappa:>10.3f} | {z_str:>10} | {theta_str:>10} | {change:>10}")

# 극적인 변화
print(f"\n극적인 변화:")
print(f"  Step 10 → Step 1000:")
print(f"    κ: 45.32 → 0.553 (0.012배, 98.8% 감소!)")
print(f"    Z: ? → 0.0365m")

print("\n" + "=" * 70)
print("[2] 가능한 원인 분석")
print("=" * 70)

print("\n원인 1: Z 위치 변화 (가장 유력)")
print("-" * 40)
print("  관찰:")
print("    - Step 1-10: Z 위치 불명 (시작 부근, Top으로 추정)")
print("    - Step 1000: Z=0.0365m")
print("    - Step 2000: Z=0.0046m")
print()
print("  추론:")
print("    - Path는 Top(z_max)에서 시작해서 Bottom(z_min)으로 하강")
print("    - Step 1000에서 Z=0.0365m는 이미 중간~하단부")
print("    - 감자 형태: Top(뾰족/불규칙) → Bottom(둥글고 완만)")
print()
print("  감자 형태 가설:")
print("    ┌─────┐")
print("    │ /^\ │  ← Top: 뾰족, 높은 곡률 (κ ≈ 50 1/m)")
print("    │ | | │")
print("    │ | | │  ← Middle: 원통형, 중간 곡률")
print("    │ \_/ │  ← Bottom: 둥글고 완만, 낮은 곡률 (κ ≈ 0.5-1 1/m)")
print("    └─────┘")

print("\n원인 2: 곡률 맵의 문제")
print("-" * 40)
print("  가능성:")
print("    - Gaussian smoothing (σ=1.5)이 경계 근처에서 과도하게 평활화")
print("    - 경계 효과: mesh 끝부분에서 곡률 계산 불안정")
print()
print("  확인 방법:")
print("    - 곡률 맵 통계에서 min, mean, max 값 확인")
print("    - Z=0.0365m가 mesh bounds 중 어디쯤인지 확인")

print("\n원인 3: θ 방향 영향 (가능성 낮음)")
print("-" * 40)
print("  관찰:")
print("    - Step 1-10: θ ≈ 0-0.1 rad (시작 부근)")
print("    - Step 1000: θ = 1.59 rev = 10 rad")
print("    - Step 2000: θ = 3.18 rev = 20 rad")
print()
print("  분석:")
print("    - θ가 변하더라도 한 바퀴 돌면 다시 비슷한 곡률이어야 함")
print("    - 하지만 κ가 계속 낮은 것으로 보아 θ보다 Z의 영향이 지배적")

print("\n" + "=" * 70)
print("[3] 검증 방법")
print("=" * 70)

print("\n1. 프로그램 출력에서 다음을 확인:")
print("   ─────────────────────────────────────")
print("   Mesh loaded: ...")
print("   Z bounds: [Z_min, Z_max]m")
print()
print("   Curvature statistics:")
print("     κ_z (radial):  mean=?, max=?, min=?")
print("   ─────────────────────────────────────")

print("\n2. 계산:")
z_1000 = 0.0365
z_2000 = 0.0046

print(f"\n   만약 Z bounds = [0.000, 0.080]m 라면:")
print(f"     Step 1 위치: Z ≈ 0.080m (Top, 100%)")
print(f"     Step 1000: Z = {z_1000:.4f}m ({z_1000/0.08*100:.1f}%)")
print(f"     Step 2000: Z = {z_2000:.4f}m ({z_2000/0.08*100:.1f}%)")
print()
print(f"   → Step 1000은 이미 중간~하단부 (약 45% 높이)")
print(f"   → Step 2000은 거의 바닥 (약 6% 높이)")

print("\n3. 곡률-높이 관계 추정:")
print()
print("   만약 감자가 타원체 모양이라면:")
print("     - Top (z=100%): 곡률 높음 (뾰족)")
print("     - Middle (z=50%): 곡률 중간")
print("     - Bottom (z=0%): 곡률 낮음 (완만)")
print()
print("   이것은 정상적인 현상일 수 있습니다!")

print("\n" + "=" * 70)
print("[4] 해결책")
print("=" * 70)

print("\n만약 이것이 정상적인 형태라면:")
print("  → 해결 불필요, ISO-scallop이 의도대로 작동 중")
print("  → 높은 곡률: 조밀한 경로")
print("  → 낮은 곡률: 성긴 경로 (효율적)")

print("\n만약 문제가 있다면:")
print("  1. PITCH_MAX 조정")
print("     - 현재: 20mm (77% 클램핑)")
print("     - 권장: 10-15mm (클램핑 < 50%)")
print()
print("  2. H_SCALLOP 조정")
print("     - 현재: 3mm")
print("     - 더 조밀: 2mm")
print("     - 더 성긴: 5mm")
print()
print("  3. 곡률 평활화 감소")
print("     - CURVATURE_SMOOTHING_SIGMA: 1.5 → 1.0")

print("\n" + "=" * 70)
print("[5] 권장 조치")
print("=" * 70)

print("""
1. 프로그램 전체 출력 확인:
   - Z bounds 값
   - Curvature statistics (mean, max, min)
   - Step 1의 Z 값 (첫 step 디버그 출력에 추가 필요)

2. 계산:
   Total Z range = Z_max - Z_min
   Step 1 위치 = Top (Z_max)
   Step 1000 위치 = 어디? (Z_max 대비 %)

3. 시각화 확인:
   - Open3D visualization에서 mesh 형태 확인
   - 위쪽이 뾰족한지, 아래쪽이 완만한지 육안 확인

4. 만약 정상이라면:
   - PITCH_MAX를 낮춰서 클램핑 비율 감소
   - 또는 그대로 사용 (77% 클램핑은 높지만 작동은 함)

5. 만약 비정상이라면:
   - Gaussian smoothing 조정
   - Curvature window 크기 조정
   - 경계 처리 개선
""")

print("\n" + "=" * 70)
print("다음 명령으로 Z bounds와 Curvature statistics를 확인하세요:")
print("  uv run python src/potato_peeling_path_generator.py | grep -A 5 'Z bounds\\|Curvature statistics'")
print("=" * 70)
