#!/usr/bin/env python3
"""
곡률-피치 관계 분석 및 파라미터 문제 진단
"""
import numpy as np

print("=" * 70)
print("ISO-scallop 파라미터 문제 분석")
print("=" * 70)

# 현재 설정
TOOL_DIAMETER = 0.010       # 10mm
PITCH_MAX_RATIO = 0.9       # 90%
H_SCALLOP = 0.0003          # 0.3mm
PITCH_MIN = 0.001           # 1mm
PITCH_MAX = TOOL_DIAMETER * PITCH_MAX_RATIO  # 9mm

print(f"\n현재 설정:")
print(f"  TOOL_DIAMETER = {TOOL_DIAMETER*1000:.1f}mm")
print(f"  PITCH_MAX = {PITCH_MAX*1000:.1f}mm (tool의 {PITCH_MAX_RATIO*100:.0f}%)")
print(f"  H_SCALLOP = {H_SCALLOP*1000:.1f}mm")

# 실제 감자에서 측정된 곡률 범위 (이전 실행 결과)
kappa_mean = 3.333
kappa_max = 163.075
kappa_min = 0.001  # 거의 0

print(f"\n실제 감자 곡률 범위:")
print(f"  평균: {kappa_mean:.3f} 1/m")
print(f"  최대: {kappa_max:.3f} 1/m")
print(f"  최소: {kappa_min:.3f} 1/m")

# ISO-scallop 공식으로 계산되는 피치
def calc_pitch(kappa, h_scallop):
    return np.sqrt(8 * h_scallop / max(kappa, 1e-6))

pitch_from_mean = calc_pitch(kappa_mean, H_SCALLOP)
pitch_from_max = calc_pitch(kappa_max, H_SCALLOP)
pitch_from_min = calc_pitch(kappa_min, H_SCALLOP)

print(f"\nISO-scallop 공식으로 계산된 이론적 피치:")
print(f"  평균 곡률 → pitch = {pitch_from_mean*1000:.2f}mm")
print(f"  최대 곡률 → pitch = {pitch_from_max*1000:.2f}mm")
print(f"  최소 곡률 → pitch = {pitch_from_min*1000:.2f}mm")

# 실제 클램핑 후 피치
pitch_actual_mean = np.clip(pitch_from_mean, PITCH_MIN, PITCH_MAX)
pitch_actual_max = np.clip(pitch_from_max, PITCH_MIN, PITCH_MAX)
pitch_actual_min = np.clip(pitch_from_min, PITCH_MIN, PITCH_MAX)

print(f"\n실제 적용되는 피치 (클램핑 후):")
print(f"  평균 곡률 → pitch = {pitch_actual_mean*1000:.2f}mm {'[CLAMPED to MAX]' if pitch_actual_mean == PITCH_MAX else ''}")
print(f"  최대 곡률 → pitch = {pitch_actual_max*1000:.2f}mm {'[CLAMPED to MIN]' if pitch_actual_max == PITCH_MIN else ''}")
print(f"  최소 곡률 → pitch = {pitch_actual_min*1000:.2f}mm {'[CLAMPED to MAX]' if pitch_actual_min == PITCH_MAX else ''}")

# 문제 진단
print("\n" + "=" * 70)
print("⚠️  문제 진단")
print("=" * 70)

print(f"\n현재 PITCH_MAX = {PITCH_MAX*1000:.1f}mm가 너무 작습니다!")
print(f"\n이유:")
print(f"  1. 평균 곡률({kappa_mean:.1f} 1/m)에서 필요한 피치: {pitch_from_mean*1000:.1f}mm")
print(f"  2. 하지만 PITCH_MAX = {PITCH_MAX*1000:.1f}mm로 제한됨")
print(f"  3. 결과: 이론적 피치의 {pitch_actual_mean/pitch_from_mean*100:.1f}%만 사용")

# 곡률 분포 분석
curvatures = np.array([0.1, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0, 100.0])
pitches_theory = np.sqrt(8 * H_SCALLOP / curvatures)
pitches_actual = np.clip(pitches_theory, PITCH_MIN, PITCH_MAX)

print(f"\n곡률별 피치 계산:")
print(f"{'곡률 (1/m)':>12} | {'이론 피치 (mm)':>16} | {'실제 피치 (mm)':>16} | {'클램핑':>10}")
print("-" * 70)

clamped_count = 0
for k, p_theory, p_actual in zip(curvatures, pitches_theory, pitches_actual):
    is_clamped = "MAX" if p_actual == PITCH_MAX else ("MIN" if p_actual == PITCH_MIN else "")
    if is_clamped:
        clamped_count += 1
    print(f"{k:12.1f} | {p_theory*1000:16.2f} | {p_actual*1000:16.2f} | {is_clamped:>10}")

print(f"\n클램핑 비율: {clamped_count}/{len(curvatures)} ({clamped_count/len(curvatures)*100:.0f}%)")

# 해결책 제안
print("\n" + "=" * 70)
print("💡 해결책")
print("=" * 70)

print("\n옵션 1: PITCH_MAX 독립적으로 설정 (권장)")
print("  - PITCH_MAX를 공구 크기와 분리")
print("  - 가공 특성에 따라 설정")
print()
print("  예시:")
print(f"  TOOL_DIAMETER = {TOOL_DIAMETER*1000:.1f}mm  # 공구 크기 (변경 없음)")
print(f"  PITCH_MAX = 0.020  # 20mm (독립 설정)")
print(f"  PITCH_MIN = 0.001  # 1mm")

# 시뮬레이션
PITCH_MAX_NEW = 0.020
pitches_new = np.clip(pitches_theory, PITCH_MIN, PITCH_MAX_NEW)
clamped_new = sum(1 for p in pitches_new if p >= PITCH_MAX_NEW or p <= PITCH_MIN)

print(f"\n  결과 (PITCH_MAX=20mm):")
print(f"    클램핑 비율: {clamped_new}/{len(curvatures)} ({clamped_new/len(curvatures)*100:.0f}%)")
print(f"    평균 곡률 피치: {min(pitch_from_mean, PITCH_MAX_NEW)*1000:.1f}mm")

print("\n옵션 2: H_SCALLOP 감소 + PITCH_MAX 증가")
print("  - 더 높은 품질 요구 + 넓은 피치 범위")
print()
H_SCALLOP_NEW = 0.0002  # 0.3mm → 0.2mm
PITCH_MAX_NEW2 = 0.015  # 15mm
pitch_from_mean_new = calc_pitch(kappa_mean, H_SCALLOP_NEW)
print(f"  H_SCALLOP = {H_SCALLOP_NEW*1000:.1f}mm")
print(f"  PITCH_MAX = {PITCH_MAX_NEW2*1000:.1f}mm")
print(f"  평균 곡률 피치: {min(pitch_from_mean_new, PITCH_MAX_NEW2)*1000:.1f}mm")

print("\n옵션 3: 공구 직경을 '참고용'으로만 사용")
print("  - PITCH_MAX_RATIO를 높여서 제약 완화")
print()
PITCH_MAX_RATIO_NEW = 2.0  # 200% (공구 직경의 2배까지 허용)
PITCH_MAX_NEW3 = TOOL_DIAMETER * PITCH_MAX_RATIO_NEW
print(f"  TOOL_DIAMETER = {TOOL_DIAMETER*1000:.1f}mm")
print(f"  PITCH_MAX_RATIO = {PITCH_MAX_RATIO_NEW} (200%)")
print(f"  PITCH_MAX = {PITCH_MAX_NEW3*1000:.1f}mm")

# 최종 권장
print("\n" + "=" * 70)
print("🎯 최종 권장 설정")
print("=" * 70)

TOOL_DIAMETER_REC = 0.010
PITCH_MIN_REC = 0.001
PITCH_MAX_REC = 0.020  # 독립 설정
H_SCALLOP_REC = 0.0003

print(f"\nTOOL_DIAMETER = {TOOL_DIAMETER_REC*1000:.1f}mm  # 공구 직경 (참고용)")
print(f"PITCH_MIN = {PITCH_MIN_REC*1000:.1f}mm  # 최소 피치 (고곡률)")
print(f"PITCH_MAX = {PITCH_MAX_REC*1000:.1f}mm  # 최대 피치 (저곡률, 독립 설정)")
print(f"H_SCALLOP = {H_SCALLOP_REC*1000:.1f}mm  # 스캘럽 높이")

print(f"\n이 설정에서:")
pitches_rec = np.clip(pitches_theory, PITCH_MIN_REC, PITCH_MAX_REC)
clamped_rec = sum(1 for p in pitches_rec if p >= PITCH_MAX_REC or p <= PITCH_MIN_REC)
print(f"  클램핑 비율: {clamped_rec}/{len(curvatures)} ({clamped_rec/len(curvatures)*100:.0f}%)")
print(f"  평균 곡률 피치: {min(pitch_from_mean, PITCH_MAX_REC)*1000:.1f}mm")
print(f"  피치 범위: {PITCH_MIN_REC*1000:.1f}~{PITCH_MAX_REC*1000:.1f}mm (20배 변화)")

print("\n" + "=" * 70)
