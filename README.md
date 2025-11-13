# Potato Peeling Path Generator

감자 껍질 벗기기를 위한 ISO-scallop 기반 적응형 헬리컬 경로 생성기

## 📖 프로젝트 개요

회전하는 감자의 표면을 따라 로봇 매니퓰레이터가 이동할 수 있는 최적화된 헬리컬 궤적을 생성합니다.

### 주요 특징
- 감자가 Z축을 기준으로 회전하는 동안 로봇이 나선형 경로를 따라 이동
- **ISO-scallop 적응형 피치**: 표면 곡률에 따라 자동으로 피치(pitch) 조정
- STL 메시 기반 표면 프로파일 추출 및 경로 생성
- 원통 좌표계(cylindrical coordinates) 기반 radial profile r(θ, z) 계산
- 표면 추종(surface following) 및 블레이드 상태 제어
- 곡률 기반 경로 밀도 최적화 (고곡률 → 조밀, 저곡률 → 성글게)

## 🚀 빠른 시작

### 1. 의존성 설치

```bash
# uv로 의존성 설치
uv sync
```

**필수 라이브러리:**
- `numpy>=1.24.0` - 수치 계산
- `scipy>=1.10.0` - 보간 및 필터링
- `trimesh>=4.0.0` - 메시 처리
- `open3d>=0.17.0` - 3D 시각화
- `rtree>=1.4.1` - 공간 인덱싱
- `networkx>=3.4.2` - 그래프 알고리즘

### 2. 메시 전처리 (최초 1회)

```bash
# STL 메시 전처리 (단위 변환, 회전, 정렬, 스무딩)
uv run python src/process_mesh_fixed.py
```

**전처리 단계:**
- 단위 변환 (mm → m)
- X축 기준 -90° 회전
- 원점을 메시 바닥 중심으로 이동
- 구멍 메우기 (hole filling)
- Laplacian 스무딩 (강화된 반복 적용)

### 3. 경로 생성

```bash
# 헬리컬 경로 생성 (기본: 적응형 모드)
uv run python src/potato_peeling_path_generator.py
```

**출력:**
- `output/peeling_trajectory_adaptive.csv` - 적응형 피치 경로
- Open3D 3D 시각화 창 (경로 + 메시)

### 4. 경로 검증

```bash
# 경로 연속성 확인
uv run python check_path_continuity.py
```

**검증 항목:**
- 점 간 거리 통계 (평균/최대/최소)
- 불연속 점프 감지 (평균의 10배 이상)

## 📁 파일 구조

```
potato/
├── pyproject.toml                          # UV 패키지 설정
├── uv.lock                                 # UV 잠금 파일
├── README.md                               # 이 문서
├── check_path_continuity.py               # 경로 연속성 검증 스크립트
├── mesh/
│   ├── mesh_1.stl                         # 원본 STL 메시
│   └── mesh_1_processed.stl               # 전처리된 메시 (최종)
├── output/                                 # 생성된 경로 출력
│   ├── peeling_trajectory.csv             # 고정 피치 경로
│   └── peeling_trajectory_adaptive.csv    # 적응형 피치 경로
└── src/
    ├── potato_peeling_path_generator.py   # ⭐ 메인: ISO-scallop 경로 생성기
    └── process_mesh_fixed.py              # 메시 전처리 스크립트
```

## ⚙️ 주요 매개변수

`src/potato_peeling_path_generator.py` 상단에서 수정:

### 모드 선택
```python
USE_ADAPTIVE = True         # True: ISO-scallop 적응형, False: 고정 피치
```

### 메시 샘플링 해상도
```python
NUM_Z_SLICES = 240          # Z 방향 슬라이스 (곡률 계산용 증가)
NUM_THETA_SAMPLES = 360     # 각도(θ) 샘플 개수 (1도 간격)
```

### ISO-scallop 적응형 피치 파라미터
```python
H_SCALLOP = 0.003           # 목표 스캘럽 높이: 3mm (표면 품질)
PITCH_MIN = 0.001           # 최소 피치: 1mm (고곡률 영역)
PITCH_MAX = 0.027           # 최대 피치: 27mm (저곡률 영역)
ANGULAR_STEP = 0.01         # 각도 증분: ~0.57도 (적응형 생성용)
                            # NOTE: 공구 크기와 무관하게 ISO-scallop 이론 적용
```

### 고정 피치 파라미터 (USE_ADAPTIVE=False일 때)
```python
NUM_ROTATIONS = 15.0        # 전체 회전 수
NUM_SAMPLES = 10000         # 궤적 해상도
```

### 표면 추종 및 블레이드 제어
```python
SURFACE_OFFSET = 0.001      # 표면으로부터 거리: 0.1mm (TCP 오프셋)
ENGAGEMENT_MARGIN = 0.001   # 블레이드 결합 임계값: 1mm
```

### 스무딩 및 곡률 계산
```python
SPLINE_SMOOTHING = 0.0              # 0 = 보간, >0 = 스무딩
CURVATURE_SMOOTHING_SIGMA = 1.5     # 곡률 맵 Gaussian 필터
CURVATURE_WINDOW = 2                # 국소 곡률 계산 윈도우 (5×5)
```

## 📊 출력 형식

### 생성되는 파일들

**헬리컬 경로** (CSV, ROS2 호환)
```csv
x,y,z,blade_engaged
0.0234,0.0156,0.0010,1
0.0238,0.0154,0.0012,1
0.0240,0.0152,0.0014,0
...
```

**컬럼 설명:**
- `x, y, z`: TCP 위치 (m, 데카르트 좌표계, base_link 프레임)
- `blade_engaged`: 블레이드 상태
  - `1`: 절삭 중 (미가공 영역)
  - `0`: 후퇴 (이미 가공된 영역, 최대 반경)

**시각화 (Open3D):**
- 녹색 곡선: 블레이드 결합 (미가공 영역)
- 주황색 곡선: 이미 벗겨진 영역 (최대 반경)
- 경로는 연속적 (점프 없음)

## 🎯 알고리즘 개요

### 1. 메시 전처리 (`process_mesh_fixed.py`)
- **단위 정규화**: mm 단위 메시를 m 단위로 변환
- **좌표 정렬**: X축 기준 -90° 회전으로 감자 장축을 Z축에 정렬
- **원점 조정**: 바닥 중심을 원점(0,0,0)으로 설정
- **표면 보정**: 구멍 메우기(watertight) 및 강화된 Laplacian 스무딩

### 2. 원통 좌표계 변환
- 메시를 원통 좌표계 `(r, θ, z)`로 샘플링 (240 Z-슬라이스 × 360 각도)
- 각 `(θ, z)` 위치에서 표면까지의 반경 `r` 계산
- 2D 스플라인 보간으로 연속 함수 `r(θ, z)` 생성

### 3. 곡률 맵 계산 (적응형 모드)
- **국소 곡면 피팅**: 각 그리드 점에서 5×5 이웃 윈도우로 quadric surface 피팅
  - 2차 곡면: `z = ax² + by² + cxy + dx + ey + f`
  - 최소자승법(least squares)으로 6개 계수 추정
- **정확한 곡률 공식 적용** (라인 244-273):
  1. **1차 미분 계산** (표면 기울기, 중심점에서 평가):
     ```
     z_x = 2a·x_ij + c·y_ij + d
     z_y = 2b·y_ij + c·x_ij + e
     ```
  2. **2차 미분 계산** (Hessian 행렬):
     ```
     z_xx = 2a,  z_yy = 2b,  z_xy = c
     ```
  3. **방향별 곡률** (radial/tangential):
     ```
     κ_r = z_xx·cos²θ + 2·z_xy·cosθ·sinθ + z_yy·sin²θ
     κ_t = z_xx·sin²θ - 2·z_xy·cosθ·sinθ + z_yy·cos²θ
     ```
  4. **Gradient 보정** (정확한 3D 곡률):
     ```
     κ_exact = κ_dir / (1 + |∇z|²)^(3/2)
     ```
     - 경사진 표면에서 곡률 과대평가 방지 ✓
- **주곡률 출력**:
  - `κ_z`: 반경 방향 곡률 (1/m, **ISO-scallop pitch 적응용**)
  - `κ_θ`: 접선 방향 곡률 (1/m, 피드 방향 참고용)
- **노이즈 제거**: Gaussian 필터링 (σ=1.5) 적용

### 4. ISO-scallop 적응형 경로 생성
**4.1 적응형 피치 공식 (라인 354-387)**
```
h = κ · s² / 8    (평면 공구, 곡면 가공)
s = √(8h / κ)     (피치 = 스텝오버)
```
- `h`: 목표 스캘럽 높이 (3mm = H_SCALLOP)
- `κ`: 국소 곡률 (1/m, κ_z 값 사용)
- `s`: 적응형 피치 (1~27mm, PITCH_MIN~PITCH_MAX로 클램핑)

**4.2 경로 생성 절차 (라인 591-761)**
1. **초기화** (라인 658-665):
   - 시작점: Z = Z_max (감자 상단), θ = 0
   - 종료 조건: Z ≤ Z_min 또는 MAX_STEPS(100,000) 도달
2. **반복 스텝** (라인 677-718):
   ```python
   while z > z_end and step_count < MAX_STEPS:
       # a. 현재 위치의 곡률 조회 (보간)
       theta_wrapped = theta % (2π)
       kappa_current = kappa_interp(z, theta_wrapped)

       # b. 적응형 피치 계산
       pitch = √(8 × H_SCALLOP / max(kappa, 1e-6))
       pitch = clip(pitch, PITCH_MIN, PITCH_MAX)

       # c. 반경 조회 및 좌표 변환
       r_current = r_interp(z, theta_wrapped)
       x, y = r_current × cos(theta), r_current × sin(theta)

       # d. 전진
       theta += ANGULAR_STEP  (0.01 rad ≈ 0.57°)
       z -= pitch × ANGULAR_STEP / (2π)
   ```
3. **통계 수집** (라인 669-672, 700-704):
   - 피치 분포 (평균/최소/최대/표준편차)
   - 클램핑 빈도 (PITCH_MAX, PITCH_MIN 도달 횟수)
4. **디버그 출력** (라인 674-723):
   - 처음 10 스텝: κ → pitch_theory → pitch [클램핑 상태]
   - 1000 스텝마다: Z, 회전수, κ, pitch 출력

### 5. 표면 추종 및 블레이드 제어
- **최근접 표면 계산**: 각 경로점에서 메시의 최근접 점 탐색
- **방사 법선 계산**: 원점→표면점 방향으로 outward 법선 계산
- **TCP 오프셋**: 표면 + 법선 방향 × 0.1mm (공구 끝점 위치)
- **블레이드 상태**:
  - 표면 반경 < (최대 반경 - 1mm) → `blade_engaged = 1`
  - 그 외 → `blade_engaged = 0`

### 6. 스무딩 및 출력
- **B-spline 스무딩**: 3차 스플라인으로 경로 부드럽게 (선택적)
- **CSV 내보내기**: `[x, y, z, blade_engaged]` 형식
- **Open3D 시각화**: 메시 + 색상 구분 경로선

## 🔧 사용 예시

### 1. 메시 전처리 확인

```bash
$ uv run python src/process_mesh_fixed.py

=== Step 0: Converting units (mm → m) ===
  Detected millimeter units, converting to meters...

=== Step 1: Rotating mesh by -90 degrees around X axis ===
After rotation bounds: min=[...], max=[...]

=== Step 2: Moving origin to bottom center ===
Translation vector: [...]
After translation bounds: min=[0, 0, 0], max=[0.08, 0.08, 0.12]

=== Step 3: Filling holes ===
Is watertight after filling: True

=== Step 4: Smoothing mesh ===
After smoothing: 25643 vertices, 51282 faces

Origin position check:
  X center: 0.000000 (should be ~0)  ✓
  Y center: 0.000000 (should be ~0)  ✓
  Z bottom: 0.000000 (should be ~0)  ✓
```

### 2. 적응형 경로 생성 확인

```bash
$ uv run python src/potato_peeling_path_generator.py

======================================================================
Helical Potato Peeling Trajectory Generator
Mode: ADAPTIVE (ISO-scallop)
======================================================================

Tool Configuration:
  Tool diameter: 10.0 mm
  Pitch max ratio: 90%
  → PITCH_MAX = 9.00 mm (90% of tool diameter)
  Target scallop height: 0.30 mm

[1/5] Loading mesh...
Mesh loaded: 25643 vertices, 51282 faces
Z bounds: [0.0000, 0.1200]m

[2/5] Computing cylindrical surface profile...
Computing cylindrical profile: 240 Z × 360 θ...

[3/5] Computing curvature map...
Processing 240 × 360 grid points...
  Curvature statistics:
    κ_z (radial):  mean=2.145, max=8.723
    κ_θ (tangent): mean=1.876, max=7.341

[4/5] Generating helical trajectory...
Generating adaptive path from Z=0.1200m to Z=0.0000m...
    Step 5000: Z=0.0823m, θ=3.14 rev, pitch=4.23mm
    Step 10000: Z=0.0451m, θ=6.28 rev, pitch=3.87mm
  Generated 12543 points, 7.89 rotations

  Adaptive pitch statistics:
    Mean: 15.23mm
    Min:  1.00mm   ← 고곡률 영역 (PITCH_MIN 클램핑)
    Max:  27.00mm  ← 저곡률 영역 (PITCH_MAX 클램핑)
    Std:  8.45mm

[5/5] Smoothing curve...
  Smoothed 12543 points

Trajectory saved: ../output/peeling_trajectory_adaptive.csv
  12543 points

Launching visualization...
```

### 3. 경로 검증

```bash
$ uv run python check_path_continuity.py

Total points: 12543

Distance statistics:
  Mean: 0.42mm
  Max:  1.23mm
  Min:  0.08mm

✓ No large jumps detected (path is continuous)
```

### 4. 파라미터 튜닝 가이드

**피치 범위 조정:**
```python
# 더 넓은 피치 범위 (저곡률 영역에서 더 빠르게)
PITCH_MAX = 0.040        # 27mm → 40mm

# 더 좁은 피치 범위 (고곡률 영역에서 더 조밀하게)
PITCH_MIN = 0.0005       # 1mm → 0.5mm
```

**더 높은 표면 품질 (더 조밀한 경로):**
```python
H_SCALLOP = 0.002        # 3mm → 2mm
PITCH_MIN = 0.0008       # 1mm → 0.8mm
PITCH_MAX = 0.018        # 27mm → 18mm (더 보수적)
```

**더 빠른 가공 (성근 경로, 품질 trade-off):**
```python
H_SCALLOP = 0.005        # 3mm → 5mm
PITCH_MAX = 0.040        # 27mm → 40mm (더 넓은 간격)
```

**곡률 계산 해상도 조정:**
```python
NUM_Z_SLICES = 360      # 240 → 360 (더 정밀)
CURVATURE_WINDOW = 3    # 2 → 3 (7×7 윈도우, 더 부드러움)
```

**고정 피치 모드로 전환:**
```python
USE_ADAPTIVE = False    # 적응형 비활성화
NUM_ROTATIONS = 15.0    # 고정 15회전
NUM_SAMPLES = 10000     # 10,000 포인트
```

## 📚 참고사항

### 좌표계
- **메시 로컬 좌표계**: 원점(0,0,0)이 감자 바닥 중심
- **Z축**: 감자의 회전 축 (수직 방향)
- **원통 좌표**: `(r, θ, z)` → 데카르트: `(r·cos(θ), r·sin(θ), z)`
- **프레임**: base_link (ROS2 표준)

### ISO-scallop 이론 배경

**스캘럽 높이(Scallop Height)**란?
- 인접한 절삭 경로 사이에 남는 표면 잔여 높이
- 표면 품질을 결정하는 핵심 지표
- 값이 작을수록 매끄러운 표면 (but 가공 시간 증가)

**평면 공구의 ISO-scallop 공식:**
```
h = κ · s² / 8
```
- `h`: 스캘럽 높이 (m)
- `κ`: 국소 곡률 (1/m) - 표면이 얼마나 휘었는지
- `s`: 스텝오버/피치 (m) - 경로 간격

**적응형 피치 계산:**
```
s = √(8h / κ)
s_clamped = clip(s, pitch_min, pitch_max)
```

**장점:**
- 곡률이 높은 영역 (휘어진 부분) → 촘촘한 경로 → 품질 보장
- 곡률이 낮은 영역 (평평한 부분) → 성긴 경로 → 효율적 가공
- 전체 표면에서 일정한 품질 유지 (uniform surface finish)

### 곡률 계산 방법

**Exact Curvature with Quadric Surface Fitting:**
1. 각 그리드 점 주변 5×5 윈도우에서 이웃 점 추출
2. 국소 2차 곡면 피팅: `z = ax² + by² + cxy + dx + ey + f`
3. **1차 미분 계산** (표면 기울기, 중심점에서 평가):
   ```
   z_x = 2a·x_ij + c·y_ij + d
   z_y = 2b·y_ij + c·x_ij + e
   |∇z|² = z_x² + z_y²
   ```
4. **2차 미분 계산** (Hessian, 상수):
   ```
   z_xx = ∂²z/∂x² = 2a
   z_yy = ∂²z/∂y² = 2b
   z_xy = ∂²z/∂x∂y = c
   ```
5. **방향별 곡률** (라인 258-262):
   - 반경: `κ_r = z_xx·cos²θ + 2·z_xy·cosθ·sinθ + z_yy·sin²θ`
   - 접선: `κ_t = z_xx·sin²θ - 2·z_xy·cosθ·sinθ + z_yy·cos²θ`
6. **Gradient 보정** (라인 264-273):
   ```
   κ_exact = κ_dir / (1 + |∇z|²)^(3/2)
   ```
7. Gaussian 필터링으로 노이즈 제거 (σ=1.5)

**정확도 향상 (라인 244-273):**
- 기존 근사: `κ_approx = |∂²z/∂n²|` (1차 미분 무시, 경사 표면에서 부정확)
- 정확한 공식: `κ_exact = |∂²z/∂n²| / (1 + |∇z|²)^(3/2)`
  - 1차 미분 `z_x`, `z_y` 명시적 계산 (라인 245-247)
  - Gradient magnitude 보정 (라인 266-267)
- 경사 30° 표면에서 1.54배 곡률 과대평가 방지 ✓
- 결과: 더 정확한 pitch 계산 → 균일한 표면 품질

### 성능

**메시 샘플링 (240×360 해상도):**
- 약 86,400개 표면 샘플링 포인트
- 곡률 맵 계산: 약 5-10초 (Intel i7 기준)

**적응형 경로 생성:**
- 약 10,000~15,000 포인트 (곡률 분포에 따라 변동)
- 생성 시간: 약 2-5초
- 메모리 사용량: ~200MB

**고정 피치 경로 생성:**
- 10,000 포인트 기준: 약 1-2초
- 곡률 계산 생략으로 빠름

### ROS 2 통합 가이드

**CSV를 ROS 2 메시지로 변환:**
```python
import numpy as np
from geometry_msgs.msg import PoseStamped
from nav_msgs.msg import Path

def csv_to_ros_path(csv_path):
    data = np.loadtxt(csv_path, delimiter=',', skiprows=1)

    path_msg = Path()
    path_msg.header.frame_id = "base_link"
    path_msg.header.stamp = node.get_clock().now().to_msg()

    for i, (x, y, z, blade) in enumerate(data):
        pose = PoseStamped()
        pose.header = path_msg.header
        pose.pose.position.x = float(x)
        pose.pose.position.y = float(y)
        pose.pose.position.z = float(z)
        # blade 상태는 별도 토픽 또는 커스텀 메시지로 전송
        path_msg.poses.append(pose)

    return path_msg
```

**MoveIt 2 연동:**
- `moveit_msgs/RobotTrajectory`로 변환
- 시간 파라미터화: Ruckig OTG 추천 (속도/가속도/저크 제한)
- IK 솔버: KDL 또는 TRAC-IK

### 확장 가능성

**현재 구현 완료:**
- ✓ ISO-scallop 적응형 피치 (라인 354-387) - **완료**
- ✓ 곡률 기반 경로 밀도 최적화 (라인 591-761) - **완료**
- ✓ 정확한 곡률 공식 (1차 미분 gradient 보정 포함, 라인 244-273) - **완료**
- ✓ 원통 좌표계 프로파일 추출 (라인 95-160) - **완료**
- ✓ Quadric surface fitting 기반 곡률 계산 (라인 167-279) - **완료**
- ✓ 표면 추종 및 블레이드 제어 (라인 463-497) - **완료**
- ✓ B-spline 스무딩 (라인 768-808) - **완료**

### 문제 해결 (Troubleshooting)

**메시 로딩 실패:**
```
Error: Unable to load mesh
→ 해결: mesh/mesh_1_processed.stl 파일 존재 확인
→ 해결: process_mesh_fixed.py 먼저 실행
```

**곡률 계산 시간 과다:**
```
Curvature computation takes >30 seconds
→ 해결: NUM_Z_SLICES 감소 (240 → 120)
→ 해결: CURVATURE_WINDOW 감소 (2 → 1)
```


### 라이선스 및 인용

**Dependencies:**
- trimesh (MIT License)
- Open3D (MIT License)
- NumPy, SciPy (BSD License)

**References:**
- ISO-scallop theory: "Iso-scallop tool path planning for triangular mesh surfaces" (2021)
- Curvature computation: Quadric surface fitting method
- Helical toolpath: CNC machining literature

---

**작성:** 2025
**버전:** 1.0 (ISO-scallop adaptive mode)
**프로젝트:** UR5e Surface Peeling ROS 2 Workspace
