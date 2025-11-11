# Potato Path Generator

감자 껍질 벗기기를 위한 헬리컬 경로 생성기

## 📖 프로젝트 개요

회전하는 감자의 표면을 따라 로봇 매니퓰레이터가 이동할 수 있는 헬리컬 궤적을 생성합니다.
- 감자가 Z축을 기준으로 회전하는 동안 로봇이 나선형 경로를 따라 이동
- STL 메시 기반 표면 프로파일 추출 및 경로 생성
- 원통 좌표계(cylindrical coordinates) 기반 radial profile r(θ, z) 계산

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
# 헬리컬 경로 생성
uv run python src/potato_peeling_path_generator.py
```

## 📁 파일 구조

```
potato/
├── pyproject.toml                      # UV 패키지 설정
├── README.md                           # 이 파일
├── mesh/
│   ├── mesh_1.stl                     # 원본 STL 메시
│   ├── mesh_1_aligned.stl             # 정렬된 메시 (중간 출력)
│   └── mesh_1_processed.stl           # 전처리된 메시 (최종)
├── output/                             # 생성된 경로 출력 디렉토리
│   └── (생성된 CSV, PLY 파일들)
└── src/
    ├── potato_peeling_path_generator.py  # ⭐ 헬리컬 경로 생성기 (메인)
    └── process_mesh_fixed.py             # 메시 전처리 스크립트
```

## ⚙️ 주요 매개변수

`src/potato_peeling_path_generator.py` 내부에서 수정:

```python
# 원통 좌표계 샘플링 해상도
num_z_slices = 120        # Z 방향 슬라이스 개수
num_theta_samples = 360   # 각도(θ) 샘플 개수

# 헬리컬 경로 매개변수
theta_start = 0           # 시작 각도 (rad)
theta_end = 8 * 2*np.pi   # 종료 각도 (8회전)
num_points = 1000         # 경로 포인트 개수
```

## 📊 출력 형식

### 생성되는 파일들

1. **원통 좌표계 프로파일** (CSV)
   - `r(θ, z)`: 각 각도와 높이에서의 반경 값
   - 표면 형상을 원통 좌표계로 표현

2. **헬리컬 경로** (CSV, ROS2 호환)
   ```csv
   x,y,z
   0.0234,0.0156,0.0010
   0.0238,0.0154,0.0012
   ...
   ```
   - `x, y, z`: TCP 위치 (m, 데카르트 좌표계)

3. **시각화 파일** (PLY)
   - Open3D로 확인 가능한 포인트 클라우드

## 🎯 알고리즘 개요

### 1. 메시 전처리
- **단위 정규화**: mm 단위 메시를 m 단위로 변환
- **좌표 정렬**: X축 회전으로 감자의 장축을 Z축에 정렬
- **원점 조정**: 바닥 중심을 원점(0,0,0)으로 설정
- **표면 보정**: 구멍 메우기 및 Laplacian 스무딩

### 2. 원통 좌표계 변환
- 메시를 원통 좌표계 `(r, θ, z)`로 샘플링
- 각 `(θ, z)` 위치에서 표면까지의 반경 `r` 계산
- 2D 스플라인 보간으로 연속 함수 `r(θ, z)` 생성

### 3. 헬리컬 경로 생성
- Z 방향으로 하강하면서 θ 증가 (나선형)
- 각 지점에서 `r(θ, z)`로 표면까지의 거리 계산
- 데카르트 좌표로 변환: `x = r·cos(θ), y = r·sin(θ), z = z`
- Gaussian 필터링으로 궤적 스무딩

### 4. 출력 및 검증
- CSV 형식으로 경로 저장
- PLY 파일로 시각화
- 원본 메시와 경로를 함께 렌더링하여 검증

## 🔧 사용 예시

### 메시 전처리 확인

```python
# src/process_mesh_fixed.py 실행 후 출력 확인
# Final bounds: min=[0, 0, 0], max=[0.08, 0.08, 0.12]
# Z bottom: 0.000000 (should be ~0)  ← 원점이 바닥에 위치
```

### 경로 생성 커스터마이징

```python
# src/potato_peeling_path_generator.py 수정
# 더 조밀한 경로가 필요한 경우:
theta_end = 16 * 2*np.pi  # 16회전
num_points = 2000          # 포인트 2배

# 더 빠른 생성이 필요한 경우:
num_z_slices = 60          # 해상도 절반
num_theta_samples = 180
```

## 📚 참고사항

### 좌표계
- **메시 로컬 좌표계**: 원점(0,0,0)이 감자 바닥 중심
- **Z축**: 감자의 회전 축 (수직 방향)
- **원통 좌표**: `(r, θ, z)` → 데카르트: `(r·cos(θ), r·sin(θ), z)`

### 성능
- 120×360 해상도: 약 43,200개 표면 샘플링 포인트
- 1000 포인트 경로 생성: 약 1-2초 소요 (하드웨어 의존)
- 스플라인 보간으로 메모리 효율적 표현

### 확장 가능성
- ROS2 통합: `geometry_msgs/PoseStamped` 형식으로 변환 가능
- 속도 프로파일: 시간 매개변수화 추가 가능
- 충돌 회피: MoveIt2와 연동하여 경로 검증 가능
