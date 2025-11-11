# Potato Path Generator

ISO-scallop 곡률 적응형 헬리컬 경로 생성기

## 🚀 빠른 시작

### 1. 의존성 설치

```bash
# uv로 의존성 설치
uv sync
```

### 2. 메시 전처리 (한 번만)

```bash
# STL 메시 전처리 (단위 변환, 정렬, 스무딩)
uv run python src/process_mesh_fixed.py
```

### 3. 경로 생성 + 시각화 (올인원)

```bash
# 경로 생성부터 Open3D 시각화까지 한 번에 실행
uv run python main.py
```

### 개별 실행 (선택)

```bash
# 경로만 생성
uv run python run_proper_helical.py

# 기존 경로 시각화만
uv run python visualize.py
uv run python visualize.py output/helical_path.csv
```

## 📁 파일 구조

```
potato/
├── pyproject.toml              # UV 패키지 설정
├── README.md                   # 이 파일
├── main.py                     # ⭐ 올인원 실행 (경로 생성 + 시각화)
├── run_proper_helical.py       # 경로 생성만
├── visualize.py                # 시각화만
├── mesh/
│   ├── mesh_1.stl             # 원본 메시
│   └── mesh_1_processed.stl   # 전처리된 메시
├── output/
│   ├── helical_path.csv       # 생성된 경로 (ROS2 형식)
│   ├── helical_path.ply       # 경로 시각화용
│   └── helical_diagnostics.csv # 진단 정보
└── src/
    ├── process_mesh_fixed.py  # 메시 전처리
    ├── mesh_utils.py          # 메시 유틸리티
    ├── curvature_fit.py       # 곡률 추정
    ├── blade_contact_simple.py # 접촉 검색
    ├── helical_path_generator.py # 메인 생성기
    ├── export_csv.py          # CSV 내보내기
    └── README_KR.md           # 상세 문서 (한글)
```

## ⚙️ 주요 매개변수

`run_proper_helical.py`에서 수정:

```python
params = HelicalPathParams(
    W_tool=0.03,              # 칼날 너비 (m)
    h_scallop=0.0003,         # 목표 스캘럽 높이 (m)
    theta_turns=8.0,          # 회전 수
    dtheta=np.deg2rad(1.5),   # 각도 증분 (rad)
    s_min=0.003,              # 최소 스텝오버 (m)
    s_max=0.015,              # 최대 스텝오버 (m)
    z_bounds=(0.001, 0.062),  # Z 범위 (m)
)
```

## 📊 출력 형식

### CSV (ROS2 호환)
```
frame_id,world
x,y,z,qx,qy,qz,qw
0.029183,0.016787,0.021973,0.500000,0.500000,0.500000,0.500000
...
```

- `x, y, z`: TCP 위치 (m)
- `qx, qy, qz, qw`: 방향 (quaternion, scalar-last)

## 🎯 알고리즘

1. **ISO-Scallop 이론** - 곡률 기반 스텝오버 계산
2. **곡률 추정** - Quadric fitting으로 k₁, k₂ 계산
3. **적응형 피치** - 스텝오버에서 피치 자동 매핑
4. **접촉 검색** - 원통 좌표계 기반 최근접 점

상세 내용: `src/README_KR.md` 참조
