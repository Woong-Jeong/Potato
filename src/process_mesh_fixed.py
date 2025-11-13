#!/usr/bin/env python3
"""
STL 메시 처리 스크립트 (수정됨):
0. 단위 변환 (mm → m)
1. X축 기준 -90도 회전
2. 원점을 메시 가장 아래 중심으로 이동
3. 구멍 메우기 (hole filling)
4. 스무딩 (smoothing)
"""

import numpy as np
import trimesh
import os

def process_mesh(input_path, output_path):
    """메시를 처리하고 저장하는 함수"""

    print(f"Loading mesh from: {input_path}")

    # 1. 메시 로드 (trimesh 사용)
    mesh = trimesh.load_mesh(input_path)
    print(f"Original mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print(f"Original bounds: min={mesh.bounds[0]}, max={mesh.bounds[1]}")

    # ★ 0. 단위 변환 (mm → m)
    print("\n=== Step 0: Converting units (mm → m) ===")
    # 경계가 30m 이상이면 밀리미터 단위로 간주
    if np.max(np.abs(mesh.bounds)) > 10.0:
        print("  Detected millimeter units, converting to meters...")
        mesh.vertices *= 0.001  # mm → m
        print(f"  After conversion bounds: min={mesh.bounds[0]}, max={mesh.bounds[1]}")
    else:
        print("  Units appear to be in meters already")

    # 2. X축 기준 -90도 회전
    print("\n=== Step 1: Rotating mesh by -90 degrees around X axis ===")
    rotation_matrix = trimesh.transformations.rotation_matrix(
        np.radians(-90),  # -90도를 라디안으로 변환
        [1, 0, 0],  # X축
        [0, 0, 0]   # 원점 기준
    )
    mesh.apply_transform(rotation_matrix)
    print(f"After rotation bounds: min={mesh.bounds[0]}, max={mesh.bounds[1]}")

    # 3. 원점을 메시 가장 아래 중심으로 이동
    print("\n=== Step 2: Moving origin to bottom center ===")
    bounds = mesh.bounds

    # X, Y는 중심으로, Z는 최소값(가장 아래)으로
    center_x = (bounds[0][0] + bounds[1][0]) / 2
    center_y = (bounds[0][1] + bounds[1][1]) / 2
    bottom_z = bounds[0][2]

    # 이동 벡터 계산
    translation_vector = np.array([-center_x, -center_y, -bottom_z])
    mesh.apply_translation(translation_vector)

    print(f"Translation vector: {translation_vector}")
    print(f"After translation bounds: min={mesh.bounds[0]}, max={mesh.bounds[1]}")

    # 4. 구멍 메우기 (trimesh 사용 - 강화)
    print("\n=== Step 3: Filling holes ===")
    print(f"Before hole filling: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")

    # watertight 확인
    print(f"Is watertight before filling: {mesh.is_watertight}")

    # trimesh의 fill_holes 사용 (더 적극적으로)
    try:
        # 먼저 작은 구멍들 메우기
        mesh.fill_holes()

        # 여전히 watertight가 아니면 더 강력한 방법 시도
        if not mesh.is_watertight:
            print("Mesh still has holes, trying more aggressive hole filling...")
            # 메시 수리 시도
            trimesh.repair.fill_holes(mesh)
            trimesh.repair.fix_normals(mesh)
            trimesh.repair.fix_inversion(mesh)

        print(f"After hole filling: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
        print(f"Is watertight after filling: {mesh.is_watertight}")
    except Exception as e:
        print(f"Warning: Hole filling encountered an issue: {e}")

    # 5. 스무딩 (강화된 Laplacian smoothing)
    print("\n=== Step 4: Smoothing mesh (강화) ===")

    # Trimesh의 smoothing 기능 사용
    try:
        # 더 강력한 Laplacian smoothing 적용 (반복 횟수와 lambda 값 증가)
        print("Applying aggressive Laplacian smoothing...")
        trimesh.smoothing.filter_laplacian(mesh, iterations=20, lamb=0.7)
        print(f"After first smoothing pass: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")

        # 추가 스무딩 패스
        trimesh.smoothing.filter_laplacian(mesh, iterations=10, lamb=0.6)
        print(f"After smoothing: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    except Exception as e:
        print(f"Note: Using alternative smoothing method: {e}")
        # 대체 방법: vertex를 이웃의 평균으로 이동 (더 많은 반복)
        print("Applying manual Laplacian smoothing with more iterations...")
        for iteration in range(10):
            vertex_adjacency = mesh.vertex_neighbors
            new_vertices = mesh.vertices.copy()
            for i, neighbors in enumerate(vertex_adjacency):
                if len(neighbors) > 0:
                    neighbor_positions = mesh.vertices[neighbors]
                    # 더 강한 스무딩 (0.5 대신 0.4)
                    new_vertices[i] = 0.4 * mesh.vertices[i] + 0.6 * np.mean(neighbor_positions, axis=0)
            mesh.vertices = new_vertices
            if (iteration + 1) % 3 == 0:
                print(f"  Smoothing iteration {iteration + 1}/10 completed")
        print(f"After smoothing: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")

    # 6. 최종 메시 저장
    mesh.export(output_path)
    print(f"\n=== Final mesh saved to: {output_path} ===")

    # 최종 통계
    print(f"Final mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
    print(f"Final bounds: min={mesh.bounds[0]}, max={mesh.bounds[1]}")

    # 원점이 실제로 메시 가장 아래 중심에 있는지 확인
    print(f"\nOrigin position check:")
    print(f"  X center: {(mesh.bounds[0][0] + mesh.bounds[1][0]) / 2:.6f} (should be ~0)")
    print(f"  Y center: {(mesh.bounds[0][1] + mesh.bounds[1][1]) / 2:.6f} (should be ~0)")
    print(f"  Z bottom: {mesh.bounds[0][2]:.6f} (should be ~0)")

    return mesh


if __name__ == "__main__":
    # 현재 스크립트 디렉토리 기준으로 경로 설정
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # 입력 및 출력 파일 경로
    input_mesh_path = os.path.join(script_dir, "..", "mesh", "mesh_1.stl")
    output_mesh_path = os.path.join(script_dir, "..", "mesh", "mesh_1_processed.stl")

    print(f"Input:  {input_mesh_path}")
    print(f"Output: {output_mesh_path}")

    # 메시 처리 실행
    processed_mesh = process_mesh(input_mesh_path, output_mesh_path)

    print("\n" + "="*60)
    print("Processing completed successfully!")
    print("="*60)
