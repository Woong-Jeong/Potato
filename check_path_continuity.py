#!/usr/bin/env python3
"""경로 연속성 진단 스크립트"""
import numpy as np

# CSV 파일 로드
try:
    data = np.loadtxt('output/peeling_trajectory_adaptive.csv', delimiter=',', skiprows=1)
    print(f"Total points: {len(data)}")
    
    # 연속된 점들 사이 거리 계산
    positions = data[:, :3]
    distances = np.linalg.norm(np.diff(positions, axis=0), axis=1)
    
    print(f"\nDistance statistics:")
    print(f"  Mean: {np.mean(distances)*1000:.2f}mm")
    print(f"  Max:  {np.max(distances)*1000:.2f}mm")
    print(f"  Min:  {np.min(distances)*1000:.2f}mm")
    
    # 큰 점프 찾기 (평균의 10배 이상)
    threshold = np.mean(distances) * 10
    jumps = np.where(distances > threshold)[0]
    
    if len(jumps) > 0:
        print(f"\n⚠️  Found {len(jumps)} large jumps (>{threshold*1000:.1f}mm):")
        for idx in jumps[:5]:  # 처음 5개만 표시
            print(f"    Point {idx} → {idx+1}: {distances[idx]*1000:.2f}mm")
    else:
        print(f"\n✓ No large jumps detected (path is continuous)")
        
except FileNotFoundError:
    print("Error: output/peeling_trajectory_adaptive.csv not found")
    print("Run the generator first!")
