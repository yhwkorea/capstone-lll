"""
test_fpylll_baseline.py — Phase 1 검증 스크립트

목적:
  1. fpylll로 dim-4 LLL 실행 → 기준 벡터 길이 측정
  2. 입력에 따른 타이밍 분산 측정 (constant-time 검증 사전 준비)
  3. 나중에 C++ 구현과 결과 비교용 레퍼런스 생성

실행:
  python3 tests/test_fpylll_baseline.py
"""

import time
import statistics
import random
from fpylll import IntegerMatrix, LLL, GSO

DIMENSION = 4
NUM_TRIALS = 1000
BIT_SIZE = 30  # 입력 격자 계수 비트 수


def random_basis(dim: int, bits: int) -> IntegerMatrix:
    A = IntegerMatrix(dim, dim)
    for i in range(dim):
        for j in range(dim):
            A[i][j] = random.randint(-(2**bits), 2**bits)
    return A


def vector_norm_sq(v) -> int:
    return sum(x * x for x in v)


def run_lll_and_measure(A: IntegerMatrix):
    """LLL 실행 후 (시간, 첫 벡터 노름²) 반환"""
    B = IntegerMatrix(A)  # copy
    t0 = time.perf_counter_ns()
    LLL.reduction(B)
    t1 = time.perf_counter_ns()
    norm_sq = vector_norm_sq(B[0])
    return t1 - t0, norm_sq


def test_correctness():
    """LLL 후 기저가 size-reduced인지 기본 확인"""
    print("=== Correctness check ===")
    A = random_basis(DIMENSION, BIT_SIZE)
    B = IntegerMatrix(A)
    LLL.reduction(B)

    M = GSO.Mat(B)
    M.update_gso()

    ok = True
    for i in range(DIMENSION):
        for j in range(i):
            mu = M.get_mu(i, j)
            if abs(mu) > 0.5 + 1e-9:
                print(f"  FAIL: mu[{i}][{j}] = {mu:.4f} > 0.5")
                ok = False
    if ok:
        print("  PASS: size-reduced condition satisfied")

    # Lovász condition (delta=0.99)
    delta = 0.99
    for i in range(1, DIMENSION):
        lhs = M.get_r(i, i)
        rhs = (delta - M.get_mu(i, i-1)**2) * M.get_r(i-1, i-1)
        if lhs < rhs - 1e-6:
            print(f"  FAIL: Lovász at i={i}: {lhs:.2f} < {rhs:.2f}")
            ok = False
    if ok:
        print("  PASS: Lovász condition satisfied (delta=0.99)")
    print()


def test_timing_distribution():
    """
    타이밍 분포 측정 — 상수 시간 구현 검증의 사전 기준선.
    C++ ct_reduce_dim4 구현 후 이 분포와 비교한다.

    이상적인 constant-time 구현은:
      - std_dev / mean < 0.05 (5% 이내)
      - 최대값 / 최솟값 < 2.0
    """
    print(f"=== Timing distribution (n={NUM_TRIALS}, dim={DIMENSION}, bits={BIT_SIZE}) ===")

    timings = []
    norms = []
    for _ in range(NUM_TRIALS):
        A = random_basis(DIMENSION, BIT_SIZE)
        elapsed, norm_sq = run_lll_and_measure(A)
        timings.append(elapsed)
        norms.append(norm_sq)

    mean_t = statistics.mean(timings)
    std_t = statistics.stdev(timings)
    min_t = min(timings)
    max_t = max(timings)

    print(f"  Time (ns): mean={mean_t:.0f}  std={std_t:.0f}  min={min_t}  max={max_t}")
    print(f"  Coefficient of variation: {std_t/mean_t*100:.1f}%")
    print(f"  Max/Min ratio: {max_t/min_t:.2f}x")
    print()

    mean_norm = statistics.mean(norms)
    print(f"  First vector norm² mean: {mean_norm:.0f}")
    print()

    # 기준선 저장 (C++ 구현 비교용)
    print("  [baseline] Save these numbers to compare with C++ ct_reduce_dim4")
    print(f"  mean_ns={mean_t:.0f}  std_ns={std_t:.0f}  mean_norm_sq={mean_norm:.0f}")


if __name__ == "__main__":
    random.seed(42)
    test_correctness()
    test_timing_distribution()
