/**
 * test_ct_lll4.cpp — Phase 1 검증
 *
 * 컴파일:
 *   g++ -O2 -std=c++17 -I../src/dim4 test_ct_lll4.cpp -o test_ct_lll4
 *
 * 실행:
 *   ./test_ct_lll4
 */

#include "ct_lll4.hpp"
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>

// ---------------------------------------------------------------------------
// Simple LCG for reproducible random bases (no dependency on libc rand quirks)
// ---------------------------------------------------------------------------
static uint64_t lcg_state = 12345678901234ULL;
static int64_t lcg_range(int64_t lo, int64_t hi) {
    lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
    uint64_t r = lcg_state >> 33;
    return lo + (int64_t)(r % (uint64_t)(hi - lo + 1));
}

static Mat4 random_basis(int bits) {
    Mat4 B;
    int64_t bound = (int64_t)1 << bits;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            B[i][j] = lcg_range(-bound, bound);
    return B;
}

// ---------------------------------------------------------------------------
// LLL size-reduced check (|mu[i][j]| <= 1/2 for all j < i)
// ---------------------------------------------------------------------------
static bool check_size_reduced(const Mat4 &B) {
    // Compute GSO directly from B for verification (not CT, just for testing)
    double gs[4][4] = {};
    double mu[4][4] = {};
    for (int i = 0; i < 4; i++) {
        for (int k = 0; k < 4; k++) gs[i][k] = (double)B[i][k];
        for (int j = 0; j < i; j++) {
            double dot = 0, gsj_sq = 0;
            for (int k = 0; k < 4; k++) { dot += gs[i][k] * gs[j][k]; gsj_sq += gs[j][k] * gs[j][k]; }
            mu[i][j] = dot / gsj_sq;
            for (int k = 0; k < 4; k++) gs[i][k] -= mu[i][j] * gs[j][k];
        }
    }
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < i; j++)
            if (std::abs(mu[i][j]) > 0.5 + 1e-6) return false;
    return true;
}

// ---------------------------------------------------------------------------
// Test 1: correctness on known small basis
// ---------------------------------------------------------------------------
static bool test_known_basis() {
    // Basis with large off-diagonal entries — should reduce cleanly
    const int64_t data[16] = {
        137, 312, -89,  44,
        -21,  15, 300, -77,
        180, -33,  67, 200,
         55,  88, -44, 111
    };
    Mat4 B = make_basis(data);
    ct_reduce_dim4(B, 10);

    bool ok = check_size_reduced(B);
    printf("[test_known_basis] %s\n", ok ? "PASS" : "FAIL");
    if (!ok) {
        printf("  Reduced basis (rows):\n");
        for (int i = 0; i < 4; i++)
            printf("    [%lld, %lld, %lld, %lld]\n",
                   (long long)B[i][0], (long long)B[i][1],
                   (long long)B[i][2], (long long)B[i][3]);
    }
    return ok;
}

// ---------------------------------------------------------------------------
// Test 2: random correctness (N random bases)
// ---------------------------------------------------------------------------
static bool test_random_correctness(int N = 200, int bits = 20) {
    int fail = 0;
    for (int t = 0; t < N; t++) {
        Mat4 B = random_basis(bits);
        ct_reduce_dim4(B, bits);
        if (!check_size_reduced(B)) fail++;
    }
    bool ok = (fail == 0);
    printf("[test_random_correctness] %s  (N=%d, bits=%d, fail=%d)\n",
           ok ? "PASS" : "FAIL", N, bits, fail);
    return ok;
}

// ---------------------------------------------------------------------------
// Test 3: timing distribution (CV and max/min ratio)
// ---------------------------------------------------------------------------
static void test_timing(int N = 1000, int bits = 30) {
    std::vector<double> times(N);
    for (int t = 0; t < N; t++) {
        Mat4 B = random_basis(bits);
        auto t0 = std::chrono::high_resolution_clock::now();
        ct_reduce_dim4(B, bits);
        auto t1 = std::chrono::high_resolution_clock::now();
        times[t] = std::chrono::duration<double, std::nano>(t1 - t0).count();
    }

    double mean = 0;
    for (double x : times) mean += x;
    mean /= N;

    double var = 0;
    for (double x : times) var += (x - mean) * (x - mean);
    double std_dev = std::sqrt(var / N);

    double mn = *std::min_element(times.begin(), times.end());
    double mx = *std::max_element(times.begin(), times.end());

    printf("[test_timing] n=%d  bits=%d\n", N, bits);
    printf("  mean=%.0f ns  std=%.0f ns  min=%.0f  max=%.0f\n", mean, std_dev, mn, mx);
    printf("  CV=%.1f%%  max/min=%.2fx\n", std_dev / mean * 100.0, mx / mn);
    printf("  (fpylll baseline: CV=32.7%%  max/min=8.26x)\n");
    printf("  Target for CT impl: CV<5%%  max/min<2.0\n");
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main() {
    printf("=== ct_reduce_dim4 — Phase 1 Tests ===\n\n");

    bool ok = true;
    ok &= test_known_basis();
    ok &= test_random_correctness();
    printf("\n");
    test_timing();
    printf("\nOverall: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
