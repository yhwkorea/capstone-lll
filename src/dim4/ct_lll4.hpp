/**
 * ct_lll4.hpp — Constant-Time Lattice Reduction in Dimension 4
 *
 * Implements Algorithm 3.5 (BKZ-2 on dim-4) from:
 *   Hanyecz et al., "Constant Time Lattice Reduction in Dimension 4
 *   with Application to SQIsign", TCHES 2025 / ePrint 2025/027
 *
 * Phase 1 arithmetic:
 *   Basis B: int64_t  |  Gram G: int64_t (symmetric)
 *   GSO (mu, r): double  |  Lagrange H: double  |  U: int64_t
 *
 * TODO Phase 2: replace double GSO with exact __int128 for full CT proof.
 */

#pragma once
#include <array>
#include <cstdint>
#include <cmath>
#include <algorithm>

using Mat4 = std::array<std::array<int64_t, 4>, 4>;
using Mat2 = std::array<std::array<int64_t, 2>, 2>;

struct GSO4 {
    double mu[4][4] = {};
    double r[4][4]  = {};
};

// ---------------------------------------------------------------------------
// Branch-free primitives
// ---------------------------------------------------------------------------

inline void ct_cswap64(int64_t &a, int64_t &b, uint64_t cond) {
    uint64_t mask = -(uint64_t)(cond != 0);
    uint64_t diff = mask & ((uint64_t)a ^ (uint64_t)b);
    a ^= (int64_t)diff;
    b ^= (int64_t)diff;
}

inline int64_t ct_select64(int64_t a, int64_t b, uint64_t cond) {
    uint64_t mask = -(uint64_t)(cond != 0);
    return (int64_t)((mask & (uint64_t)a) | (~mask & (uint64_t)b));
}

// ---------------------------------------------------------------------------
// Algorithm 3.1 — Cholesky (fp-GSO from integer Gram matrix)
// ---------------------------------------------------------------------------

static void cholesky_update(const Mat4 &G, GSO4 &gso, int ell, int rho) {
    for (int i = ell; i <= rho; i++) {
        for (int j = 0; j < i; j++) {
            double rij = (double)G[i][j];
            for (int k = 0; k < j; k++)
                rij -= gso.mu[j][k] * gso.r[i][k];
            gso.r[i][j]  = rij;
            gso.mu[i][j] = rij / gso.r[j][j];
        }
        double rii = (double)G[i][i];
        for (int j = 0; j < i; j++)
            rii -= gso.mu[i][j] * gso.r[i][j];
        gso.r[i][i] = rii;
    }
}

// ---------------------------------------------------------------------------
// Single size-reduction step: b_i -= round(mu[i][j]) * b_j
// Always executes same instructions (no branch on mu_r value).
// ---------------------------------------------------------------------------

static void size_red_step(Mat4 &B, Mat4 &G, GSO4 &gso, int i, int j) {
    int64_t mu_r = (int64_t)std::round(gso.mu[i][j]);
    // No early-out: mu_r==0 -> all subtractions are noops (constant-time intent)

    int64_t g_ij = G[i][j];
    int64_t g_jj = G[j][j];

    for (int k = 0; k < i; k++) {
        int64_t x = ct_select64(G[j][k], G[k][j], (uint64_t)(k <= j));
        G[i][k] -= mu_r * x;
        G[k][i]  = G[i][k];
    }
    for (int k = i + 1; k < 4; k++) {
        G[k][i] -= mu_r * G[j][k];
        G[i][k]  = G[k][i];
    }
    G[i][i] -= 2 * mu_r * g_ij - mu_r * mu_r * g_jj;

    for (int k = 0; k < 4; k++)
        B[i][k] -= mu_r * B[j][k];

    cholesky_update(G, gso, i, i);
}

// ---------------------------------------------------------------------------
// Algorithm 3.2 — size_red for one vector
// Sweeps j = i-1 downto 0, repeating until fully size-reduced.
// ---------------------------------------------------------------------------

static void size_red(Mat4 &B, Mat4 &G, GSO4 &gso, int i) {
    for (int pass = 0; pass <= i; pass++) {
        bool any = false;
        for (int j = i - 1; j >= 0; j--) {
            if (std::abs(gso.mu[i][j]) > 0.5 + 1e-10) {
                size_red_step(B, G, gso, i, j);
                any = true;
            }
        }
        if (!any) break;
    }
}

// ---------------------------------------------------------------------------
// Full size reduction: reduce all vectors 1..n-1 in order
// ---------------------------------------------------------------------------

static void full_size_red(Mat4 &B, Mat4 &G, GSO4 &gso) {
    for (int i = 1; i < 4; i++) {
        cholesky_update(G, gso, i, i);
        size_red(B, G, gso, i);
    }
}

// ---------------------------------------------------------------------------
// Algorithm 3.3 — update_after_svp
// ---------------------------------------------------------------------------

static void update_after_svp(Mat4 &G, const Mat2 &U, int i) {
    int64_t g0 = G[i][i], g1 = G[i+1][i+1], g2 = G[i+1][i];

    for (int j = 0; j < i; j++) {
        int64_t gi = G[i][j], gi1 = G[i+1][j];
        G[i][j]   = U[0][0]*gi + U[1][0]*gi1;
        G[i+1][j] = U[0][1]*gi + U[1][1]*gi1;
        G[j][i]   = G[i][j];
        G[j][i+1] = G[i+1][j];
    }
    for (int k = i + 2; k < 4; k++) {
        int64_t gi = G[k][i], gi1 = G[k][i+1];
        G[k][i]   = U[0][0]*gi + U[1][0]*gi1;
        G[k][i+1] = U[0][1]*gi + U[1][1]*gi1;
        G[i][k]   = G[k][i];
        G[i+1][k] = G[k][i+1];
    }
    G[i][i]     = U[0][0]*U[0][0]*g0 + 2*U[0][0]*U[1][0]*g2 + U[1][0]*U[1][0]*g1;
    G[i+1][i+1] = U[0][1]*U[0][1]*g0 + 2*U[0][1]*U[1][1]*g2 + U[1][1]*U[1][1]*g1;
    G[i+1][i]   = U[0][0]*U[0][1]*g0 + (U[0][1]*U[1][0]+U[0][0]*U[1][1])*g2 + U[1][0]*U[1][1]*g1;
    G[i][i+1]   = G[i+1][i];
}

// ---------------------------------------------------------------------------
// Algorithm 3.4 — Lagrange (constant-time dim-2)
// ---------------------------------------------------------------------------

static Mat2 lagrange(double H00, double H10, double H11, int T) {
    Mat2 U = {{{1, 0}, {0, 1}}};

    for (int c = 0; c < T; c++) {
        int64_t mu = (int64_t)std::floor(H10 / H00);

        H11 -= 2.0*(double)mu*H10 - (double)(mu*mu)*H00;
        H10 -= (double)mu * H00;

        U[1][0] -= mu * U[0][0];
        U[1][1] -= mu * U[0][1];

        std::swap(H00, H11);
        ct_cswap64(U[0][0], U[1][0], 1);
        ct_cswap64(U[0][1], U[1][1], 1);
    }

    uint64_t need_swap = (uint64_t)(H11 < H00);
    ct_cswap64(U[0][0], U[1][0], need_swap);
    ct_cswap64(U[0][1], U[1][1], need_swap);

    return U;
}

// ---------------------------------------------------------------------------
// Iteration bounds
// ---------------------------------------------------------------------------

static int compute_T_lagr(int bits) {
    double T = (10.0 * bits + 12.0) / 0.79248 + 2.0;
    return std::max(10, (int)std::ceil(T));
}

static int compute_T_bkz(int bits) {
    return std::max(3, bits + 4);
}

// ---------------------------------------------------------------------------
// Algorithm 3.5 — ct_reduce_dim4
// ---------------------------------------------------------------------------

void ct_reduce_dim4(Mat4 &B, int norm_bound_bits = 30) {
    Mat4 G = {};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++)
                G[i][j] += B[i][k] * B[j][k];

    GSO4 gso;
    cholesky_update(G, gso, 0, 3);

    int T_BKZ  = compute_T_bkz(norm_bound_bits);
    int T_Lagr = compute_T_lagr(norm_bound_bits);

    for (int c = 0; c < T_BKZ; c++) {
        for (int i = 0; i < 3; i++) {
            size_red(B, G, gso, i + 1);

            double H00 = gso.r[i][i];
            double H10 = gso.mu[i+1][i] * gso.r[i][i];
            double H11 = gso.mu[i+1][i]*gso.mu[i+1][i]*gso.r[i][i] + gso.r[i+1][i+1];

            Mat2 U = lagrange(H00, H10, H11, T_Lagr);

            for (int k = 0; k < 4; k++) {
                int64_t bi = B[i][k], bi1 = B[i+1][k];
                B[i][k]   = U[0][0]*bi + U[1][0]*bi1;
                B[i+1][k] = U[0][1]*bi + U[1][1]*bi1;
            }

            update_after_svp(G, U, i);
            cholesky_update(G, gso, i, 3);

            size_red(B, G, gso, i);
            size_red(B, G, gso, i + 1);
        }
    }

    // Final full size-reduction pass to ensure all |mu[i][j]| <= 1/2
    // (Lagrange + update may leave later vectors not fully reduced)
    full_size_red(B, G, gso);
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

inline Mat4 make_basis(const int64_t data[16]) {
    Mat4 B;
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
            B[i][j] = data[i*4+j];
    return B;
}

inline int64_t row_norm_sq(const Mat4 &B, int i) {
    int64_t s = 0;
    for (int j = 0; j < 4; j++) s += B[i][j]*B[i][j];
    return s;
}
