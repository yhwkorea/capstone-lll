# Capstone Research: Constant-Time LLL Lattice Reduction

> **Goal:** Extend constant-time LLL-type lattice reduction from dimension 4 to dimension 8/16,  
> and integrate it into post-quantum signature schemes (FALCON, SQIsign).

---

## Research Question

The state-of-the-art constant-time lattice reduction (TCHES 2025, [ePrint 2025/027](https://eprint.iacr.org/2025/027))  
is explicitly limited to **dimension 4**. Higher-dimensional post-quantum schemes such as FALCON and NTRU  
require fast, side-channel-resistant basis reduction in larger dimensions.

**Can we build a constant-time LLL/Lagrange-type reduction for dim ≥ 8,  
and demonstrate its practical use inside a post-quantum signature scheme?**

---

## Project Structure

```
capstone-lll/
├── docs/
│   ├── papers.md       ← Survey of related work (2024–2026)
│   └── roadmap.md      ← 12-month implementation plan
├── src/                ← C++ / Python implementation (to be added)
└── benchmarks/         ← Timing measurements, side-channel tests
```

---

## Key Papers

| Paper | Venue | Relevance |
|-------|-------|-----------|
| [Constant Time Lattice Reduction in Dim 4 (2025/027)](https://eprint.iacr.org/2025/027) | TCHES 2025 | **Direct base — the gap we extend** |
| [Towards a Modern LLL Implementation (2025/774)](https://eprint.iacr.org/2025/774) | ASIACRYPT 2025 | BLASter: 100x faster LLL, reference impl |
| [Another L Makes It Better (2024/1681)](https://eprint.iacr.org/2024/1681) | ALENEX 2025 | L4 algorithm: 16% shorter vectors |
| [Complete Analysis of BKZ (J. Cryptology 2025)](https://link.springer.com/article/10.1007/s00145-024-09527-0) | J. Cryptology | First rigorous dynamic BKZ analysis |
| [Predicting Module-Lattice Reduction (2025/1904)](https://eprint.iacr.org/2025/1904) | ASIACRYPT 2025 | Module-BKZ block size bounds |

Full survey → [`docs/papers.md`](docs/papers.md)

---

## Timeline

| Phase | Period | Milestone |
|-------|--------|-----------|
| Background | Month 1–2 | Reproduce dim-4 constant-time LLL |
| Design | Month 3–5 | dim-8 branch-free reduction design |
| Integration | Month 6–9 | Plug into FALCON subroutine |
| Evaluation | Month 10–11 | Timing channel measurement, benchmarks |
| Writing | Month 12 | Paper draft (target: TCHES / CCS) |

Full roadmap → [`docs/roadmap.md`](docs/roadmap.md)

---

## Target Venues

- **TCHES** (IACR Transactions on Cryptographic Hardware and Embedded Systems) — implementation-friendly
- **CCS** (ACM Conference on Computer and Communications Security)
- **EUROCRYPT / ASIACRYPT** — if theoretical bound improvement is achieved
