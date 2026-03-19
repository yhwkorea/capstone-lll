# LLL 관련 논문 서베이 (2024–2026)

> 조사 기준: CRYPTO, EUROCRYPT, ASIACRYPT, CCS, TCHES, IEEE S&P, IACR ePrint  
> 조사일: 2026-03-19

---

## 1. 알고리즘 최적화

### Towards a Modern LLL Implementation
- **저자:** Léo Ducas, Ludo N. Pulles 등
- **발표:** ASIACRYPT 2025 / ePrint 2025/774
- **핵심:** BLAS 기반 "BLASter" 구현체. fplll 대비 100x 이상 속도 향상.  
  segmented deep-LLL로 BKZ-15 수준 품질 달성.
- **링크:** https://eprint.iacr.org/2025/774

### Another L Makes It Better? (L4 Algorithm)
- **저자:** Université de Picardie Jules Verne 팀
- **발표:** SIAM ALENEX 2025 / ePrint 2024/1681
- **핵심:** Lagrange 쌍별 리덕션 + LLL 융합 → 평균 16% 짧은 벡터.  
  BKZ-24 전처리 시 최단 벡터 노름 평균 3% 감소.
- **링크:** https://eprint.iacr.org/2024/1681

### A Complete Analysis of the BKZ Lattice Reduction Algorithm
- **발표:** Journal of Cryptology, 2025
- **핵심:** 최초 엄밀한 BKZ 동적 분석. 근사 오라클 일반화, 실행 중 기저 품질 보장.
- **링크:** https://link.springer.com/article/10.1007/s00145-024-09527-0

### A Parameter Study for LLL and BKZ with Application to SVP
- **발표:** arXiv 2025 (2502.05160)
- **핵심:** LLL/BKZ 파라미터 설정이 SVP 성능에 미치는 영향 체계적 분석.
- **링크:** https://arxiv.org/html/2502.05160v1

### Predicting Module-Lattice Reduction
- **저자:** Léo Ducas, Engelberts, de Perthuis
- **발표:** ASIACRYPT 2025 / ePrint 2025/1904
- **핵심:** 모듈 격자에서 BKZ 평균 케이스 분석.  
  원분체 종류에 따라 블록크기 d−1+o(1) 추가 필요 vs 부분지수 속도향상 분기.
- **링크:** https://eprint.iacr.org/2025/1904

---

## 2. 특정 암호 응용 (본 연구와 직접 연관)

### ★ Constant Time Lattice Reduction in Dimension 4 with Application to SQIsign
- **저자:** Otto Hanyecz, Alexander Karenin, Elena Kirshanova, Péter Kutas, Sina Schaeffler
- **발표:** TCHES 2025 / ePrint 2025/027
- **핵심:** 최초 상수 시간 LLL류 알고리즘 (dim-4, 정수 격자).  
  SQIsign 후양자 서명 서브루틴으로 직접 적용.  
  출력 벡터 길이 보장 포함, Minkowski 리덕션 기저 출력 실험 확인.
- **Gap:** 차원 4에 명시적으로 한정됨. 저자들이 dim ≥ 8 확장을 미해결 과제로 언급.
- **링크:** https://eprint.iacr.org/2025/027

### Lattice Reduction via Dense Sublattices: A Cryptanalytic No-Go
- **발표:** IACR Communications in Cryptology 2025 / ePrint 2025/1694
- **핵심:** 밀집 부분격자 기반 리덕션의 한계 분석.  
  k=2,3에서 BKZ보다 약간 낫지만 "Dimensions for Free" BKZ보다 비효율.
- **링크:** https://eprint.iacr.org/2025/1694

---

## 3. LWE/격자 기반 암호 공격

### Solving LWE with Independent Hints about Secret and Errors
- **발표:** ePrint 2025/1128
- **핵심:** LLL 리덕션을 행렬 곱으로 대체하여 격자 구성 복잡도 감소.  
  CRYSTALS-KYBER 공격 가속 실험 검증.
- **링크:** https://eprint.iacr.org/2025/1128

### Refined Modelling of the Primal Attack, and Variants Against Module-LWE
- **발표:** ePrint 2025/2195
- **핵심:** BKZ 시뮬레이터 일반화. Module-LWE Primal Attack 추정기 설계.
- **링크:** https://eprint.iacr.org/2025/2195

---

## 4. RSA 및 RSA 유사 암호 공격

### A New Generalized Lattice Attack Against RSA-Like Cryptosystems
- **발표:** ePrint 2025/1559
- **핵심:** 소형/대형 비밀 지수 통합 일반화 공격 프레임워크.
- **링크:** https://eprint.iacr.org/2025/1559

### Improving RSA Cryptanalysis: Continued Fractions + Coppersmith
- **발표:** ePrint 2025/1281
- **핵심:** 연분수(d < N^0.25) + Coppersmith 격자(d < N^0.292) 결합 분석.
- **링크:** https://eprint.iacr.org/2025/1281

### A note on the analysis of Herrmann-May Lattices for Small Exponent RSA
- **발표:** Scientific Reports, 2025
- **핵심:** Herrmann-May 경계를 d < N^0.292256으로 수정 (기존보다 엄밀히 낮음).
- **링크:** https://www.nature.com/articles/s41598-025-10019-9

### Attacking an RSA-like Cryptosystem Using Continued Fractions and Lattices
- **발표:** ePrint 2025/1739
- **핵심:** RSA 유사 암호시스템에 연분수 + LLL 격자 결합 적용.
- **링크:** https://eprint.iacr.org/2025/1739

---

## 동향 요약

| 분야 | 주요 동향 |
|------|-----------|
| LLL 구현 최적화 | BLAS 기반 BLASter (100x), segmented deep-LLL |
| BKZ 이론 분석 | 최초 엄밀 동적 분석 (J. Cryptology 2025) |
| 모듈 격자 | 원분체 종류별 module-BKZ 블록크기 분기 |
| **사이드채널 방어** | **SQIsign용 최초 상수 시간 dim-4 LLL (TCHES 2025) ← 본 연구 출발점** |
| LWE 공격 | LLL → 행렬 곱 대체로 Kyber 공격 가속 |
| RSA 공격 | Coppersmith+연분수 결합, Herrmann-May 경계 수정 |
