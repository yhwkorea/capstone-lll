# 12개월 연구 로드맵

## 연구 목표

TCHES 2025 (ePrint 2025/027)에서 차원 4에 한정된 상수 시간 격자 리덕션을  
차원 8 → 16으로 확장하고, FALCON 서명 체계 서브루틴에 통합한다.

---

## Phase 1 — 배경 학습 및 재현 (Month 1–2)

### 목표
- LLL/Lagrange 리덕션 이론 습득
- dim-4 상수 시간 LLL (2025/027) 직접 재현

### 세부 작업
- [ ] Nguyen-Stehlé floating-point LLL 논문 정독
- [ ] fplll 라이브러리 빌드 및 코드 분석
- [ ] ePrint 2025/027 알고리즘 C++ 재현 (dim-4)
- [ ] 타이밍 측정 환경 구성 (perf, clock_gettime, constant-time 검증 도구)
- [ ] Minkowski 리덕션 조건 검증 테스트 작성

### 체크포인트
- dim-4 상수 시간 구현체가 논문 표의 벤치마크 ±10% 이내 재현

---

## Phase 2 — dim-8 설계 및 구현 (Month 3–5)

### 목표
- 분기(branch) 없는 dim-8 Lagrange-type 리덕션 알고리즘 설계
- 상수 시간 조건 유지하면서 리덕션 품질 보장

### 핵심 도전
- dim-4에서 사용한 exhaustive case analysis가 dim-8에서 폭발적으로 증가
- SIMD/브랜치리스 비교 연산으로 대체 전략 필요
- Gram-Schmidt 계수 업데이트를 정수/고정소수점으로 관리

### 세부 작업
- [ ] dim-4 → dim-8 케이스 분석 확장 설계 문서 작성
- [ ] 브랜치리스 swap/conditional move 연산 라이브러리 구현
- [ ] dim-8 리덕션 C++ 구현 (초안)
- [ ] 단위 테스트: 임의 격자 10,000개에 대해 리덕션 품질 검증
- [ ] 타이밍 균일성 검증 (입력 의존 분기 없음 확인)

### 체크포인트
- dim-8 구현체가 일반 LLL 대비 리덕션 품질 동등 수준
- Valgrind memcheck / ctgrind로 상수 시간 조건 통과

---

## Phase 3 — FALCON 서브루틴 통합 (Month 6–9)

### 목표
- FALCON 참조 구현(C)에 dim-8 상수 시간 리덕션 교체 적용
- 서명 생성/검증 정확도 유지 확인

### 세부 작업
- [ ] FALCON 참조 구현 코드 분석 (ntru_solve, LDL 분해 부분)
- [ ] 격자 리덕션 호출 지점 식별 및 교체 인터페이스 설계
- [ ] 통합 빌드 및 KAT(Known Answer Test) 통과 확인
- [ ] FALCON 파라미터 셋 512, 1024 모두 테스트
- [ ] dim-16 확장 프로토타입 (선택, 시간 여유 시)

### 체크포인트
- FALCON KAT 100% 통과
- 서명 생성 속도 기존 대비 측정 (느려도 무방, 상수 시간이 목표)

---

## Phase 4 — 사이드채널 평가 및 벤치마크 (Month 10–11)

### 목표
- 타이밍 사이드채널 부재 실험적 검증
- 성능 벤치마크 (dim, 격자 크기별)

### 세부 작업
- [ ] dudect (leakage detection tool) 적용 — t-test 기반 타이밍 누수 검출
- [ ] Flush+Reload 캐시 타이밍 분석 (가능 범위 내)
- [ ] x86-64, ARM (Raspberry Pi / M-series Mac) 이식성 테스트
- [ ] dim-4 (논문) vs dim-8 (본 구현) 속도 비교표 작성
- [ ] 격자 리덕션 품질 통계 (SVP 근사비 분포)

### 체크포인트
- dudect p-value > 0.05 (타이밍 누수 없음)
- 벤치마크 테이블 완성

---

## Phase 5 — 논문 작성 (Month 12)

### 목표
- TCHES 또는 CCS 투고용 논문 초안 완성

### 구성
1. Introduction — 문제 정의, 기여 요약
2. Background — LLL, 상수 시간 요건, FALCON
3. Algorithm — dim-8 브랜치리스 리덕션 설계
4. Security Analysis — 상수 시간 증명 (형식적 or 실험적)
5. Implementation — C++ 구현 세부사항
6. Evaluation — 벤치마크, 사이드채널 실험
7. Conclusion — dim-16 확장 전망

### 투고 타겟
- 1순위: **TCHES** (구현 논문에 우호적, 격자/PQC 트랙)
- 2순위: **CCS** (보안 + 구현 결합 논문)

---

## 기술 스택

| 역할 | 도구 |
|------|------|
| 핵심 구현 | C++ (C++17), 일부 Python (SageMath) |
| 격자 라이브러리 참조 | fplll, NTL |
| 상수 시간 검증 | ctgrind, dudect |
| FALCON 참조 구현 | https://falcon-sign.info |
| 벤치마크 | Google Benchmark, perf |
| 빌드 | CMake |

---

## 참고 자료

- [ePrint 2025/027](https://eprint.iacr.org/2025/027) — 직접 기반 논문 (반드시 완독)
- [ePrint 2025/774](https://eprint.iacr.org/2025/774) — BLASter, 최신 LLL 구현 참고
- [ePrint 2024/1681](https://eprint.iacr.org/2024/1681) — L4 알고리즘 (품질 개선 참고)
- [FALCON spec](https://falcon-sign.info/falcon.pdf) — FALCON 격자 구조
- Nguyen & Stehlé, "Floating-Point LLL Revisited" — 수치 안정성 기초
