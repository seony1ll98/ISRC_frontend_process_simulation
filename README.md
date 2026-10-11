# ISRC Front-end Process Simulation

서울대학교 반도체공동연구소(ISRC) 반도체 기본공정교육에서 CMOS 소자 제작 공정을 실습한 뒤, 각 단위 공정을 물리 모델로 다시 계산해 보는 프로젝트입니다.
실습에서 다룬 장비와 공정 순서를 기준으로 모듈을 나누고, 모듈마다 **모델 구현, 교과서 및 독립 구현과의 검증, 한계 분석**을 함께 정리합니다.

## 모듈

| # | 공정 | 모델 | 상태 | 보고서 |
|---|---|---|---|---|
| 01 | 열산화 (Thermal Oxidation) | Deal-Grove, Massoud | 완료 | [docs/01_oxidation.md](docs/01_oxidation.md) |
| 02 | 이온주입 및 확산 (Ion Implantation, Drive-in) | Gaussian 분포, Fick 확산(해석해와 수치해), Masetti 이동도 | 완료 | [docs/02_implant_diffusion.md](docs/02_implant_diffusion.md) |
| 03 | 포토리소그래피 (Photolithography) | Fresnel 회절(근접 노광), Abbe 부분 결맞음 결상, 공정 창 | 완료 | [docs/03_lithography.md](docs/03_lithography.md) |
| 04 | 건식 식각 (Dry Etching) | 식각률, 선택비, 이방성 | 예정 | |
| 05 | 증착 및 금속배선 (CVD, Metallization) | Arrhenius 율속 영역, 면저항 | 예정 | |
| 06 | 공정 통합 (Process Integration) | CMOS 단면 구성, Vth 추정 | 예정 | |

## 01. 열산화 주요 결과

![thickness vs time](figures/oxidation/fig1_thickness_vs_time.png)

- 온도, 시간, 분위기(dry/wet), 결정 방향에 따른 산화막 두께를 계산하고, 목표 두께에 필요한 공정 시간을 역계산합니다.
- 공개 산화 계산기(University of Alberta nanoFAB)와 16개 조건을 비교해 최대 오차 0.063 %를 확인했습니다.
- 10 nm 게이트 산화막을 Deal-Grove 모델로만 설계하면 실제(Massoud 모델 기준) 두께가 약 2배가 됩니다. 얇은 산화막에서는 모델 선택이 공정 결과를 좌우합니다.

## 02. 이온주입 및 확산 주요 결과

![well anneal](figures/implant_diffusion/fig3_well_anneal.png)

- 주입 분포, 열처리 후 분포, 접합 깊이, 면저항을 한 흐름으로 계산합니다(n-well drive-in 설계 예시).
- 확산 해석해와 Crank-Nicolson 수치해가 최대 0.008 % 이내로 일치하고, 이동도 모델은 BYU 계산기와 P, As 14개 조건에서 일치합니다.
- 같은 dose라도 drive-in이 길어질수록 면저항이 낮아집니다(1651 → 755 Ω/sq). 분포가 퍼지며 최대 농도가 낮아지고 이동도가 높아지기 때문입니다.

## 03. 포토리소그래피 주요 결과

![process window](figures/lithography/fig5_process_window.png)

- 근접 노광은 Fresnel 회절로, 축소 투영 노광은 부분 결맞음 결상(Abbe)으로 계산하고, 레지스트 contrast와 비교해 형성 가능한 최소 패턴과 공정 창을 구합니다.
- 근접 노광의 인쇄 폭은 간격과 무관하게 W/√(λg) 하나로 정리되며, CD를 ±10 % 안으로 관리하려면 패턴이 교재의 최소 선폭 √(λg)의 약 2.1배 이상이어야 합니다.
- g-line, i-line, KrF 모두 레지스트 기준 최소 half pitch가 k₁ ≈ 0.40에서 나옵니다. k₁이 장비가 아니라 조명과 레지스트로 정해지는 공정 계수임을 보여줍니다.

## 구조

```
src/isrc_sim/        물리 모델 (constants.py에 계수와 출처)
scripts/             그래프와 결과 표 생성 스크립트
tests/               pytest 검증
data/reference/      검증용 참조 데이터
figures/, results/   생성된 그래프와 표
docs/                모듈별 보고서
```

## 실행

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest
.venv/bin/python scripts/run_oxidation.py
.venv/bin/python scripts/run_implant_diffusion.py
.venv/bin/python scripts/run_lithography.py
```
