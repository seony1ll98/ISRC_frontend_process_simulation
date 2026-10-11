# 검증용 참조 데이터

## nanofab_sithox_2026-10-11.csv

University of Alberta nanoFAB의 [Silicon thermal oxidation calculator](https://toolbox.nanofab.ualberta.ca/sithox/index.php) 출력값입니다. 2026-10-11에 조회했습니다.
이 계산기는 Deal-Grove 모델과 Massoud 모델을 구현한 공개 도구이며, 이 프로젝트와 독립적으로 구현된 결과와 비교하려고 사용했습니다.

| 열 | 의미 |
|---|---|
| ambient | dry(O₂) 또는 wet(H₂O) |
| pp_atm | 산화제 분압 [atm] |
| model | Deal-Grove 또는 Massoud |
| orientation | Si 결정 방향 |
| xi_A | 초기 산화막 두께 [Å] |
| T_C | 산화 온도 [°C] |
| time_min | 산화 시간 [min] |
| thickness_A | 계산기가 출력한 최종 산화막 두께 [Å] (초기 산화막 포함) |

## byu_resistivity_2026-10-11.csv

Brigham Young University Cleanroom의 [Resistivity & Mobility Calculator](https://cleanroom.byu.edu/ResistivityCal) 출력값입니다. 2026-10-11에 웹 페이지의 계산 함수를 직접 실행해 얻었습니다.
이 계산기는 Masetti(1983) 이동도 모델을 사용합니다. 정공(B)의 경우 원 논문의 μ_min2 = 0 대신 μ_min2 = μ_min1로 계산하므로, 이 프로젝트와는 As, P에서만 일치하고 B에서는 차이가 납니다(docs/02_implant_diffusion.md 참고).

| 열 | 의미 |
|---|---|
| dopant | As, B, P |
| N_cm3 | 도핑 농도 [cm⁻³] |
| mobility_cm2_Vs | 이동도 [cm²/Vs] |
| resistivity_ohm_cm | 비저항 [Ω·cm] |
