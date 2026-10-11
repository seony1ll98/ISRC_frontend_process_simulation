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
