# FAR23_LOADS

## Ch3 Weight Envelope 사용법

'그 긴거 시작'

### 1. 사전 준비

Ch3는 Ch5의 wing geometry 결과에서 `MAC`, `XMACLE` 값을 가져와 structural limit을 계산한다.

따라서 Ch3를 실행하기 전에 Ch5를 먼저 실행해야 한다.

---
outputs/cases/ch5_geometry_output.json
---

필요


### 2. Ch3 입력 DB 확인

Ch3 입력 DB는 다음 파일이다.

---
inputs/weight_db.json
---

이 파일에는 Ch3에서 사용할 weight component DB가 들어간다.

구성은 다음 기준으로 나뉜다.

---
C = 0 : Empty weight item
C = 1 : Minimum weight item
C = 2 : Discretionary load item
---

주의할 점:

- `BALLAST`는 `weight_db.json`에 직접 넣지 않는다.
- `BALLAST`는 Ch3 ballast case 생성 코드가 자동으로 계산해서 Ch4용 DB에 추가한다.
- `MAC`, `XMACLE`은 `weight_db.json`에 직접 쓰지 않고 Ch5 output에서 읽는다.


### 3. Ch3 Weight Envelope 계산

다음 코드를 실행하여 Ch3 envelope를 계산한다.

---
chapters/ch03_weight_envelope.py
---

생성되는 주요 파일:

---
outputs/cases/ch3_weight_envelope_output.json
outputs/reports/ch3_weight_envelope_report.txt
---


`load_points`에는 그래프에서 선택할 수 있는 point 번호가 들어간다. 이 번호는 ballast case를 만들 때 사용한다.


### 4. Ch3 Envelope PNG 생성

Ch3 output을 바탕으로 PNG 그래프를 생성한다.

---
visualization/weight_envelope_plot.py
---

코드 실행

생성되는 파일:

---
outputs/figures/ch3_weight_envelope.png
---

그래프에는 각 load point 번호가 표시된다.

---
예시:
[7] 5TH PERSON
[8] BAGGAGE / FUEL TO FULL
---

같은 최종 full-load 점이 forward edge와 aft edge에 동시에 존재하는 경우, 하나의 point 번호로 합쳐서 표시한다.


### 5. Ballast Case 설정

그래프에서 원하는 load point를 고른 뒤, 다음 파일에서 point 번호를 설정한다.

---
inputs/ballast_case_config.json
---

주요 설정값:

---
target             : 맞추고 싶은 structural limit point
selected_point_id  : 그래프에서 선택한 load point 번호
ballast_z          : ballast의 Z 위치
output_input_path  : 생성될 Ch4용 입력 DB
output_report_path : ballast 계산 요약 report
---

현재 사용 가능한 target 예시는 다음과 같다.

---
FWDGROSS
AFTGROSS
FWDRED
---

### 6. Ballast Case 생성

설정한 point 번호를 기준으로 ballast를 계산하고, Ch4에서 사용할 수 있는 DB를 생성한다.

---
chapters/ch03_make_ballast_case.py
---

코드 실행

생성되는 파일:

---
inputs/weight_db_for_ch4.json
outputs/reports/ch3_aft_gross_ballast_case_report.txt
---

`weight_db_for_ch4.json`에는 선택한 load 구성과 자동 계산된 `BALLAST`가 포함된다.

예를 들어 Aft Gross 조건에서는 다음 개념으로 ballast가 계산된다.

---
목표 조건:
WL = target weight
XL = target XCG

선택한 load point:
WA = selected point weight
XA = selected point XCG

Ballast weight:
WB = WL - WA

Ballast X location:
XB = (WL * XL - WA * XA) / WB
---

즉 ballast는 임의 위치에 추가되는 것이 아니라, 최종 총중량과 CG가 target point에 맞도록 자동 계산된다.


### 7. Ch4와 연결

Ch3에서 생성한 다음 파일을 Ch4 입력 DB로 사용한다.

---
inputs/weight_db_for_ch4.json
---

이 파일은 선택한 load point에 ballast가 추가된 case DB이다. Ch4는 이 DB를 일반 component 목록으로 읽어서 weight, CG, inertia를 계산하면 된다.
