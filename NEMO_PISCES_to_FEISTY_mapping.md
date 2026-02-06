# NEMO-PISCES → FEISTY 입력 변수 매칭 가이드

## 1. FEISTY 모델이 필요한 입력 변수 (8개)

FEISTY 글로벌 런(`Global_run_illustration.R`)은 `Input_global.csv` 파일에서 아래 8개 변수를 읽어
`setupVertical2()` 함수에 전달합니다.

| # | CSV 컬럼 | 의미 | FEISTY 최종 단위 |
|---|----------|------|-----------------|
| 1 | `Tp` | 표층 수온 (0-100m 평균) | **°C** |
| 2 | `Tm` | 중층 수온 (500-1500m 평균) | **°C** |
| 3 | `Tb` | 저층 수온 (해저면) | **°C** |
| 4 | `depth` | 해저 수심 | **m** |
| 5 | `photic` | 유광층(euphotic zone) 깊이 | **m** |
| 6 | `szprod` | 소형 동물플랑크톤 최대 생산율 (r·Rmax) | **g WW m⁻² yr⁻¹** |
| 7 | `lzprod` | 대형 동물플랑크톤 최대 생산율 (r·Rmax) | **g WW m⁻² yr⁻¹** |
| 8 | `dfbot` | 해저면 도달 detrital flux | **g WW m⁻² yr⁻¹** |

> **WW = Wet Weight**. 원래 FEISTY 논문에서는 COBALT 출력(g C)을 wet weight:carbon = **9:1** 비율로 변환하였음.

---

## 2. NEMO-PISCES 출력 설정 현황 (output_sy/ XML 기반)

### 2-0. 출력 파일별 사용 가능한 변수

| 출력 파일 | 시간 해상도 | FEISTY에 필요한 변수 |
|----------|-----------|---------------------|
| **grid_T** | 5일 평균 | `thetao` (3D 수온), **`sbt`** (해저면 수온) |
| **ptrc_T** | 월평균 | `NCHL`, `DCHL` (엽록소) |
| **ptrc_T** | 연평균 | **`ZOO`**, **`ZOO2`** (동물플랑크톤 농도) |
| **diad_T** | 연평균 | **`Heup`** (유광층 깊이), **`EPC100`** (100m export flux), `GRAZ1`, `GRAZ2` |

### 주의: GRAZ1/GRAZ2 vs Closure Term

- `GRAZ1` = 미소동물플랑크톤이 **먹이를 섭식하는 양** (grazing BY zoo)
- `GRAZ2` = 중형동물플랑크톤이 **먹이를 섭식하는 양** (grazing BY zoo)
- FEISTY가 필요한 것 = 동물플랑크톤이 **상위 포식자에게 잡아먹히는 양** (loss OF zoo)
- → **GRAZ1/GRAZ2는 직접 사용 불가**, closure term(quadratic mortality)으로 계산 필요

---

## 3. 변수별 매칭 및 단위 변환 상세

### 3-1. Tp — 표층 수온 (0-100m 평균)

| 항목 | 내용 |
|------|------|
| NEMO 변수 | `thetao` (grid_T 파일, 5일 평균) |
| PISCES 단위 | **°C** |
| FEISTY 단위 | **°C** |
| 단위 변환 | **불필요** |
| 처리 방법 | 0-100m 깊이 레벨에 대해 e3t(셀 두께) 가중 평균 |

```python
# 개념 코드
Tp = Σ(thetao[z] × e3t[z]) / Σ(e3t[z])   for z where depth ≤ 100m
```

### 3-2. Tm — 중층 수온 (500-1500m 평균)

| 항목 | 내용 |
|------|------|
| NEMO 변수 | `thetao` (grid_T 파일, 5일 평균) |
| PISCES 단위 | **°C** |
| FEISTY 단위 | **°C** |
| 단위 변환 | **불필요** |
| 처리 방법 | 500-1500m 깊이 레벨에 대해 e3t 가중 평균 |

```python
# 개념 코드
Tm = Σ(thetao[z] × e3t[z]) / Σ(e3t[z])   for z where 500m ≤ depth ≤ 1500m
```

### 3-3. Tb — 저층 수온 (해저면)

| 항목 | 내용 |
|------|------|
| NEMO 변수 | **`sbt`** (grid_T 파일, 5일 평균) — 해저면 수온 직접 출력! |
| PISCES 단위 | **°C** |
| FEISTY 단위 | **°C** |
| 단위 변환 | **불필요** |
| 처리 방법 | **직접 사용** (계산 불필요) |

### 3-4. depth — 해저 수심

| 항목 | 내용 |
|------|------|
| NEMO 소스 | `domain_cfg.nc` 또는 `mesh_mask.nc`의 bathymetry |
| 단위 | **m** |
| 단위 변환 | **불필요** |
| 대안 | e3t 수직 합산으로 추정: `depth ≈ Σ(e3t[z])` (해양 격자점에서) |

> `domain_cfg.nc` 파일 위치를 확인해야 함

### 3-5. photic — 유광층 깊이

| 항목 | 내용 |
|------|------|
| PISCES 변수 | **`Heup`** (yearly diad_T 파일) |
| PISCES 단위 | **m** |
| FEISTY 단위 | **m** |
| 단위 변환 | **불필요** |
| 처리 방법 | **직접 사용** (계산 불필요) |

### 3-6. szprod — 소형 동물플랑크톤 최대 생산율

| 항목 | 내용 |
|------|------|
| PISCES 변수 | `ZOO` (yearly ptrc_T) + namelist 파라미터 `mzrat` |
| ZOO 단위 | **mmol C /m³** (field_def_nemo-pisces.xml에서 확인) |
| FEISTY 단위 | **g WW /m² /yr** |

**생물학적 의미**: 소형 동물플랑크톤이 상위 포식자에게 잡아먹히는 손실량(closure term)
= FEISTY에서의 최대 생산율(r·Rmax)

**PISCES closure term (quadratic mortality):**
```
closure_ZOO(x,y,z) = mzrat × ZOO(x,y,z)²
```
- `mzrat = 0.02 d⁻¹ (mmol C/m³)⁻¹` (namelist_pisces_ref)

**단위 변환 과정:**

```
Step 1: Closure term 계산 (각 격자점, 각 깊이)
  ZOO 단위가 이미 mmol C/m³ (PISCES 기본 단위)
  closure = mzrat × ZOO²
  = 0.02 [d⁻¹·(mmol C/m³)⁻¹] × ZOO² [(mmol C/m³)²]
  = [mmol C /m³ /d]

Step 2: 수직 적분 (연직 전체 합산, e3t = 셀 두께)
  closure_2D = Σ(closure × e3t)    [mmol C /m² /d]

Step 3: mmol C → g C
  × 12/1000 = [g C /m² /d]         (C 원자량 = 12 g/mol)

Step 4: 일 → 년
  × 365 = [g C /m² /yr]

Step 5: g C → g WW (wet weight)
  × 9 = [g WW /m² /yr]             (WW:C = 9:1)
```

**요약 공식:**
```
szprod = Σ_z(0.02 × ZOO² × e3t) × 12/1000 × 365 × 9
       = Σ_z(0.02 × ZOO² × e3t) × 39.42
```

### 3-7. lzprod — 대형 동물플랑크톤 최대 생산율

| 항목 | 내용 |
|------|------|
| PISCES 변수 | `ZOO2` (yearly ptrc_T) + namelist 파라미터 `mzrat2` |
| ZOO2 단위 | **mmol C /m³** |
| FEISTY 단위 | **g WW /m² /yr** |

**계산 방식은 szprod와 동일 + 온도 보정:**

```
Step 1-5: szprod와 동일
  raw_lzprod = Σ_z(0.01 × ZOO2² × e3t) × 39.42
```
- `mzrat2 = 0.01 d⁻¹ (mmol C/m³)⁻¹` (namelist_pisces_ref)

**추가: 온도 보정 (FEISTY 원논문 방법)**

대형 동물플랑크톤은 어류에 의해 top-down 제어되므로, closure term ≠ 최대 생산량.
온도별 보정 비율을 적용하여 최대 생산량을 추정:

```
lzprod = raw_lzprod / ratio(Tp)
```

| 표층 수온 (Tp) | ratio | 의미 |
|---------------|-------|------|
| 30°C | 0.91 | closure term은 최대 생산의 91% |
| 20°C | 0.89 | closure term은 최대 생산의 89% |
| 10°C | 0.85 | closure term은 최대 생산의 85% |
| 0°C  | 0.78 | closure term은 최대 생산의 78% |

중간 온도는 선형 보간(linear interpolation) 적용.

### 3-8. dfbot — 해저면 Detrital Flux

| 항목 | 내용 |
|------|------|
| PISCES 변수 | `EPC100` (yearly diad_T) |
| EPC100 단위 | **mol C /m² /s** (field_def에서 확인) |
| FEISTY 단위 | **g WW /m² /yr** |

**100m flux → 해저면 flux 감쇄 추정 (Martin curve):**

```
Step 1: Martin curve 감쇄
  dfbot_raw = EPC100 × (depth / 100)^(-0.858)    [mol C /m² /s]

Step 2: mol C → g C
  × 12                                            [g C /m² /s]

Step 3: 초 → 년
  × 86400 × 365                                   [g C /m² /yr]

Step 4: g C → g WW
  × 9                                             [g WW /m² /yr]
```

**요약 공식:**
```
dfbot = EPC100 × (depth/100)^(-0.858) × 12 × 86400 × 365 × 9
      = EPC100 × (depth/100)^(-0.858) × 3.407×10⁹
```

> **주의**: Martin curve 지수(-0.858)는 일반적 값. FEISTY 원논문은 COBALT에서 직접 해저면 flux를 가져왔으므로 이 감쇄 추정은 근사값.

---

## 4. 전체 요약 표

| # | FEISTY 변수 | FEISTY 단위 | NEMO-PISCES 변수 | PISCES 단위 | 변환 공식 | 난이도 |
|---|------------|------------|-----------------|------------|----------|--------|
| 1 | `Tp` | °C | `thetao` (0-100m 가중평균) | °C | **변환 불필요** | 쉬움 |
| 2 | `Tm` | °C | `thetao` (500-1500m 가중평균) | °C | **변환 불필요** | 쉬움 |
| 3 | `Tb` | °C | **`sbt`** (직접 사용) | °C | **변환 불필요** | 쉬움 |
| 4 | `depth` | m | `domain_cfg` / `e3t` 합산 | m | **변환 불필요** | 쉬움 |
| 5 | `photic` | m | **`Heup`** (직접 사용) | m | **변환 불필요** | 쉬움 |
| 6 | `szprod` | g WW/m²/yr | `ZOO` + `mzrat=0.02` | mmol C/m³ | `Σ(0.02×ZOO²×e3t) × 39.42` | 중간 |
| 7 | `lzprod` | g WW/m²/yr | `ZOO2` + `mzrat2=0.01` | mmol C/m³ | `Σ(0.01×ZOO2²×e3t) × 39.42 / ratio(Tp)` | 중간 |
| 8 | `dfbot` | g WW/m²/yr | `EPC100` + `depth` | mol C/m²/s | `EPC100 × (depth/100)^(-0.858) × 3.407e9` | 중간 |

> 단위 변환 상수: `39.42 = 12/1000 × 365 × 9` , `3.407e9 = 12 × 86400 × 365 × 9`

---

## 5. 아직 확인 필요한 사항

| 항목 | 상태 | 필요한 조치 |
|------|------|-----------|
| `domain_cfg.nc` 파일 위치 | **미확인** | 서버에서 find 명령어로 확인 필요 |
| grid_T 실제 출력 nc 파일 | **미확인** | 본인 NEMO 실행 결과 파일 위치 확인 필요 |
| `mzrat`, `mzrat2` 실제 사용값 | namelist_ref 값 확인 완료 | namelist_cfg에서 오버라이드 여부 확인 필요 |

---

## 6. 참고 문헌

- FEISTY input sources: `quarto/data/Input_global_sources.txt`
- PISCES-v2: Aumont et al. (2015), *Geosci. Model Dev.*, 8, 2465-2513
- Martin curve: Martin et al. (1987), *Deep-Sea Research*, 34(2), 267-285
- Morel & Berthon (1989), *Limnology and Oceanography*, 34(8), 1545-1562
- WOA 2018: Locarnini et al. (2018), NOAA Atlas NESDIS 81
- NEMO output XML: `/data01/labdisk/sungjin/NEMO_output/output_sy/file_def_nemo-oce.xml`
- PISCES output XML: `/data01/labdisk/sungjin/NEMO_output/output_sy/file_def_nemo-pisces.xml`
- PISCES field definitions: `/data01/labdisk/sungjin/NEMO_output/output_sy/field_def_nemo-pisces.xml`
