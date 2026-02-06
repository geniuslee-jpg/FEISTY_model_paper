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

## 2. NEMO-PISCES 출력 변수 → FEISTY 입력 매칭

### 2-1. 수온 (Tp, Tm, Tb)

| FEISTY | NEMO 변수 | 파일 | 처리 방법 |
|--------|----------|------|----------|
| `Tp` | `votemper` (또는 `thetao`) | grid_T output | 0-100m 깊이의 가중 평균 |
| `Tm` | `votemper` | grid_T output | 500-1500m 깊이의 가중 평균 |
| `Tb` | `votemper` | grid_T output | 해저면(최하층) 값 |

- **단위 변환**: 불필요 (NEMO 출력 = °C, FEISTY 입력 = °C)
- **주의**: zenodo 데이터에는 grid_T 파일이 없음. NEMO를 직접 돌리거나 World Ocean Atlas 2018 사용 필요.

### 2-2. 해저 수심 (depth)

| FEISTY | NEMO 소스 | 파일 |
|--------|----------|------|
| `depth` | `mbathy`, `bottom_level`, 또는 격자 정보 | `domain_cfg.nc` 또는 `mesh_mask.nc` |

- **단위 변환**: 불필요 (m → m)
- **대안**: ETOPO 전지구 수심 데이터 사용 가능

### 2-3. 유광층 깊이 (photic)

| FEISTY | PISCES 변수 | 파일 |
|--------|------------|------|
| `photic` | `heup` (진단 출력 시) | diad_T output |

- **단위 변환**: 불필요 (m → m)
- **현재 zenodo 데이터에 `heup` 없음**. 대안으로 Chl-a에서 계산:

```
PISCES 표층 Chl-a: Csur = NCHL[표층] + DCHL[표층]   (mg Chl-a /m³)

Morel & Berthon (1989) 공식:
  Ctot = 40.6 × Csur^0.459                           (mg Chl-a /m²)
  z_eu = 568.2 × Ctot^(-0.746)                       (m)
```

### 2-4. 소형 동물플랑크톤 생산율 (szprod)

| FEISTY | PISCES 변수 | 파일 |
|--------|------------|------|
| `szprod` | `ZOO` (microzooplankton 농도) + `mzrat` (사망률 파라미터) | ptrc_T.nc + namelist_pisces |

**생물학적 의미**: 소형 동물플랑크톤(0.2-2.0mm)이 상위 포식자에게 잡아먹히는 손실량 = FEISTY에서의 최대 생산율(r·Rmax)

**PISCES에서 이 손실량은 quadratic mortality (closure term)로 표현됨:**

```
closure_ZOO(x,y,z) = mzrat × ZOO(x,y,z)²
```

- `mzrat = 0.02 d⁻¹ (mmol C/m³)⁻¹` (namelist_pisces_ref 값)
- `ZOO` = microzooplankton 농도 (mol C/L in PISCES output)

**단위 변환 과정:**

```
Step 1: PISCES 농도 단위 통일
  ZOO [mol C/L] × 1000 = ZOO [mmol C/m³]

Step 2: Closure term 계산 (각 격자점, 각 깊이)
  closure = mzrat × ZOO²
  = 0.02 [d⁻¹·(mmol C/m³)⁻¹] × ZOO² [(mmol C/m³)²]
  = [mmol C /m³ /d]

Step 3: 수직 적분 (연직 전체 합산)
  closure_2D = Σ(closure × dz)    [mmol C /m² /d]

Step 4: mmol C → g C
  × 12/1000 = [g C /m² /d]         (C 원자량 = 12 g/mol)

Step 5: 일 → 년
  × 365 = [g C /m² /yr]

Step 6: g C → g WW (wet weight)
  × 9 = [g WW /m² /yr]             (WW:C = 9:1)
```

**요약 공식:**
```
szprod = Σ_z(mzrat × ZOO² × dz) × 12/1000 × 365 × 9
       = Σ_z(0.02 × ZOO² × dz) × 39.42
```
(여기서 39.42 = 12/1000 × 365 × 9)

### 2-5. 대형 동물플랑크톤 생산율 (lzprod)

| FEISTY | PISCES 변수 | 파일 |
|--------|------------|------|
| `lzprod` | `ZOO2` (mesozooplankton 농도) + `mzrat2` (사망률 파라미터) | ptrc_T.nc + namelist_pisces |

**계산 방식은 szprod와 동일하되, 온도 보정이 추가됨:**

```
Step 1-6: szprod와 동일한 과정으로 closure term 계산
  raw_lzprod = Σ_z(mzrat2 × ZOO2² × dz) × 12/1000 × 365 × 9
             = Σ_z(0.01 × ZOO2² × dz) × 39.42
```

- `mzrat2 = 0.01 d⁻¹ (mmol C/m³)⁻¹` (namelist_pisces_ref 값)

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

### 2-6. 해저면 Detrital Flux (dfbot)

| FEISTY | PISCES 변수 | 파일 |
|--------|------------|------|
| `dfbot` | `EPC100` (100m 깊이의 export flux) | diad_T.nc |

- `EPC100` = Export Production of Carbon at 100m (mol C /m² /s)

**직접적인 해저면 flux가 아니므로 감쇄 추정 필요:**

```
Step 1: 100m flux에서 해저면 flux로 감쇄 (Martin curve)
  dfbot_raw = EPC100 × (depth / 100)^(-0.858)    [mol C /m² /s]

Step 2: mol C → g C
  × 12 = [g C /m² /s]

Step 3: 초 → 년
  × 86400 × 365 = [g C /m² /yr]

Step 4: g C → g WW
  × 9 = [g WW /m² /yr]
```

**요약 공식:**
```
dfbot = EPC100 × (depth/100)^(-0.858) × 12 × 86400 × 365 × 9
      = EPC100 × (depth/100)^(-0.858) × 3.40e9
```
(여기서 3.40e9 = 12 × 86400 × 365 × 9)

> **주의**: Martin curve 지수(-0.858)는 일반적인 값. 원래 FEISTY 논문은 COBALT에서 직접 해저면 flux를 가져왔으므로 이 감쇄 추정은 근사값임.

---

## 3. 전체 요약 표

| # | FEISTY 변수 | FEISTY 단위 | NEMO-PISCES 소스 | PISCES 원래 단위 | 변환 공식 |
|---|------------|------------|-----------------|-----------------|----------|
| 1 | `Tp` | °C | `votemper` (0-100m avg) | °C | **변환 불필요** |
| 2 | `Tm` | °C | `votemper` (500-1500m avg) | °C | **변환 불필요** |
| 3 | `Tb` | °C | `votemper` (bottom) | °C | **변환 불필요** |
| 4 | `depth` | m | `domain_cfg` / bathymetry | m | **변환 불필요** |
| 5 | `photic` | m | `NCHL + DCHL` (표층) | mg Chl-a /m³ | Morel & Berthon 공식 |
| 6 | `szprod` | g WW/m²/yr | `ZOO` + `mzrat` | mol C/L | `Σ(mzrat×ZOO²×dz) × 39.42` |
| 7 | `lzprod` | g WW/m²/yr | `ZOO2` + `mzrat2` | mol C/L | `Σ(mzrat2×ZOO2²×dz) × 39.42 / ratio(Tp)` |
| 8 | `dfbot` | g WW/m²/yr | `EPC100` | mol C/m²/s | `EPC100 × (depth/100)^(-0.858) × 3.40e9` |

---

## 4. 현재 zenodo 데이터 가용성

| 변수 | zenodo 데이터로 가능? | 비고 |
|------|---------------------|------|
| `Tp`, `Tm`, `Tb` | **불가** | grid_T 물리 출력 없음 → WOA 2018 또는 NEMO 재실행 필요 |
| `depth` | **불가** | domain_cfg 없음 → ETOPO 또는 ORCA2 격자 파일 필요 |
| `photic` | **가능** | ptrc_T.nc의 `NCHL + DCHL` 표층값으로 계산 |
| `szprod` | **가능** | ptrc_T.nc의 `ZOO` + namelist `mzrat=0.02` |
| `lzprod` | **가능** | ptrc_T.nc의 `ZOO2` + namelist `mzrat2=0.01` (+ Tp 보정 필요) |
| `dfbot` | **가능** | diad_T.nc의 `EPC100` + depth 필요 (Martin curve 감쇄) |

---

## 5. 참고 문헌

- FEISTY input sources: `quarto/data/Input_global_sources.txt`
- PISCES-v2: Aumont et al. (2015), *Geosci. Model Dev.*, 8, 2465-2513
- Martin curve: Martin et al. (1987), *Deep-Sea Research*, 34(2), 267-285
- Morel & Berthon (1989), *Limnology and Oceanography*, 34(8), 1545-1562
- WOA 2018: Locarnini et al. (2018), NOAA Atlas NESDIS 81
