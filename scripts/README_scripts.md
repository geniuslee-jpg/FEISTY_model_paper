# NEMO-FEISTY 데이터 처리 스크립트

## 워크플로우

### Step 1: 데이터 탐색
```bash
cd /data01/labdisk/sungjin/NEMO-FEISTY
python /path/to/scripts/01_inspect_nemo_output.py .
```
- 각 NetCDF 파일의 변수명, 차원, 단위 확인
- 출력 결과를 저장해둘 것: `python 01_inspect_nemo_output.py > inspect_result.txt`

### Step 2: 변수명 매핑 수정
`02_extract_feisty_input.py` 상단 CONFIG 섹션에서 실제 변수명으로 수정:
- `VAR_TEMP`: 온도 변수명 (예: `votemper`, `thetao`)
- `VAR_SZPROD`: 소형 동물플랑크톤 생산 변수명
- `VAR_LZPROD`: 대형 동물플랑크톤 생산 변수명
- `VAR_DFBOT`: 해저 detrital flux 변수명
- `VAR_PHOTIC`: 유광층 깊이 변수명

### Step 3: 데이터 추출
```bash
python /path/to/scripts/02_extract_feisty_input.py
```
- 10년 평균 계산
- `quarto/data/Input_NEMO_FEISTY.csv` 생성
- 검증 플롯 `Input_NEMO_FEISTY_check.png` 생성

### Step 4: FEISTY 실행
기존 `Global_run_illustration.R`에서 입력 파일만 교체:
```r
glob <- read.csv(file = "data/Input_NEMO_FEISTY.csv")
```

## 필요 패키지 (Python)
```
xarray
netCDF4
pandas
numpy
matplotlib (선택, 검증 플롯용)
```

## FEISTY Input 형식
| 변수 | 설명 | 단위 |
|------|------|------|
| lon | 경도 | degrees |
| lat | 위도 | degrees |
| lzprod | 대형 동물플랑크톤 생산 | g ww m⁻² yr⁻¹ |
| szprod | 소형 동물플랑크톤 생산 | g ww m⁻² yr⁻¹ |
| dfbot | 해저 detrital flux | g m⁻² yr⁻¹ |
| photic | 유광층 깊이 | m |
| depth | 해저 깊이 | m |
| Tb | 저층 온도 | °C |
| Tm | 중층 온도 (500-1500m) | °C |
| Tp | 표층 온도 (0-100m) | °C |
