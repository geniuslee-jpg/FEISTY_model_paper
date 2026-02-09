"""
02_extract_feisty_input.py
==========================
NEMO-FEISTY NetCDF 출력에서 FEISTY R 패키지 입력 형식(Input_global.csv)으로 변환.

10년 시간 평균을 계산하여 climatology 형태의 입력 데이터를 생성.

사용법:
    python 02_extract_feisty_input.py

출력:
    ../quarto/data/Input_NEMO_FEISTY.csv

※ 01_inspect_nemo_output.py 결과를 보고 아래 CONFIG 섹션의 변수명을 확인/수정할 것.
"""

import os
import sys
import numpy as np
import xarray as xr
import pandas as pd

# ============================================================================
# CONFIG - 이 부분을 01_inspect 결과에 맞게 수정하세요
# ============================================================================

DATA_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY"

# 시간 범위 설정 (연도 기준, model year)
# 50년 데이터 중 마지막 10년 사용 (예: year 3-12 또는 1-10)
# 현재 12년 돌리고 있으므로 year 1-10 사용
YEAR_START = 1   # 시작 연도 (model year)
YEAR_END = 10    # 종료 연도 (model year)

# 출력 파일 경로
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "quarto", "data")
OUTPUT_FILE = "Input_NEMO_FEISTY.csv"

# --- 파일명 패턴 ---
# 월평균 데이터 (biogeochemistry)
FILE_PTRC_1M = "ORCA2_1m_00010101_00501231_ptrc_T.nc"     # passive tracers
FILE_DIAD_1M = "ORCA2_1m_00010101_00501231_diad_T.nc"     # diagnostics
# 연평균 데이터
FILE_PTRC_1Y = "ORCA2_1y_00010101_00501231_ptrc_T.nc"
FILE_DIAD_1Y = "ORCA2_1y_00010101_00501231_diad_T.nc"
# Grid 데이터
FILE_GRID_T = "ORCA2_5d_00010101_00501231_grid_T.nc"       # temperature
FILE_STATIC = "ORCA2_1y_00010101_00501231_grid_T_static.nc" # static fields

# --- 변수명 (01_inspect 결과 보고 수정) ---
# NEMO-PISCES 표준 변수명을 기본으로 설정, 실제 파일 확인 후 수정 필요

# 온도 변수 (grid_T 파일)
VAR_TEMP = "votemper"       # 또는 "thetao", "toce" 등

# 좌표/격자 변수 (grid_T 또는 static 파일)
VAR_LON = "nav_lon"         # 경도
VAR_LAT = "nav_lat"         # 위도
VAR_DEPTH_AXIS = "deptht"   # 깊이 축 (좌표) - 또는 "depth", "gdept_1d"
VAR_BATHY = "mbathy"        # 해저 깊이 - 또는 "bathy_level", "bottom_level"
# static 파일에서 실제 수심(m): "e3t_0"의 합 또는 별도 변수

# PISCES 생지화학 변수 (ptrc_T 또는 diad_T 파일)
# 동물플랑크톤 관련 - diad_T에 production 관련 변수가 있을 수 있음
# 일반적인 PISCES 변수:
VAR_ZOO = "ZOO"             # microzooplankton biomass (mmol C/m3)
VAR_ZOO2 = "ZOO2"           # mesozooplankton biomass (mmol C/m3)

# 동물플랑크톤 생산 관련 (diad_T)
# PISCES 진단 변수 후보: GRAZ1 (grazing by ZOO), GRAZ2 (grazing by ZOO2),
# ZPROD (ZOO production), ZPROD2 (ZOO2 production) 등
VAR_SZPROD = None            # 소형 zoo 생산 - inspect 후 설정
VAR_LZPROD = None            # 대형 zoo 생산 - inspect 후 설정

# Detrital flux 관련 (diad_T)
# PISCES 진단 변수 후보: EPC100 (export at 100m), EPCAL100 (CaCO3 export),
# Cflx (carbon flux), POC_FLUX_BOTTOM 등
VAR_DFBOT = None             # 해저 detrital flux - inspect 후 설정

# 유광층 깊이 (diad_T)
# PISCES 진단 변수 후보: ZEU, hmld, Heup 등
VAR_PHOTIC = None            # 유광층 깊이 - inspect 후 설정

# --- 단위 변환 상수 ---
# PISCES는 보통 mmol C/m3 단위, FEISTY는 g wet weight/m2/yr
# C_to_WW: mmol C -> g wet weight (12 mg/mmol * 9 ww/C ratio / 1000 mg/g = 0.108)
MMOLC_TO_GWW = 12e-3 * 9  # = 0.108 g ww per mmol C

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def select_time_range(ds, year_start, year_end):
    """시간 차원에서 특정 연도 범위 선택.

    NEMO 모델의 시간 좌표가 다양할 수 있으므로 여러 방법 시도.
    """
    time_dim = None
    for tname in ['time_counter', 'time', 't']:
        if tname in ds.dims:
            time_dim = tname
            break

    if time_dim is None:
        print("WARNING: 시간 차원을 찾을 수 없습니다. 전체 데이터 사용.")
        return ds

    time_var = ds[time_dim]
    n_time = len(time_var)

    # 방법 1: 시간 좌표가 datetime인 경우
    try:
        years = pd.DatetimeIndex(time_var.values).year
        mask = (years >= year_start) & (years <= year_end)
        if mask.any():
            print(f"  시간 선택 (datetime): year {year_start}-{year_end}, "
                  f"{mask.sum()}/{n_time} timesteps")
            return ds.isel({time_dim: mask})
    except Exception:
        pass

    # 방법 2: cftime 객체인 경우
    try:
        import cftime
        times = time_var.values
        if hasattr(times[0], 'year'):
            years = np.array([t.year for t in times])
            mask = (years >= year_start) & (years <= year_end)
            if mask.any():
                print(f"  시간 선택 (cftime): year {year_start}-{year_end}, "
                      f"{mask.sum()}/{n_time} timesteps")
                return ds.isel({time_dim: mask})
    except Exception:
        pass

    # 방법 3: 파일이 월별이라면 인덱스로 선택
    # year 1 = index 0-11 (month), year_start-year_end
    idx_start = (year_start - 1) * 12
    idx_end = year_end * 12
    if idx_end > n_time:
        idx_end = n_time
        print(f"  WARNING: 요청 범위가 데이터를 초과. {n_time} timesteps까지만 사용.")

    print(f"  시간 선택 (index): timestep {idx_start}-{idx_end-1} "
          f"(year {year_start}-{year_end}), {idx_end-idx_start}/{n_time} timesteps")
    return ds.isel({time_dim: slice(idx_start, idx_end)})


def select_time_range_yearly(ds, year_start, year_end):
    """연평균 파일에서 시간 범위 선택."""
    time_dim = None
    for tname in ['time_counter', 'time', 't']:
        if tname in ds.dims:
            time_dim = tname
            break

    if time_dim is None:
        return ds

    n_time = len(ds[time_dim])

    # 연평균이라면 year 1 = index 0
    idx_start = year_start - 1
    idx_end = year_end
    if idx_end > n_time:
        idx_end = n_time

    print(f"  시간 선택 (yearly): index {idx_start}-{idx_end-1}, "
          f"{idx_end-idx_start}/{n_time} years")
    return ds.isel({time_dim: slice(idx_start, idx_end)})


def depth_average(da, depth_coord, z_min, z_max):
    """특정 깊이 범위의 가중 평균 계산.

    Parameters
    ----------
    da : xr.DataArray
        3D+ 데이터 (time, depth, y, x)
    depth_coord : str
        깊이 좌표 이름
    z_min, z_max : float
        깊이 범위 (m, 양수)
    """
    depths = da[depth_coord].values
    mask = (depths >= z_min) & (depths <= z_max)

    if not mask.any():
        print(f"  WARNING: 깊이 {z_min}-{z_max}m 범위에 데이터 없음.")
        # 가장 가까운 깊이 사용
        idx = np.argmin(np.abs(depths - (z_min + z_max) / 2))
        return da.isel({depth_coord: idx})

    da_sub = da.isel({depth_coord: mask})

    # 깊이 축 기준 평균
    return da_sub.mean(dim=depth_coord)


def get_bottom_values(da, depth_coord, mask_var=None):
    """각 격자점에서 최저층(해저) 값 추출.

    Parameters
    ----------
    da : xr.DataArray
        3D+ 데이터 (time, depth, y, x)
    depth_coord : str
        깊이 좌표 이름
    """
    # NaN이 아닌 마지막 깊이 인덱스 찾기
    # 시간 평균 후 처리
    if 'time_counter' in da.dims:
        da_mean = da.mean(dim='time_counter')
    elif 'time' in da.dims:
        da_mean = da.mean(dim='time')
    else:
        da_mean = da

    # 각 (y,x) 격자에서 유효한 마지막 깊이층의 값
    valid = da_mean.notnull()
    # argmax on reversed depth gives last valid index
    n_depth = len(da_mean[depth_coord])

    # 방법: 뒤집어서 첫 번째 유효값 찾기
    valid_rev = valid.isel({depth_coord: slice(None, None, -1)})
    last_valid_rev = valid_rev.argmax(dim=depth_coord)
    last_valid_idx = n_depth - 1 - last_valid_rev

    # advanced indexing으로 bottom 값 추출
    result = da_mean.isel({depth_coord: last_valid_idx})
    return result


def compute_depth_from_grid(ds_static, depth_coord):
    """Static 파일에서 해저 깊이(m) 계산.

    여러 방법 시도:
    1. 'bathy_level' 또는 직접적인 bathymetry 변수
    2. 깊이축의 마지막 유효 레벨에서 깊이 추출
    """
    # 직접적인 수심 변수 확인
    for vname in ['bathy_meter', 'Bathymetry', 'bathymetry', 'ht_0', 'ht']:
        if vname in ds_static:
            print(f"  해저 깊이 변수 발견: {vname}")
            return ds_static[vname].squeeze()

    # mbathy (깊이 레벨 인덱스) + 깊이 좌표로 계산
    for vname in ['mbathy', 'bottom_level']:
        if vname in ds_static and depth_coord in ds_static:
            depths = ds_static[depth_coord].values
            bathy_idx = ds_static[vname].squeeze().values.astype(int)
            bathy_m = np.where(bathy_idx > 0,
                               depths[np.clip(bathy_idx - 1, 0, len(depths) - 1)],
                               np.nan)
            return xr.DataArray(bathy_m,
                                dims=ds_static[vname].squeeze().dims,
                                coords=ds_static[vname].squeeze().coords)

    print("  WARNING: 해저 깊이 변수를 찾을 수 없습니다.")
    return None


def compute_photic_depth_from_chl(chl_surf):
    """표층 Chl-a (mg/m3)로부터 유광층 깊이 계산.

    Morel & Berthon (1989):
        C_tot = 40.6 * C_sur^0.459
        Z_eu = 568.2 * C_tot^(-0.746)
    """
    c_tot = 40.6 * np.power(np.maximum(chl_surf, 1e-6), 0.459)
    z_eu = 568.2 * np.power(c_tot, -0.746)
    return z_eu


# ============================================================================
# MAIN EXTRACTION
# ============================================================================

def main():
    print("=" * 70)
    print("NEMO-FEISTY 출력 → FEISTY 입력 변환")
    print(f"데이터 경로: {DATA_DIR}")
    print(f"시간 범위: Year {YEAR_START} - {YEAR_END}")
    print("=" * 70)

    # ------------------------------------------------------------------
    # 1. 좌표 및 격자 정보 로드 (static 파일)
    # ------------------------------------------------------------------
    print("\n[1] 격자/좌표 정보 로드...")
    static_path = os.path.join(DATA_DIR, FILE_STATIC)
    ds_static = xr.open_dataset(static_path, decode_times=False)
    print(f"  Static 파일 변수: {list(ds_static.data_vars)}")

    # 경위도
    if VAR_LON in ds_static:
        lon2d = ds_static[VAR_LON].squeeze().values
        lat2d = ds_static[VAR_LAT].squeeze().values
    else:
        # 좌표에서 찾기
        lon2d = ds_static.coords[VAR_LON].values
        lat2d = ds_static.coords[VAR_LAT].values

    print(f"  격자 형태: {lon2d.shape}")

    # 해저 깊이
    depth_bottom = compute_depth_from_grid(ds_static, VAR_DEPTH_AXIS)

    ds_static.close()

    # ------------------------------------------------------------------
    # 2. 온도 처리 (grid_T 파일) → Tp, Tm, Tb
    # ------------------------------------------------------------------
    print("\n[2] 온도 데이터 처리...")
    grid_path = os.path.join(DATA_DIR, FILE_GRID_T)

    # 5d 파일은 클 수 있으므로 lazy loading
    ds_grid = xr.open_dataset(grid_path, chunks={'time_counter': 'auto'},
                              decode_times=False)
    print(f"  Grid_T 변수: {list(ds_grid.data_vars)}")

    # 시간 범위 선택
    # 5d 데이터: 1년에 73 timestep (365/5)
    time_dim = 'time_counter' if 'time_counter' in ds_grid.dims else 'time'
    n_time = len(ds_grid[time_dim])
    print(f"  총 timestep: {n_time}")

    # 5일 간격이면 연도별 인덱스 계산
    steps_per_year = 73  # 365/5 = 73 (NEMO에서는 보통 72 또는 73)
    idx_start = (YEAR_START - 1) * steps_per_year
    idx_end = YEAR_END * steps_per_year
    if idx_end > n_time:
        idx_end = n_time
        print(f"  WARNING: 요청 범위 초과. {n_time} timesteps까지 사용.")

    print(f"  5d 시간 선택: index {idx_start}-{idx_end-1}")
    temp = ds_grid[VAR_TEMP].isel({time_dim: slice(idx_start, idx_end)})

    # 깊이 축 확인
    depth_coord = VAR_DEPTH_AXIS
    if depth_coord not in temp.dims:
        for d in temp.dims:
            if 'depth' in d.lower():
                depth_coord = d
                break
    depths = temp[depth_coord].values
    print(f"  깊이 레벨: {len(depths)} levels, range {depths[0]:.1f}-{depths[-1]:.1f} m")

    # Tp: 표층 0-100m 평균
    print("  Tp (0-100m) 계산...")
    Tp = depth_average(temp, depth_coord, 0, 100).mean(dim=time_dim).compute()
    print(f"    Tp range: {float(Tp.min()):.2f} ~ {float(Tp.max()):.2f} °C")

    # Tm: 중층 500-1500m 평균
    print("  Tm (500-1500m) 계산...")
    if depths[-1] >= 500:
        Tm = depth_average(temp, depth_coord, 500, 1500).mean(dim=time_dim).compute()
        print(f"    Tm range: {float(Tm.min()):.2f} ~ {float(Tm.max()):.2f} °C")
    else:
        print("  WARNING: 깊이가 500m 미만. Tm을 최저층 온도로 대체.")
        Tm = temp.isel({depth_coord: -1}).mean(dim=time_dim).compute()

    # Tb: 저층 (해저) 온도
    print("  Tb (bottom) 계산...")
    temp_tmean = temp.mean(dim=time_dim).compute()
    Tb = get_bottom_values(temp_tmean, depth_coord)
    print(f"    Tb range: {float(Tb.min()):.2f} ~ {float(Tb.max()):.2f} °C")

    ds_grid.close()

    # ------------------------------------------------------------------
    # 3. 생지화학 변수 처리 (ptrc_T, diad_T)
    # ------------------------------------------------------------------
    print("\n[3] 생지화학 변수 처리...")

    # --- 월평균 또는 연평균 데이터 사용 ---
    # 먼저 연평균 시도, 없으면 월평균 사용
    diad_path = os.path.join(DATA_DIR, FILE_DIAD_1Y)
    ptrc_path = os.path.join(DATA_DIR, FILE_PTRC_1Y)
    use_yearly = True

    if not os.path.exists(diad_path):
        diad_path = os.path.join(DATA_DIR, FILE_DIAD_1M)
        use_yearly = False
    if not os.path.exists(ptrc_path):
        ptrc_path = os.path.join(DATA_DIR, FILE_PTRC_1M)
        use_yearly = False

    ds_diad = xr.open_dataset(diad_path, chunks='auto', decode_times=False)
    ds_ptrc = xr.open_dataset(ptrc_path, chunks='auto', decode_times=False)

    print(f"  DIAD 변수: {list(ds_diad.data_vars)}")
    print(f"  PTRC 변수: {list(ds_ptrc.data_vars)}")

    # 시간 범위 선택
    if use_yearly:
        ds_diad = select_time_range_yearly(ds_diad, YEAR_START, YEAR_END)
        ds_ptrc = select_time_range_yearly(ds_ptrc, YEAR_START, YEAR_END)
    else:
        ds_diad = select_time_range(ds_diad, YEAR_START, YEAR_END)
        ds_ptrc = select_time_range(ds_ptrc, YEAR_START, YEAR_END)

    # ---------------------------------------------------------------
    # 동물플랑크톤 생산 (lzprod, szprod)
    # ---------------------------------------------------------------
    # 옵션 A: diad_T에 직접적인 생산 변수가 있는 경우
    # 옵션 B: ptrc_T의 동물플랑크톤 biomass에서 추정
    # 실제 변수명은 01_inspect 결과를 보고 수정

    lzprod = None
    szprod = None
    dfbot = None
    photic = None

    time_dim_bio = 'time_counter' if 'time_counter' in ds_diad.dims else 'time'

    if VAR_LZPROD and VAR_LZPROD in ds_diad:
        print(f"  대형 zoo 생산: {VAR_LZPROD}")
        lzprod = ds_diad[VAR_LZPROD].mean(dim=time_dim_bio).squeeze().compute()
    elif VAR_LZPROD and VAR_LZPROD in ds_ptrc:
        print(f"  대형 zoo 생산: {VAR_LZPROD} (ptrc)")
        lzprod = ds_ptrc[VAR_LZPROD].mean(dim=time_dim_bio).squeeze().compute()

    if VAR_SZPROD and VAR_SZPROD in ds_diad:
        print(f"  소형 zoo 생산: {VAR_SZPROD}")
        szprod = ds_diad[VAR_SZPROD].mean(dim=time_dim_bio).squeeze().compute()
    elif VAR_SZPROD and VAR_SZPROD in ds_ptrc:
        print(f"  소형 zoo 생산: {VAR_SZPROD} (ptrc)")
        szprod = ds_ptrc[VAR_SZPROD].mean(dim=time_dim_bio).squeeze().compute()

    # ---------------------------------------------------------------
    # Detrital flux (dfbot)
    # ---------------------------------------------------------------
    if VAR_DFBOT and VAR_DFBOT in ds_diad:
        print(f"  Detrital flux: {VAR_DFBOT}")
        dfbot = ds_diad[VAR_DFBOT].mean(dim=time_dim_bio).squeeze().compute()
    elif VAR_DFBOT and VAR_DFBOT in ds_ptrc:
        print(f"  Detrital flux: {VAR_DFBOT} (ptrc)")
        dfbot = ds_ptrc[VAR_DFBOT].mean(dim=time_dim_bio).squeeze().compute()

    # ---------------------------------------------------------------
    # 유광층 깊이 (photic)
    # ---------------------------------------------------------------
    if VAR_PHOTIC and VAR_PHOTIC in ds_diad:
        print(f"  유광층 깊이: {VAR_PHOTIC}")
        photic = ds_diad[VAR_PHOTIC].mean(dim=time_dim_bio).squeeze().compute()
    elif VAR_PHOTIC and VAR_PHOTIC in ds_ptrc:
        print(f"  유광층 깊이: {VAR_PHOTIC} (ptrc)")
        photic = ds_ptrc[VAR_PHOTIC].mean(dim=time_dim_bio).squeeze().compute()

    ds_diad.close()
    ds_ptrc.close()

    # ------------------------------------------------------------------
    # 4. 변수 확인 및 대체 방법
    # ------------------------------------------------------------------
    print("\n[4] 변수 상태 확인...")

    missing_vars = []
    if lzprod is None:
        missing_vars.append("lzprod (대형 동물플랑크톤 생산)")
    if szprod is None:
        missing_vars.append("szprod (소형 동물플랑크톤 생산)")
    if dfbot is None:
        missing_vars.append("dfbot (해저 detrital flux)")
    if photic is None:
        missing_vars.append("photic (유광층 깊이)")
    if depth_bottom is None:
        missing_vars.append("depth (해저 깊이)")

    if missing_vars:
        print("  *** 아래 변수를 찾지 못했습니다 ***")
        for v in missing_vars:
            print(f"    - {v}")
        print("\n  >> 01_inspect_nemo_output.py 결과를 확인하고")
        print("  >> 이 스크립트 상단 CONFIG 섹션의 변수명을 수정하세요.")
        print("  >> 또는 'placeholder' 값을 사용하여 일단 진행합니다.\n")

    # Placeholder: 누락된 변수는 NaN으로 채움
    ny, nx = lon2d.shape if lon2d.ndim == 2 else (len(lon2d), 1)
    shape_2d = (ny, nx)

    if lzprod is None:
        lzprod_vals = np.full(shape_2d, np.nan)
    else:
        lzprod_vals = np.array(lzprod.values).reshape(shape_2d)

    if szprod is None:
        szprod_vals = np.full(shape_2d, np.nan)
    else:
        szprod_vals = np.array(szprod.values).reshape(shape_2d)

    if dfbot is None:
        dfbot_vals = np.full(shape_2d, np.nan)
    else:
        dfbot_vals = np.array(dfbot.values).reshape(shape_2d)

    if photic is None:
        photic_vals = np.full(shape_2d, np.nan)
    else:
        photic_vals = np.array(photic.values).reshape(shape_2d)

    if depth_bottom is None:
        depth_vals = np.full(shape_2d, np.nan)
    else:
        depth_vals = np.array(depth_bottom.values).reshape(shape_2d)

    Tp_vals = np.array(Tp.values).reshape(shape_2d)
    Tm_vals = np.array(Tm.values).reshape(shape_2d)
    Tb_vals = np.array(Tb.values).reshape(shape_2d)

    if lon2d.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon2d, lat2d)

    # ------------------------------------------------------------------
    # 5. DataFrame 생성 및 바다 격자만 선택
    # ------------------------------------------------------------------
    print("\n[5] DataFrame 생성...")

    df = pd.DataFrame({
        'lon': lon2d.ravel(),
        'lat': lat2d.ravel(),
        'lzprod': lzprod_vals.ravel(),
        'szprod': szprod_vals.ravel(),
        'dfbot': dfbot_vals.ravel(),
        'photic': photic_vals.ravel(),
        'depth': depth_vals.ravel(),
        'Tb': Tb_vals.ravel(),
        'Tm': Tm_vals.ravel(),
        'Tp': Tp_vals.ravel(),
    })

    # 육지 제거 (온도가 NaN인 격자)
    n_total = len(df)
    df = df.dropna(subset=['Tp'])
    n_ocean = len(df)
    print(f"  전체 격자: {n_total}, 바다 격자: {n_ocean}")

    # ------------------------------------------------------------------
    # 6. CSV 저장
    # ------------------------------------------------------------------
    print("\n[6] CSV 저장...")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outpath = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    df.to_csv(outpath, index=False)
    print(f"  저장 완료: {outpath}")
    print(f"  격자 수: {len(df)}")
    print(f"\n  데이터 요약:")
    print(df.describe().to_string())

    # ------------------------------------------------------------------
    # 7. 빠른 검증 플롯 (선택)
    # ------------------------------------------------------------------
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 3, figsize=(18, 12))
        plot_vars = ['Tp', 'Tm', 'Tb', 'depth', 'photic', 'lzprod', 'szprod', 'dfbot']

        for ax, vname in zip(axes.ravel(), plot_vars):
            sc = ax.scatter(df['lon'], df['lat'], c=df[vname],
                            s=0.5, cmap='viridis', edgecolors='none')
            ax.set_title(vname)
            plt.colorbar(sc, ax=ax)

        axes.ravel()[-1].axis('off')
        plt.tight_layout()
        fig_path = os.path.join(OUTPUT_DIR, "Input_NEMO_FEISTY_check.png")
        plt.savefig(fig_path, dpi=150)
        plt.close()
        print(f"\n  검증 플롯: {fig_path}")
    except ImportError:
        print("\n  matplotlib 미설치 - 검증 플롯 생략")

    print("\n" + "=" * 70)
    print("완료!")
    print("=" * 70)


if __name__ == "__main__":
    main()
