"""
02_extract_feisty_input.py
==========================
NEMO-PISCES NetCDF 출력에서 FEISTY R 패키지 입력 형식(Input_global.csv)으로 변환.

10년 평균을 계산하여 climatology 형태의 입력 데이터를 생성.

변수 매핑:
  - Tp (0-100m 수온)    ← grid_T 5d: thetao, depth-averaged 0-100m
  - Tm (500-1500m 수온)  ← grid_T 5d: thetao, depth-averaged 500-1500m
  - Tb (해저 수온)       ← grid_T 5d: sbt (sea bottom temperature)
  - depth (해저 깊이)    ← grid_T 5d: thetao 3D mask에서 재구성 (deepest valid level)
  - photic (유광층)      ← diad_T 1y: Heup (euphotic layer depth)
  - szprod (소형 zoo)    ← diad_T 1y: GRAZ1 (micro-zoo grazing, 수직적분)
  - lzprod (대형 zoo)    ← diad_T 1y: GRAZ2 (meso-zoo grazing, 수직적분)
  - dfbot (해저 detritus) ← diad_T 1y: EPC100 + Martin curve

NOTE: static depth 파일은 NEMO에서 time_counter=0으로 출력되어 값이 비어있음.
      따라서 thetao의 3D wet mask로부터 해저 깊이를 재구성함.

사용법:
    python 02_extract_feisty_input.py [DATA_DIR]
"""

import os
import sys
import numpy as np
import xarray as xr
import pandas as pd

# ============================================================================
# CONFIG
# ============================================================================

DATA_DIR = sys.argv[1] if len(sys.argv) > 1 else "/data01/labdisk/sungjin/NEMO-FEISTY/subset"

# 출력 경로
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "..", "quarto", "data")
OUTPUT_FILE = "Input_NEMO_FEISTY.csv"

# 파일명 - CDO 전처리 후 subset 파일 (00_preprocess_cdo.sh로 생성)
FILE_GRID_T  = "grid_T_10yr.nc"           # thetao, sbt (5d, 730 steps)
FILE_DIAD_1Y = "diad_T_1y_10yr.nc"        # Heup, EPC100, GRAZ1, GRAZ2 (1y, 10 steps)

# 단위 변환 상수 (noleap calendar)
SEC_PER_YEAR = 365 * 86400  # 31,536,000 s
MOL_C_TO_G_C = 12.0         # g C per mol C
WW_TO_C_RATIO = 9.0         # g wet weight per g C

# mol C/m2/s → g ww/m2/yr
MOLC_M2S_TO_GWW_M2YR = SEC_PER_YEAR * MOL_C_TO_G_C * WW_TO_C_RATIO  # = 3,405,888,000

# Martin curve exponent for POC flux attenuation
MARTIN_B = 0.858


# ============================================================================
# MAIN
# ============================================================================

def main():
    print("=" * 70)
    print("NEMO-PISCES → FEISTY Input 변환")
    print(f"데이터: {DATA_DIR}")
    print("=" * 70)

    # ==================================================================
    # 1. grid_T: 좌표, 깊이축, 온도 (Tp, Tm, Tb) + bathymetry 재구성
    # ==================================================================
    print("\n[1] grid_T 로드...")
    ds_grid = xr.open_dataset(os.path.join(DATA_DIR, FILE_GRID_T),
                              decode_times=False)

    # 좌표
    lon2d = ds_grid['nav_lon'].values   # (y=148, x=180)
    lat2d = ds_grid['nav_lat'].values
    print(f"  격자: {lon2d.shape}")

    # 깊이 레벨
    deptht = ds_grid['deptht'].values
    print(f"  깊이 레벨 ({len(deptht)}): {np.round(deptht, 1)}")

    n_time = len(ds_grid['time_counter'])
    print(f"  시간: {n_time} steps")

    # 레벨 두께 계산
    if 'deptht_bounds' in ds_grid:
        dz = ds_grid['deptht_bounds'].values  # (31, 2)
        layer_thickness = dz[:, 1] - dz[:, 0]
        # 각 레벨 하한 (bathymetry 추정용)
        level_bottom = dz[:, 1]
    else:
        # 중간점 방식
        layer_thickness = np.zeros(len(deptht))
        layer_thickness[0] = (deptht[0] + deptht[1]) / 2.0
        for k in range(1, len(deptht) - 1):
            layer_thickness[k] = (deptht[k+1] - deptht[k-1]) / 2.0
        layer_thickness[-1] = layer_thickness[-2]
        # 레벨 하한 근사
        level_bottom = deptht + layer_thickness / 2.0
        print("  (deptht_bounds 없음 → 중간점 방식)")

    print(f"  layer thickness: {np.round(layer_thickness, 1)}")

    # ------------------------------------------------------------------
    # 1a. Bathymetry 재구성: thetao의 deepest valid level
    #     static depth 파일이 비어있으므로 (time_counter=0 문제)
    #     thetao 3D mask에서 해저 깊이를 추정
    # ------------------------------------------------------------------
    print("\n  [1a] Bathymetry 재구성 (thetao 3D mask)...")

    # 첫 번째 timestep의 thetao로 3D mask 결정
    thetao_t0 = ds_grid['thetao'].isel(time_counter=0).values  # (deptht, y, x)
    # FillValue(1e20) → NaN
    thetao_t0 = np.where(np.abs(thetao_t0) > 1e10, np.nan, thetao_t0)

    # 각 (y, x)에서 유효한 가장 깊은 레벨의 하한 = bathymetry
    valid_mask = ~np.isnan(thetao_t0)  # (31, y, x)

    # 가장 깊은 유효 레벨 인덱스 찾기
    ny, nx = lon2d.shape
    depth_bot = np.full((ny, nx), np.nan)

    for k in range(len(deptht) - 1, -1, -1):
        # 아직 depth가 할당되지 않은 격자 중 이 레벨에서 유효한 곳
        update = np.isnan(depth_bot) & valid_mask[k]
        depth_bot[update] = level_bottom[k]

    ocean_mask = ~np.isnan(depth_bot)
    print(f"    바다 격자: {ocean_mask.sum()}")
    print(f"    depth range: {np.nanmin(depth_bot):.1f} ~ {np.nanmax(depth_bot):.1f} m")

    # ------------------------------------------------------------------
    # 1b. Tb: sbt (sea bottom temperature)
    # ------------------------------------------------------------------
    print("\n  [1b] Tb (sbt)...")
    Tb = ds_grid['sbt'].mean(dim='time_counter').values
    Tb = np.where(np.abs(Tb) > 1e10, np.nan, Tb)
    print(f"    range: {np.nanmin(Tb):.2f} ~ {np.nanmax(Tb):.2f} °C")

    # ------------------------------------------------------------------
    # 1c. Tp: thetao 0-100m 깊이 가중 평균
    # ------------------------------------------------------------------
    print("\n  [1c] Tp (thetao 0-100m)...")
    idx_tp = np.where(deptht <= 100)[0]
    print(f"    levels: {np.round(deptht[idx_tp], 1)}")

    thetao = ds_grid['thetao']
    weights_tp = xr.DataArray(layer_thickness[idx_tp], dims=['deptht'])
    temp_tp = thetao.isel(deptht=idx_tp)
    Tp = (temp_tp * weights_tp).sum(dim='deptht') / weights_tp.sum()
    Tp = Tp.mean(dim='time_counter').values
    Tp = np.where(np.abs(Tp) > 1e10, np.nan, Tp)
    print(f"    range: {np.nanmin(Tp):.2f} ~ {np.nanmax(Tp):.2f} °C")

    # ------------------------------------------------------------------
    # 1d. Tm: thetao 500-1500m 깊이 가중 평균
    # ------------------------------------------------------------------
    print("\n  [1d] Tm (thetao 500-1500m)...")
    idx_tm = np.where((deptht >= 500) & (deptht <= 1500))[0]
    print(f"    levels: {np.round(deptht[idx_tm], 1)}")

    if len(idx_tm) > 0:
        weights_tm = xr.DataArray(layer_thickness[idx_tm], dims=['deptht'])
        temp_tm = thetao.isel(deptht=idx_tm)
        Tm = (temp_tm * weights_tm).sum(dim='deptht') / weights_tm.sum()
        Tm = Tm.mean(dim='time_counter').values
        Tm = np.where(np.abs(Tm) > 1e10, np.nan, Tm)
    else:
        print("    WARNING: 500-1500m 레벨 없음 → 최하층 사용")
        Tm = thetao.isel(deptht=-1).mean(dim='time_counter').values
        Tm = np.where(np.abs(Tm) > 1e10, np.nan, Tm)

    # 얕은 해역 (500m 미만): Tm NaN → Tb로 대체
    Tm = np.where(np.isnan(Tm) & ocean_mask, Tb, Tm)
    print(f"    range: {np.nanmin(Tm):.2f} ~ {np.nanmax(Tm):.2f} °C")

    ds_grid.close()

    # ==================================================================
    # 2. Biogeochemistry from diad_T (yearly)
    # ==================================================================
    print("\n[2] 생지화학 (diad_T 1y)...")
    ds_diad = xr.open_dataset(os.path.join(DATA_DIR, FILE_DIAD_1Y),
                              decode_times=False)

    n_yr = len(ds_diad['time_counter'])
    print(f"  {n_yr} years")

    # --- Heup: 유광층 깊이 (m) ---
    print("  photic (Heup)...")
    photic = ds_diad['Heup'].mean(dim='time_counter').values
    photic = np.where(np.abs(photic) > 1e10, np.nan, photic)
    print(f"    range: {np.nanmin(photic):.1f} ~ {np.nanmax(photic):.1f} m")

    # --- EPC100 → dfbot (Martin curve) ---
    print("  dfbot (EPC100 → Martin curve → bottom)...")
    epc100 = ds_diad['EPC100'].mean(dim='time_counter').values  # mol C/m2/s
    epc100 = np.where(np.abs(epc100) > 1e10, np.nan, epc100)

    epc100_gww = epc100 * MOLC_M2S_TO_GWW_M2YR  # g ww/m2/yr at 100m
    depth_for_martin = np.maximum(depth_bot, 100.0)
    dfbot = epc100_gww * np.power(depth_for_martin / 100.0, -MARTIN_B)
    print(f"    EPC100 at 100m: {np.nanmin(epc100_gww):.2f} ~ {np.nanmax(epc100_gww):.2f} g ww/m2/yr")
    print(f"    dfbot at bottom: {np.nanmin(dfbot):.4f} ~ {np.nanmax(dfbot):.2f} g ww/m2/yr")

    # --- GRAZ1 → szprod, GRAZ2 → lzprod (수직적분) ---
    print("  szprod/lzprod (GRAZ1/GRAZ2 수직적분)...")
    graz1_3d = ds_diad['GRAZ1'].mean(dim='time_counter')  # (deptht, y, x)
    graz2_3d = ds_diad['GRAZ2'].mean(dim='time_counter')

    dz_weights = xr.DataArray(layer_thickness, dims=['deptht'],
                              coords={'deptht': ds_diad['deptht'].values})

    graz1_int = (graz1_3d * dz_weights).sum(dim='deptht').values
    graz2_int = (graz2_3d * dz_weights).sum(dim='deptht').values

    graz1_int = np.where(np.abs(graz1_int) > 1e10, np.nan, graz1_int)
    graz2_int = np.where(np.abs(graz2_int) > 1e10, np.nan, graz2_int)

    szprod = graz1_int * MOLC_M2S_TO_GWW_M2YR
    lzprod = graz2_int * MOLC_M2S_TO_GWW_M2YR

    print(f"    szprod: {np.nanmin(szprod):.2f} ~ {np.nanmax(szprod):.2f} g ww/m2/yr")
    print(f"    lzprod: {np.nanmin(lzprod):.2f} ~ {np.nanmax(lzprod):.2f} g ww/m2/yr")

    ds_diad.close()

    # ==================================================================
    # 3. DataFrame 생성 및 CSV 저장
    # ==================================================================
    print("\n[3] DataFrame → CSV...")

    df = pd.DataFrame({
        'lon':    lon2d.ravel(),
        'lat':    lat2d.ravel(),
        'lzprod': lzprod.ravel(),
        'szprod': szprod.ravel(),
        'dfbot':  dfbot.ravel(),
        'photic': photic.ravel(),
        'depth':  depth_bot.ravel(),
        'Tb':     Tb.ravel(),
        'Tm':     Tm.ravel(),
        'Tp':     Tp.ravel(),
    })

    # 육지 제거 (ocean_mask 기반, 즉 thetao가 유효한 곳)
    n_total = len(df)
    df = df.dropna(subset=['Tp', 'depth'])
    n_ocean = len(df)
    print(f"  전체: {n_total}, 바다: {n_ocean}")

    # 음수 보정
    for col in ['lzprod', 'szprod', 'dfbot']:
        n_neg = (df[col] < 0).sum()
        if n_neg > 0:
            print(f"  {col}: 음수 {n_neg}개 → 0")
            df[col] = df[col].clip(lower=0)

    # 저장
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outpath = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    df.to_csv(outpath, index=False)
    print(f"\n  저장: {outpath}")
    print(f"  격자 수: {len(df)}")

    # 요약
    print(f"\n{'='*60}")
    print("데이터 요약")
    print(f"{'='*60}")
    print(df.describe().round(3).to_string())

    # 기존 Input_global.csv와 비교
    ref_path = os.path.join(OUTPUT_DIR, "Input_global.csv")
    if os.path.exists(ref_path):
        df_ref = pd.read_csv(ref_path)
        print(f"\n{'='*60}")
        print("기존 Input_global.csv (참고)")
        print(f"{'='*60}")
        print(f"  격자 수: {len(df_ref)} (NEMO: {len(df)})")
        print(df_ref.describe().round(3).to_string())

    # ==================================================================
    # 4. 검증 플롯
    # ==================================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 4, figsize=(20, 8))
        plot_info = [
            ('Tp', 'Tp (0-100m) [°C]', 'RdYlBu_r'),
            ('Tm', 'Tm (500-1500m) [°C]', 'RdYlBu_r'),
            ('Tb', 'Tb (bottom) [°C]', 'RdYlBu_r'),
            ('depth', 'Seafloor depth [m]', 'viridis_r'),
            ('photic', 'Euphotic depth [m]', 'YlGn_r'),
            ('lzprod', 'lzprod [g ww/m²/yr]', 'YlOrRd'),
            ('szprod', 'szprod [g ww/m²/yr]', 'YlOrRd'),
            ('dfbot', 'dfbot [g ww/m²/yr]', 'YlOrBr'),
        ]

        for ax, (var, title, cmap) in zip(axes.ravel(), plot_info):
            vals = df[var].values
            finite = vals[np.isfinite(vals)]
            if len(finite) == 0:
                ax.set_title(f"{title} (NO DATA)")
                continue
            vmin, vmax = np.nanpercentile(finite, [2, 98])
            sc = ax.scatter(df['lon'], df['lat'], c=vals, s=0.3,
                            cmap=cmap, vmin=vmin, vmax=vmax,
                            edgecolors='none', rasterized=True)
            ax.set_title(title, fontsize=10)
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

        fig.suptitle('NEMO-PISCES → FEISTY Input (10-year mean)',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()

        fig_path = os.path.join(OUTPUT_DIR, "Input_NEMO_FEISTY_check.png")
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n  검증 플롯: {fig_path}")

    except ImportError as e:
        print(f"\n  matplotlib 미설치 - 플롯 생략 ({e})")

    print("\n" + "=" * 70)
    print("완료! FEISTY 실행:")
    print('  glob <- read.csv("data/Input_NEMO_FEISTY.csv")')
    print("=" * 70)


if __name__ == "__main__":
    main()
