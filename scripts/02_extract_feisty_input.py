"""
02_extract_feisty_input.py
==========================
NEMO-PISCES NetCDF 출력에서 FEISTY R 패키지 입력 형식(Input_global.csv)으로 변환.

10년 평균을 계산하여 climatology 형태의 입력 데이터를 생성.

변수 매핑 (ncdump -h 결과 확인 완료):
  - Tp (0-100m 수온)    ← grid_T 5d: thetao, depth-averaged 0-100m
  - Tm (500-1500m 수온)  ← grid_T 5d: thetao, depth-averaged 500-1500m
  - Tb (해저 수온)       ← grid_T 5d: sbt (sea bottom temperature)
  - depth (해저 깊이)    ← grid_T_static: depth
  - photic (유광층)      ← diad_T 1y: Heup (euphotic layer depth)
  - szprod (소형 zoo)    ← diad_T 1y: GRAZ1 (micro-zoo grazing, 수직적분)
  - lzprod (대형 zoo)    ← diad_T 1y: GRAZ2 (meso-zoo grazing, 수직적분)
  - dfbot (해저 detritus) ← diad_T 1y: EPC100 + Martin curve

단위 변환:
  GRAZ1/GRAZ2: mol C/m3/s → 수직적분(×dz) → mol C/m2/s → ×sec_per_yr×12×9 → g ww/m2/yr
  EPC100: mol C/m2/s → ×sec_per_yr×12×9 → g ww/m2/yr (at 100m) → Martin curve → bottom

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
FILE_STATIC  = "static_depth.nc"           # depth
FILE_DIAD_1Y = "diad_T_1y_10yr.nc"        # Heup, EPC100, GRAZ1, GRAZ2 (1y, 10 steps)
FILE_DIAD_1M = "diad_T_1m_10yr.nc"        # 같은 변수 월별 (1m, 120 steps)

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
    print(f"(CDO 전처리된 subset 파일 사용)")
    print("=" * 70)

    # ==================================================================
    # 1. Static: 좌표 + 해저 깊이
    # ==================================================================
    print("\n[1] Static 파일 (좌표, 수심)...")
    ds_static = xr.open_dataset(os.path.join(DATA_DIR, FILE_STATIC),
                                decode_times=False)

    lon2d = ds_static['nav_lon'].values   # (y=148, x=180)
    lat2d = ds_static['nav_lat'].values
    depth_bot = ds_static['depth'].values  # sea floor depth (m), FillValue=1e20

    # FillValue → NaN (육지)
    depth_bot = np.where(depth_bot > 1e10, np.nan, depth_bot)
    ocean_mask = ~np.isnan(depth_bot)

    print(f"  격자: {lon2d.shape}, 바다 격자: {ocean_mask.sum()}")
    ds_static.close()

    # ==================================================================
    # 2. Temperature: Tp, Tm, Tb from grid_T (5d)
    # ==================================================================
    print("\n[2] 온도 데이터 (grid_T 5d)...")
    ds_grid = xr.open_dataset(os.path.join(DATA_DIR, FILE_GRID_T),
                              decode_times=False)

    # 깊이 레벨 확인
    deptht = ds_grid['deptht'].values
    print(f"  깊이 레벨 ({len(deptht)}): {np.round(deptht, 1)}")

    # CDO 전처리 후 이미 10년만 포함 → 전체 사용
    n_time = len(ds_grid['time_counter'])
    t_start = 0
    t_end = n_time
    print(f"  시간: {n_time} steps (전체 사용, CDO로 이미 10년 추출됨)")

    # --- Tb: sbt (sea bottom temperature) - 2D, 바로 사용 ---
    print("  Tb (sbt)...")
    sbt = ds_grid['sbt'].isel(time_counter=slice(t_start, t_end))
    Tb = sbt.mean(dim='time_counter').values
    Tb = np.where(np.abs(Tb) > 1e10, np.nan, Tb)
    print(f"    range: {np.nanmin(Tb):.2f} ~ {np.nanmax(Tb):.2f} °C")

    # --- Tp: thetao 0-100m depth-average ---
    print("  Tp (thetao 0-100m)...")
    # 0-100m에 해당하는 깊이 인덱스
    idx_tp = np.where(deptht <= 100)[0]
    print(f"    depth levels used: {deptht[idx_tp]}")

    thetao = ds_grid['thetao'].isel(time_counter=slice(t_start, t_end))

    # 레벨 두께 계산 (deptht_bounds 있으면 사용, 없으면 중간점 방식)
    if 'deptht_bounds' in ds_grid:
        dz = ds_grid['deptht_bounds'].values  # (31, 2)
        layer_thickness = dz[:, 1] - dz[:, 0]
    else:
        # 중간점 방식: 각 레벨의 두께 = 인접 레벨 중간점 간 거리
        layer_thickness = np.zeros(len(deptht))
        layer_thickness[0] = (deptht[0] + deptht[1]) / 2.0
        for k in range(1, len(deptht) - 1):
            layer_thickness[k] = (deptht[k+1] - deptht[k-1]) / 2.0
        layer_thickness[-1] = layer_thickness[-2]  # 마지막 레벨
        print("    (deptht_bounds 없음 → 중간점 방식으로 두께 계산)")
    print(f"    layer thickness (0-100m): {np.round(layer_thickness[idx_tp], 1)}")

    # 가중 평균: Tp = Σ(T_k × dz_k) / Σ(dz_k) for k in 0-100m
    weights_tp = xr.DataArray(layer_thickness[idx_tp], dims=['deptht'])
    temp_tp = thetao.isel(deptht=idx_tp)
    Tp = (temp_tp * weights_tp).sum(dim='deptht') / weights_tp.sum()
    Tp = Tp.mean(dim='time_counter').values
    Tp = np.where(np.abs(Tp) > 1e10, np.nan, Tp)
    print(f"    range: {np.nanmin(Tp):.2f} ~ {np.nanmax(Tp):.2f} °C")

    # --- Tm: thetao 500-1500m depth-average ---
    print("  Tm (thetao 500-1500m)...")
    idx_tm = np.where((deptht >= 500) & (deptht <= 1500))[0]
    print(f"    depth levels used: {np.round(deptht[idx_tm], 1)}")

    if len(idx_tm) > 0:
        weights_tm = xr.DataArray(layer_thickness[idx_tm], dims=['deptht'])
        temp_tm = thetao.isel(deptht=idx_tm)
        Tm = (temp_tm * weights_tm).sum(dim='deptht') / weights_tm.sum()
        Tm = Tm.mean(dim='time_counter').values
        Tm = np.where(np.abs(Tm) > 1e10, np.nan, Tm)
    else:
        print("    WARNING: 500-1500m 레벨 없음. 최하층 사용.")
        Tm = thetao.isel(deptht=-1).mean(dim='time_counter').values
        Tm = np.where(np.abs(Tm) > 1e10, np.nan, Tm)

    # 얕은 해역에서 Tm이 NaN인 경우 Tb로 대체
    Tm = np.where(np.isnan(Tm) & ocean_mask, Tb, Tm)
    print(f"    range: {np.nanmin(Tm):.2f} ~ {np.nanmax(Tm):.2f} °C")

    ds_grid.close()

    # ==================================================================
    # 3. Biogeochemistry from diad_T (yearly)
    # ==================================================================
    print("\n[3] 생지화학 (diad_T 1y)...")
    ds_diad = xr.open_dataset(os.path.join(DATA_DIR, FILE_DIAD_1Y),
                              decode_times=False)

    n_yr = len(ds_diad['time_counter'])
    print(f"  연평균: {n_yr} years (전체 사용, CDO로 이미 10년 추출됨)")
    ds_diad_sel = ds_diad  # 전체 사용

    # --- Heup: 유광층 깊이 (m) ---
    print("  photic (Heup)...")
    photic = ds_diad_sel['Heup'].mean(dim='time_counter').values
    photic = np.where(np.abs(photic) > 1e10, np.nan, photic)
    print(f"    range: {np.nanmin(photic):.1f} ~ {np.nanmax(photic):.1f} m")

    # --- EPC100 → dfbot (Martin curve) ---
    print("  dfbot (EPC100 → Martin curve → bottom)...")
    epc100 = ds_diad_sel['EPC100'].mean(dim='time_counter').values  # mol C/m2/s
    epc100 = np.where(np.abs(epc100) > 1e10, np.nan, epc100)

    # mol C/m2/s → g ww/m2/yr at 100m
    epc100_gww = epc100 * MOLC_M2S_TO_GWW_M2YR

    # Martin curve: F(z_bot) = F(100m) × (z_bot/100)^(-b)
    # depth_bot < 100m → 그냥 100m 값 사용
    depth_for_martin = np.maximum(depth_bot, 100.0)
    dfbot = epc100_gww * np.power(depth_for_martin / 100.0, -MARTIN_B)
    print(f"    EPC100 (g ww/m2/yr at 100m): {np.nanmin(epc100_gww):.2f} ~ {np.nanmax(epc100_gww):.2f}")
    print(f"    dfbot (at seafloor): {np.nanmin(dfbot):.4f} ~ {np.nanmax(dfbot):.2f}")

    # --- GRAZ1 → szprod, GRAZ2 → lzprod (수직적분) ---
    print("  szprod (GRAZ1 수직적분)...")
    graz1_3d = ds_diad_sel['GRAZ1'].mean(dim='time_counter')  # (deptht, y, x) mol/m3/s
    graz2_3d = ds_diad_sel['GRAZ2'].mean(dim='time_counter')

    # 수직적분: Σ GRAZ(k) × dz(k) → mol C/m2/s
    # layer_thickness는 위에서 deptht_bounds로 계산
    dz_weights = xr.DataArray(layer_thickness, dims=['deptht'],
                              coords={'deptht': ds_diad['deptht'].values})

    graz1_int = (graz1_3d * dz_weights).sum(dim='deptht').values  # mol C/m2/s
    graz2_int = (graz2_3d * dz_weights).sum(dim='deptht').values

    # FillValue 처리
    graz1_int = np.where(np.abs(graz1_int) > 1e10, np.nan, graz1_int)
    graz2_int = np.where(np.abs(graz2_int) > 1e10, np.nan, graz2_int)

    # mol C/m2/s → g ww/m2/yr
    szprod = graz1_int * MOLC_M2S_TO_GWW_M2YR
    lzprod = graz2_int * MOLC_M2S_TO_GWW_M2YR

    print(f"    szprod: {np.nanmin(szprod):.2f} ~ {np.nanmax(szprod):.2f} g ww/m2/yr")
    print(f"  lzprod (GRAZ2 수직적분)...")
    print(f"    lzprod: {np.nanmin(lzprod):.2f} ~ {np.nanmax(lzprod):.2f} g ww/m2/yr")

    ds_diad.close()

    # ==================================================================
    # 4. DataFrame 생성 및 CSV 저장
    # ==================================================================
    print("\n[4] DataFrame 생성...")

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

    # 육지 제거 (Tp가 NaN인 곳)
    n_total = len(df)
    df = df.dropna(subset=['Tp'])
    n_ocean = len(df)
    print(f"  전체: {n_total}, 바다: {n_ocean}")

    # 음수 값 보정 (생산량은 0 이상이어야 함)
    for col in ['lzprod', 'szprod', 'dfbot']:
        n_neg = (df[col] < 0).sum()
        if n_neg > 0:
            print(f"  WARNING: {col}에 음수 {n_neg}개 → 0으로 처리")
            df[col] = df[col].clip(lower=0)

    # 저장
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    outpath = os.path.join(OUTPUT_DIR, OUTPUT_FILE)
    df.to_csv(outpath, index=False)
    print(f"\n  저장: {outpath}")
    print(f"  격자 수: {len(df)}")

    # 요약 통계
    print(f"\n{'='*60}")
    print("데이터 요약")
    print(f"{'='*60}")
    print(df.describe().round(3).to_string())

    # 기존 Input_global.csv와 비교
    ref_path = os.path.join(OUTPUT_DIR, "Input_global.csv")
    if os.path.exists(ref_path):
        df_ref = pd.read_csv(ref_path)
        print(f"\n{'='*60}")
        print("기존 Input_global.csv (참고용 비교)")
        print(f"{'='*60}")
        print(f"  격자 수: {len(df_ref)} (NEMO: {len(df)})")
        print(df_ref.describe().round(3).to_string())

    # ==================================================================
    # 5. 검증 플롯
    # ==================================================================
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        fig = plt.figure(figsize=(20, 14))
        gs = gridspec.GridSpec(3, 3, hspace=0.35, wspace=0.3)

        plot_info = [
            ('Tp', 'Tp (0-100m) [°C]', 'RdYlBu_r'),
            ('Tm', 'Tm (500-1500m) [°C]', 'RdYlBu_r'),
            ('Tb', 'Tb (bottom) [°C]', 'RdYlBu_r'),
            ('depth', 'Seafloor depth [m]', 'viridis_r'),
            ('photic', 'Euphotic depth [m]', 'YlGn_r'),
            ('lzprod', 'lzprod (GRAZ2) [g ww/m²/yr]', 'YlOrRd'),
            ('szprod', 'szprod (GRAZ1) [g ww/m²/yr]', 'YlOrRd'),
            ('dfbot', 'dfbot (seafloor) [g ww/m²/yr]', 'YlOrBr'),
        ]

        for i, (var, title, cmap) in enumerate(plot_info):
            ax = fig.add_subplot(gs[i // 3, i % 3])
            vals = df[var].values
            vmin, vmax = np.nanpercentile(vals[np.isfinite(vals)], [2, 98])
            sc = ax.scatter(df['lon'], df['lat'], c=vals,
                            s=0.3, cmap=cmap, vmin=vmin, vmax=vmax,
                            edgecolors='none', rasterized=True)
            ax.set_title(title, fontsize=10)
            ax.set_xlim(-180, 180)
            ax.set_ylim(-90, 90)
            plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

        # 마지막 칸: 비교 히스토그램
        ax = fig.add_subplot(gs[2, 2])
        if os.path.exists(ref_path):
            df_ref = pd.read_csv(ref_path)
            ax.hist(df_ref['lzprod'], bins=50, alpha=0.5, label='Original (COBALT)',
                    density=True, color='blue')
            ax.hist(df['lzprod'].dropna(), bins=50, alpha=0.5, label='NEMO-PISCES',
                    density=True, color='red')
            ax.set_xlabel('lzprod [g ww/m²/yr]')
            ax.set_title('lzprod 분포 비교')
            ax.legend(fontsize=8)
        else:
            ax.text(0.5, 0.5, 'No reference\nInput_global.csv', ha='center',
                    va='center', transform=ax.transAxes)

        fig.suptitle('NEMO-PISCES → FEISTY Input (10-year mean)',
                     fontsize=14, fontweight='bold')

        fig_path = os.path.join(OUTPUT_DIR, "Input_NEMO_FEISTY_check.png")
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"\n  검증 플롯: {fig_path}")

    except ImportError as e:
        print(f"\n  matplotlib 미설치 - 검증 플롯 생략 ({e})")

    print("\n" + "=" * 70)
    print("완료! FEISTY 실행:")
    print('  glob <- read.csv("data/Input_NEMO_FEISTY.csv")')
    print("=" * 70)


if __name__ == "__main__":
    main()
