#!/usr/bin/env python3
"""
NEMO-PISCES → FEISTY Input CSV 변환
domcfg 기반 bathymetry + tmask + e3t_0 사용
"""

import os
import xarray as xr
import numpy as np
import pandas as pd

# =============================================================
# CONFIG
# =============================================================
BASE_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY"
SUBSET_DIR = f"{BASE_DIR}/02.subset"
OUT_DIR = f"{BASE_DIR}/03.feisty_input"
os.makedirs(OUT_DIR, exist_ok=True)

grid_file = f"{SUBSET_DIR}/grid_T_10yr.nc"
diad_file = f"{SUBSET_DIR}/diad_T_1m_10yr.nc"
zoo_file = f"{BASE_DIR}/01.org/ORCA2_1y_00010101_00501231_ptrc_T.nc"
domcfg_file = f"{BASE_DIR}/01.org/ORCA_R2_zps_domcfg.nc"

# 단위 변환: mol C/m2/s → g wet weight/m2/yr
CONV = 365 * 86400 * 12.0 * 9.0  # = 3,405,888,000

# Closure term 단위 변환: mmol C/m3/day → (적분 후) g ww/m2/yr
# = 1e-3(mmol→mol) × 365(day→yr) × 12(mol C→g C) × 9(g C→g ww)
CONV_CLOSURE = 1e-3 * 365 * 12.0 * 9.0  # = 39.42

# PISCES 2차 사망률 계수 (namelist_pisces_ref)
MZRAT = 0.02    # microzooplankton (/day/(mmol/m3))
MZRAT2 = 0.01   # mesozooplankton  (/day/(mmol/m3))

# Martin curve exponent
MARTIN_B = 0.858

print("=" * 60)
print("NEMO-PISCES → FEISTY Input (domcfg version)")
print(f"grid:   {grid_file}")
print(f"diad:   {diad_file}")
print(f"zoo:    {zoo_file}")
print(f"domcfg: {domcfg_file}")
print(f"MZRAT={MZRAT}, MZRAT2={MZRAT2}")
print("=" * 60)

# =============================================================
# [0] 파일 로드
# =============================================================
print("\n[0] 파일 로드...")
ds_grid = xr.open_dataset(grid_file, decode_times=False)
ds_diad = xr.open_dataset(diad_file, decode_times=False)
ds_dom = xr.open_dataset(domcfg_file)

# 시간 차원 자동 감지
time_dim_grid = [d for d in ds_grid["thetao"].dims if "time" in d][0]
time_dim_diad = [d for d in ds_diad["Heup"].dims if "time" in d][0]

print(f"  grid time: {time_dim_grid}, steps: {len(ds_grid[time_dim_grid])}")
print(f"  diad time: {time_dim_diad}, steps: {len(ds_diad[time_dim_diad])}")

lon2d = ds_grid["nav_lon"].values
lat2d = ds_grid["nav_lat"].values
deptht = ds_grid["deptht"].values
ny, nx = lon2d.shape
nz = len(deptht)

print(f"  격자: ({ny}, {nx}), 총 {ny*nx}")
print(f"  deptht ({nz} levels): {np.round(deptht, 1)}")

# =============================================================
# [1] domcfg: bottom_level, e3t_0, tmask
# =============================================================
print("\n[1] domcfg 로드...")
bottom_level = ds_dom["bottom_level"].values  # (y, x), 0=land
e3t_0 = ds_dom["e3t_0"].values               # (z, y, x) 실제 셀 두께

print(f"  bottom_level shape: {bottom_level.shape}")
print(f"  bottom_level range: {bottom_level.min()} ~ {bottom_level.max()}")
print(f"  e3t_0 shape: {e3t_0.shape}")

# ocean mask
ocean_mask = bottom_level > 0
n_ocean = ocean_mask.sum()
print(f"  바다 격자: {n_ocean}")

# tmask 생성: k < bottom_level → 유효(1), 아니면 0
k3d = np.arange(nz)[:, None, None]  # (z, 1, 1)
tmask = (k3d < bottom_level[None, :, :]).astype(np.float64)  # (z, y, x)

print(f"  tmask shape: {tmask.shape}")
print(f"  tmask per level: {[int(tmask[k].sum()) for k in range(nz)]}")

ds_dom.close()

# =============================================================
# [2] Bathymetry (depth) 계산
# =============================================================
print("\n[2] Bathymetry 계산...")

# 누적 두께 → 각 레벨 하단 깊이
cumsum_e3t = np.cumsum(e3t_0, axis=0)  # (z, y, x)

# depth_bot = cumsum_e3t[bottom_level - 1] (fancy indexing)
yy, xx = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
bl_idx = np.clip(bottom_level - 1, 0, nz - 1)
depth_bot = cumsum_e3t[bl_idx, yy, xx]
depth_bot[~ocean_mask] = np.nan

print(f"  depth range: {np.nanmin(depth_bot):.1f} ~ {np.nanmax(depth_bot):.1f} m")
print(f"  depth mean:  {np.nanmean(depth_bot):.1f} m")

# =============================================================
# [3] Tb: sbt (sea bottom temperature)
# =============================================================
print("\n[3] Tb (sbt)...")

Tb = ds_grid["sbt"].mean(dim=time_dim_grid).values
Tb = np.where(np.abs(Tb) > 1e10, np.nan, Tb)
Tb[~ocean_mask] = np.nan

print(f"  valid: {np.count_nonzero(~np.isnan(Tb))}")
print(f"  range: {np.nanmin(Tb):.2f} ~ {np.nanmax(Tb):.2f} C")

# =============================================================
# [4] thetao 시간평균 + tmask 적용
# =============================================================
print("\n[4] thetao 시간평균...")

thetao_mean = ds_grid["thetao"].mean(dim=time_dim_grid).values  # (z, y, x)
thetao_mean = np.where(np.abs(thetao_mean) > 1e10, np.nan, thetao_mean)

# tmask 적용: 해저 아래 → NaN
thetao_masked = np.where(tmask == 1, thetao_mean, np.nan)

print(f"  thetao_masked shape: {thetao_masked.shape}")

# =============================================================
# [5] Tp: 0-100m 깊이가중 평균 (e3t_0 사용)
# =============================================================
print("\n[5] Tp (thetao 0-100m)...")

idx_tp = np.where(deptht <= 100)[0]
print(f"  levels: {idx_tp} → deptht = {np.round(deptht[idx_tp], 1)}")

# weight = e3t_0 × tmask (해저 아래는 두께 0)
e3t_tp = e3t_0[idx_tp] * tmask[idx_tp]  # (n_lev, y, x)
theta_tp = thetao_masked[idx_tp]

weight_sum = np.nansum(e3t_tp, axis=0)
weighted_theta = np.nansum(theta_tp * e3t_tp, axis=0)

Tp = np.where(weight_sum > 0, weighted_theta / weight_sum, np.nan)
Tp[~ocean_mask] = np.nan

print(f"  valid: {np.count_nonzero(~np.isnan(Tp))}")
print(f"  range: {np.nanmin(Tp):.2f} ~ {np.nanmax(Tp):.2f} C")

# =============================================================
# [6] Tm: 500-1500m 깊이가중 평균
# =============================================================
print("\n[6] Tm (thetao 500-1500m)...")

idx_tm = np.where((deptht >= 500) & (deptht <= 1500))[0]
print(f"  levels: {idx_tm} → deptht = {np.round(deptht[idx_tm], 1)}")

if len(idx_tm) > 0:
    e3t_tm = e3t_0[idx_tm] * tmask[idx_tm]
    theta_tm = thetao_masked[idx_tm]
    weight_sum_tm = np.nansum(e3t_tm, axis=0)
    weighted_theta_tm = np.nansum(theta_tm * e3t_tm, axis=0)
    Tm = np.where(weight_sum_tm > 0, weighted_theta_tm / weight_sum_tm, np.nan)
else:
    Tm = np.full((ny, nx), np.nan)

# 얕은 해역 (500m 미만): Tm NaN → Tb로 대체
n_shallow = int((np.isnan(Tm) & ocean_mask).sum())
Tm = np.where(np.isnan(Tm) & ocean_mask, Tb, Tm)
Tm[~ocean_mask] = np.nan

print(f"  얕은 해역 Tb 대체: {n_shallow}개")
print(f"  range: {np.nanmin(Tm):.2f} ~ {np.nanmax(Tm):.2f} C")

ds_grid.close()

# =============================================================
# [7] Heup → photic (PISCES 유광층 깊이)
# =============================================================
print("\n[7] photic (Heup)...")

photic = ds_diad["Heup"].mean(dim=time_dim_diad).values
photic = np.where(np.abs(photic) > 1e10, np.nan, photic)
photic[~ocean_mask] = np.nan

print(f"  range: {np.nanmin(photic):.1f} ~ {np.nanmax(photic):.1f} m")

# =============================================================
# [8] EPC100 → dfbot (Martin curve)
# =============================================================
print("\n[8] dfbot (EPC100 → Martin curve)...")

epc100 = ds_diad["EPC100"].mean(dim=time_dim_diad).values  # mol C/m2/s
epc100 = np.where(np.abs(epc100) > 1e10, np.nan, epc100)

# mol C/m2/s → g ww/m2/yr at 100m
epc100_gww = epc100 * CONV
print(f"  EPC100 (mol/m2/s): {np.nanmin(epc100):.2e} ~ {np.nanmax(epc100):.2e}")
print(f"  EPC100 (gww/m2/yr): {np.nanmin(epc100_gww):.2f} ~ {np.nanmax(epc100_gww):.2f}")

# Martin curve: F(z_bot) = F(100m) × (z_bot/100)^(-b)
depth_martin = np.maximum(depth_bot, 100.0)
dfbot = epc100_gww * np.power(depth_martin / 100.0, -MARTIN_B)
print(f"  dfbot range: {np.nanmin(dfbot):.4f} ~ {np.nanmax(dfbot):.2f}")

# =============================================================
# [9] Closure term: mzrat × ZOO², mzrat2 × ZOO2²
#     Reference(COBALT)의 density-dependent closure term과 동일 개념
#     = 상위 포식자(물고기)에 의한 동물플랑크톤 에너지 손실
# =============================================================
print("\n[9] szprod/lzprod (closure term: mzrat × ZOO²)...")

ds_zoo = xr.open_dataset(zoo_file, decode_times=False)
time_dim_zoo = [d for d in ds_zoo["ZOO"].dims if "time" in d][0]

zoo_mean = ds_zoo["ZOO"].mean(dim=time_dim_zoo).values    # (z, y, x) mmol C/m3
zoo2_mean = ds_zoo["ZOO2"].mean(dim=time_dim_zoo).values

zoo_mean = np.where(np.abs(zoo_mean) > 1e10, np.nan, zoo_mean)
zoo2_mean = np.where(np.abs(zoo2_mean) > 1e10, np.nan, zoo2_mean)

print(f"  ZOO  shape: {zoo_mean.shape}")
print(f"  ZOO  range: {np.nanmin(zoo_mean):.4f} ~ {np.nanmax(zoo_mean):.4f} mmol/m3")
print(f"  ZOO2 range: {np.nanmin(zoo2_mean):.4f} ~ {np.nanmax(zoo2_mean):.4f} mmol/m3")

# closure term = mzrat × ZOO² (mmol C/m3/day)
closure_zoo = MZRAT * zoo_mean**2
closure_zoo2 = MZRAT2 * zoo2_mean**2

# 수직적분: sum(closure × e3t_0 × tmask) → mmol C/m2/day
closure_zoo_int = np.nansum(closure_zoo * e3t_0 * tmask, axis=0)
closure_zoo2_int = np.nansum(closure_zoo2 * e3t_0 * tmask, axis=0)

# mmol C/m2/day → g ww/m2/yr
szprod = closure_zoo_int * CONV_CLOSURE
lzprod_raw = closure_zoo2_int * CONV_CLOSURE

# lzprod 온도별 보정 (Reference: Stock et al. 2014, 2017)
# ratio = loss_to_fish / max_production → lzprod = closure / ratio
# 0°C→0.78, 10°C→0.85, 20°C→0.89, 30°C→0.91
T_KNOTS = np.array([0.0, 10.0, 20.0, 30.0])
R_KNOTS = np.array([0.78, 0.85, 0.89, 0.91])
ratio_map = np.interp(np.clip(Tp, 0, 30), T_KNOTS, R_KNOTS)
lzprod = lzprod_raw / ratio_map

szprod[~ocean_mask] = np.nan
lzprod[~ocean_mask] = np.nan

print(f"  szprod: {np.nanmin(szprod):.2f} ~ {np.nanmax(szprod):.2f} gww/m2/yr")
print(f"  lzprod (보정 전): {np.nanmin(lzprod_raw[ocean_mask]):.2f} ~ {np.nanmax(lzprod_raw[ocean_mask]):.2f}")
print(f"  lzprod (보정 후): {np.nanmin(lzprod):.2f} ~ {np.nanmax(lzprod):.2f} gww/m2/yr")

ds_zoo.close()
ds_diad.close()

# =============================================================
# [10] DataFrame 생성
# =============================================================
print("\n[10] DataFrame 생성...")

arrays = {
    "lon": lon2d, "lat": lat2d, "depth": depth_bot,
    "Tp": Tp, "Tm": Tm, "Tb": Tb,
    "photic": photic, "dfbot": dfbot,
    "szprod": szprod, "lzprod": lzprod,
}
print("  Shape 검증:")
for name, arr in arrays.items():
    print(f"    {name:8s}: {arr.shape}")
    assert arr.shape == (ny, nx), f"{name} shape mismatch!"

df = pd.DataFrame({
    "lon":    lon2d.ravel(),
    "lat":    lat2d.ravel(),
    "lzprod": lzprod.ravel(),
    "szprod": szprod.ravel(),
    "dfbot":  dfbot.ravel(),
    "photic": photic.ravel(),
    "depth":  depth_bot.ravel(),
    "Tb":     Tb.ravel(),
    "Tm":     Tm.ravel(),
    "Tp":     Tp.ravel(),
})

print(f"  전체: {len(df)}")

# 육지 제거
df = df.dropna(subset=["Tp", "depth"])
print(f"  바다: {len(df)}")

# 음수 보정
for col in ["lzprod", "szprod", "dfbot"]:
    n_neg = (df[col] < 0).sum()
    if n_neg > 0:
        print(f"  {col}: 음수 {n_neg}개 → 0")
        df[col] = df[col].clip(lower=0)

# 온도 하한 보정: 해수 어는점 = ~-1.8°C (35 PSU 기준)
T_MIN = -1.8
for col in ["Tb", "Tm", "Tp"]:
    n_below = (df[col] < T_MIN).sum()
    if n_below > 0:
        print(f"  {col}: {T_MIN}°C 미만 {n_below}개 → {T_MIN}°C로 클리핑")
        df[col] = df[col].clip(lower=T_MIN)

# NaN 잔여 확인
print("\n  NaN 개수:")
for col in df.columns:
    n_nan = df[col].isna().sum()
    if n_nan > 0:
        print(f"    {col}: {n_nan}")

# =============================================================
# [11] 저장
# =============================================================
out_file = f"{OUT_DIR}/Input_NEMO_FEISTY.csv"
df.to_csv(out_file, index=False)

print(f"\n[11] 저장: {out_file}")
print(f"  격자 수: {len(df)}")

# =============================================================
# [12] 요약 통계
# =============================================================
print(f"\n{'='*60}")
print("데이터 요약")
print(f"{'='*60}")
print(df.describe().round(3).to_string())

# =============================================================
# [13] 검증 플롯
# =============================================================
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 5, figsize=(25, 8))
    plot_info = [
        ("Tp",     "Tp [C]",             "RdYlBu_r"),
        ("Tm",     "Tm [C]",             "RdYlBu_r"),
        ("Tb",     "Tb [C]",             "RdYlBu_r"),
        ("depth",  "Depth [m]",          "viridis_r"),
        ("photic", "Euphotic [m]",       "YlGn_r"),
        ("lzprod", "lzprod [gww/m2/yr]", "YlOrRd"),
        ("szprod", "szprod [gww/m2/yr]", "YlOrRd"),
        ("dfbot",  "dfbot [gww/m2/yr]",  "YlOrBr"),
        ("lon",    "Longitude",          "twilight"),
        ("lat",    "Latitude",           "coolwarm"),
    ]

    for ax, (var, title, cmap) in zip(axes.ravel(), plot_info):
        vals = df[var].values
        finite = vals[np.isfinite(vals)]
        if len(finite) == 0:
            ax.set_title(f"{title} (NO DATA)")
            continue
        vmin, vmax = np.nanpercentile(finite, [2, 98])
        sc = ax.scatter(df["lon"], df["lat"], c=vals, s=0.3,
                        cmap=cmap, vmin=vmin, vmax=vmax,
                        edgecolors="none", rasterized=True)
        ax.set_title(title, fontsize=9)
        ax.set_xlim(-180, 180)
        ax.set_ylim(-90, 90)
        plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("NEMO-PISCES -> FEISTY Input (10yr mean, domcfg)", fontsize=13)
    plt.tight_layout()

    fig_path = f"{OUT_DIR}/Input_NEMO_FEISTY_check.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  검증 플롯: {fig_path}")

except ImportError as e:
    print(f"\n  matplotlib 없음 - 플롯 생략 ({e})")

print("\n" + "=" * 60)
print("완료!")
print("=" * 60)
