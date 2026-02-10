#!/usr/bin/env python3
"""
NEMO-PISCES → FEISTY Input CSV 변환 (1도 격자 버전)
CDO remapbil로 regrid된 파일 + domcfg 사용
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

grid_file = f"{SUBSET_DIR}/grid_T_10yr_1deg.nc"
diad_file = f"{SUBSET_DIR}/diad_T_1m_10yr_1deg.nc"
bl_file   = f"{SUBSET_DIR}/bottom_level_1deg.nc"
e3t_file  = f"{SUBSET_DIR}/e3t_0_1deg.nc"
# domcfg 원본은 01.org/ORCA_R2_zps_domcfg.nc → CDO가 regrid해서 위 파일들 생성

CONV = 365 * 86400 * 12.0 * 9.0
MARTIN_B = 0.858

print("=" * 60)
print("NEMO-PISCES → FEISTY Input (1deg regrid version)")
print(f"grid: {grid_file}")
print(f"diad: {diad_file}")
print(f"bl:   {bl_file}")
print(f"e3t:  {e3t_file}")
print("=" * 60)

# =============================================================
# [0] 파일 로드
# =============================================================
print("\n[0] 파일 로드...")
ds_grid = xr.open_dataset(grid_file, decode_times=False)
ds_diad = xr.open_dataset(diad_file, decode_times=False)
ds_bl   = xr.open_dataset(bl_file)
ds_e3t  = xr.open_dataset(e3t_file)

time_dim_grid = [d for d in ds_grid["thetao"].dims if "time" in d][0]
time_dim_diad = [d for d in ds_diad["Heup"].dims if "time" in d][0]

print(f"  grid time: {time_dim_grid}, steps: {len(ds_grid[time_dim_grid])}")
print(f"  diad time: {time_dim_diad}, steps: {len(ds_diad[time_dim_diad])}")

# 좌표: CDO regrid 후 lon/lat 1D
if "nav_lon" in ds_grid:
    lon1d = ds_grid["nav_lon"].values
    lat1d = ds_grid["nav_lat"].values
elif "lon" in ds_grid:
    lon1d = ds_grid["lon"].values
    lat1d = ds_grid["lat"].values
else:
    lon1d = ds_grid.coords[list(ds_grid.dims)[-1]].values
    lat1d = ds_grid.coords[list(ds_grid.dims)[-2]].values

# 1D → 2D meshgrid
if lon1d.ndim == 1:
    lon2d, lat2d = np.meshgrid(lon1d, lat1d)
else:
    lon2d, lat2d = lon1d, lat1d

# depth axis
depth_dim = [d for d in ds_grid["thetao"].dims if "depth" in d or d == "z"][0]
deptht = ds_grid[depth_dim].values
ny, nx = lon2d.shape
nz = len(deptht)

print(f"  격자: ({ny}, {nx}), 총 {ny*nx}")
print(f"  deptht ({nz} levels)")
print(f"  lon: {lon2d.min():.1f} ~ {lon2d.max():.1f}")
print(f"  lat: {lat2d.min():.1f} ~ {lat2d.max():.1f}")

# =============================================================
# [1] bottom_level, e3t_0, tmask
# =============================================================
print("\n[1] bottom_level, e3t_0 로드...")

bottom_level = ds_bl["bottom_level"].values
if bottom_level.ndim == 3:
    bottom_level = bottom_level[0]  # time dim 제거
bottom_level = np.round(bottom_level).astype(int)  # remapnn 후 float일 수 있음

e3t_0 = ds_e3t["e3t_0"].values
if e3t_0.ndim == 4:
    e3t_0 = e3t_0[0]  # time dim 제거

print(f"  bottom_level shape: {bottom_level.shape}, range: {bottom_level.min()} ~ {bottom_level.max()}")
print(f"  e3t_0 shape: {e3t_0.shape}")

ocean_mask = bottom_level > 0
n_ocean = ocean_mask.sum()
print(f"  바다 격자: {n_ocean}")

# tmask
k3d = np.arange(nz)[:, None, None]
tmask = (k3d < bottom_level[None, :, :]).astype(np.float64)

print(f"  tmask per level (처음 5): {[int(tmask[k].sum()) for k in range(min(5, nz))]}")

ds_bl.close()
ds_e3t.close()

# =============================================================
# [2] Bathymetry
# =============================================================
print("\n[2] Bathymetry 계산...")

cumsum_e3t = np.cumsum(e3t_0, axis=0)
yy, xx = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
bl_idx = np.clip(bottom_level - 1, 0, nz - 1)
depth_bot = cumsum_e3t[bl_idx, yy, xx]
depth_bot[~ocean_mask] = np.nan

print(f"  depth: {np.nanmin(depth_bot):.1f} ~ {np.nanmax(depth_bot):.1f} m, mean: {np.nanmean(depth_bot):.1f}")

# =============================================================
# [3] Tb
# =============================================================
print("\n[3] Tb (sbt)...")

Tb = ds_grid["sbt"].mean(dim=time_dim_grid).values
Tb = np.where(np.abs(Tb) > 1e10, np.nan, Tb)
Tb[~ocean_mask] = np.nan
print(f"  range: {np.nanmin(Tb):.2f} ~ {np.nanmax(Tb):.2f} C")

# =============================================================
# [4] thetao 시간평균 + tmask
# =============================================================
print("\n[4] thetao 시간평균...")

thetao_mean = ds_grid["thetao"].mean(dim=time_dim_grid).values
thetao_mean = np.where(np.abs(thetao_mean) > 1e10, np.nan, thetao_mean)
thetao_masked = np.where(tmask == 1, thetao_mean, np.nan)

# =============================================================
# [5] Tp: 0-100m
# =============================================================
print("\n[5] Tp (0-100m)...")

idx_tp = np.where(deptht <= 100)[0]
e3t_tp = e3t_0[idx_tp] * tmask[idx_tp]
theta_tp = thetao_masked[idx_tp]
weight_sum = np.nansum(e3t_tp, axis=0)
weighted_theta = np.nansum(theta_tp * e3t_tp, axis=0)
Tp = np.where(weight_sum > 0, weighted_theta / weight_sum, np.nan)
Tp[~ocean_mask] = np.nan
print(f"  range: {np.nanmin(Tp):.2f} ~ {np.nanmax(Tp):.2f} C")

# =============================================================
# [6] Tm: 500-1500m
# =============================================================
print("\n[6] Tm (500-1500m)...")

idx_tm = np.where((deptht >= 500) & (deptht <= 1500))[0]
if len(idx_tm) > 0:
    e3t_tm = e3t_0[idx_tm] * tmask[idx_tm]
    theta_tm = thetao_masked[idx_tm]
    weight_sum_tm = np.nansum(e3t_tm, axis=0)
    weighted_theta_tm = np.nansum(theta_tm * e3t_tm, axis=0)
    Tm = np.where(weight_sum_tm > 0, weighted_theta_tm / weight_sum_tm, np.nan)
else:
    Tm = np.full((ny, nx), np.nan)

n_shallow = int((np.isnan(Tm) & ocean_mask).sum())
Tm = np.where(np.isnan(Tm) & ocean_mask, Tb, Tm)
Tm[~ocean_mask] = np.nan
print(f"  얕은 해역 Tb 대체: {n_shallow}개")
print(f"  range: {np.nanmin(Tm):.2f} ~ {np.nanmax(Tm):.2f} C")

ds_grid.close()

# =============================================================
# [7] photic
# =============================================================
print("\n[7] photic (Heup)...")

photic = ds_diad["Heup"].mean(dim=time_dim_diad).values
photic = np.where(np.abs(photic) > 1e10, np.nan, photic)
photic[~ocean_mask] = np.nan
print(f"  range: {np.nanmin(photic):.1f} ~ {np.nanmax(photic):.1f} m")

# =============================================================
# [8] dfbot (Martin curve)
# =============================================================
print("\n[8] dfbot (EPC100 → Martin curve)...")

epc100 = ds_diad["EPC100"].mean(dim=time_dim_diad).values
epc100 = np.where(np.abs(epc100) > 1e10, np.nan, epc100)
epc100_gww = epc100 * CONV

depth_martin = np.maximum(depth_bot, 100.0)
dfbot = epc100_gww * np.power(depth_martin / 100.0, -MARTIN_B)
print(f"  dfbot: {np.nanmin(dfbot):.4f} ~ {np.nanmax(dfbot):.2f}")

# =============================================================
# [9] szprod/lzprod
# =============================================================
print("\n[9] szprod/lzprod (GRAZ1/GRAZ2)...")

graz1_mean = ds_diad["GRAZ1"].mean(dim=time_dim_diad).values
graz2_mean = ds_diad["GRAZ2"].mean(dim=time_dim_diad).values
graz1_mean = np.where(np.abs(graz1_mean) > 1e10, np.nan, graz1_mean)
graz2_mean = np.where(np.abs(graz2_mean) > 1e10, np.nan, graz2_mean)

graz1_int = np.nansum(graz1_mean * e3t_0 * tmask, axis=0)
graz2_int = np.nansum(graz2_mean * e3t_0 * tmask, axis=0)

szprod = graz1_int * CONV
lzprod = graz2_int * CONV
szprod[~ocean_mask] = np.nan
lzprod[~ocean_mask] = np.nan

print(f"  szprod: {np.nanmin(szprod):.2f} ~ {np.nanmax(szprod):.2f}")
print(f"  lzprod: {np.nanmin(lzprod):.2f} ~ {np.nanmax(lzprod):.2f}")

ds_diad.close()

# =============================================================
# [10] DataFrame
# =============================================================
print("\n[10] DataFrame 생성...")

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
df = df.dropna(subset=["Tp", "depth"])
print(f"  바다: {len(df)}")

for col in ["lzprod", "szprod", "dfbot"]:
    n_neg = (df[col] < 0).sum()
    if n_neg > 0:
        print(f"  {col}: 음수 {n_neg}개 → 0")
        df[col] = df[col].clip(lower=0)

print("\n  NaN 개수:")
for col in df.columns:
    n_nan = df[col].isna().sum()
    if n_nan > 0:
        print(f"    {col}: {n_nan}")

# =============================================================
# [11] 저장
# =============================================================
out_file = f"{OUT_DIR}/Input_NEMO_FEISTY_1deg.csv"
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

    fig.suptitle("NEMO-PISCES -> FEISTY Input (1deg regrid)", fontsize=13)
    plt.tight_layout()
    fig_path = f"{OUT_DIR}/Input_NEMO_FEISTY_1deg_check.png"
    plt.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  검증 플롯: {fig_path}")

except ImportError as e:
    print(f"\n  matplotlib 없음 ({e})")

print("\n" + "=" * 60)
print("완료!")
print("=" * 60)
