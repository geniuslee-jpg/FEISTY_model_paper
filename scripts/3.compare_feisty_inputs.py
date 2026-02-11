#!/usr/bin/env python3
"""
FEISTY Input 비교 검증 플롯
- 레퍼런스 (Input_global.csv) vs NEMO-PISCES (Input_NEMO_FEISTY.csv)
- 양쪽 모두 1°×1° (360×180) 정규 격자로 binning 후 Robinson 투영
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic_2d
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.util import add_cyclic_point

# =============================================================
# 파일 경로
# =============================================================
BASE_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY/03.feisty_input"
ref_file  = f"{BASE_DIR}/Input_global.csv"
nemo_file = f"{BASE_DIR}/Input_NEMO_FEISTY.csv"

print("Loading data...")
df_ref = pd.read_csv(ref_file)
df_nemo = pd.read_csv(nemo_file)

df_ref["lon"] = np.where(df_ref["lon"] > 180, df_ref["lon"] - 360, df_ref["lon"])
df_ref = df_ref.dropna(subset=["Tp", "depth"])

print(f"  Reference: {len(df_ref)} points")
print(f"  NEMO:      {len(df_nemo)} points")

# =============================================================
# 1° × 1° 정규 격자 정의
# =============================================================
lon_edges = np.arange(-180, 181, 1)  # 361 edges → 360 bins
lat_edges = np.arange(-90, 91, 1)    # 181 edges → 180 bins
lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2  # -179.5 ~ 179.5
lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2  # -89.5 ~ 89.5
LON, LAT = np.meshgrid(lon_centers, lat_centers)

def ref_to_grid(df, var):
    """Reference (~1° 정규격자) → 1°×1° binning."""
    vals = df[var].values
    mask = np.isfinite(vals)
    result = binned_statistic_2d(
        df["lon"].values[mask], df["lat"].values[mask], vals[mask],
        statistic="mean", bins=[lon_edges, lat_edges]
    )
    return result.statistic.T  # (nlat, nlon)

def nemo_to_grid(df, var):
    """NEMO ORCA2 (~2° 비정규격자) → 1°×1° nearest 보간 + 육지 마스킹."""
    vals = df[var].values
    mask = np.isfinite(vals)
    points = np.column_stack([df["lon"].values[mask], df["lat"].values[mask]])
    # nearest-neighbor 보간 (빈 셀 없이 채움)
    grid = griddata(points, vals[mask], (LON, LAT), method="nearest")
    # 육지 마스킹: 가장 가까운 데이터 포인트가 2.5° 이상이면 NaN
    tree = cKDTree(points)
    dist, _ = tree.query(np.column_stack([LON.ravel(), LAT.ravel()]))
    grid[dist.reshape(LON.shape) > 2.5] = np.nan
    return grid

print("  Regridding to 1°×1°...")

# =============================================================
# 투영 설정
# =============================================================
PROJ = ccrs.Robinson(central_longitude=0)
DATA_CRS = ccrs.PlateCarree()

def add_map_features(ax):
    """Robinson 투영 축에 해안선 + 육지 추가."""
    ax.set_global()
    ax.coastlines(linewidth=0.4, color="k")
    ax.add_feature(cfeature.LAND, facecolor="#d0d0d0", edgecolor="none", zorder=2)
    ax.gridlines(draw_labels=False, linewidth=0.2, color="grey",
                 alpha=0.5, linestyle="--")

def plot_on_robinson(ax, grid, cmap, vmin, vmax):
    """1°격자 데이터를 Robinson 투영에 pcolormesh."""
    # cyclic point 추가 (경도 180° 줄 방지)
    data_c, lon_c = add_cyclic_point(grid, coord=lon_centers)
    pc = ax.pcolormesh(lon_c, lat_centers, data_c,
                       transform=DATA_CRS,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="auto", rasterized=True)
    add_map_features(ax)
    return pc

# =============================================================
# 비교 변수 + 컬러바 범위
# =============================================================
variables = [
    ("Tp",     "Tp [\u00b0C]",           "RdYlBu_r", None, None),
    ("Tm",     "Tm [\u00b0C]",           "RdYlBu_r", None, None),
    ("Tb",     "Tb [\u00b0C]",           "RdYlBu_r", None, None),
    ("depth",  "Depth [m]",              "viridis_r", None, None),
    ("photic", "Euphotic depth [m]",     "YlGn_r",    0,    300),
    ("szprod", "szprod [gww/m\u00b2/yr]","YlOrRd",    None, None),
    ("lzprod", "lzprod [gww/m\u00b2/yr]","YlOrRd",    None, None),
    ("dfbot",  "dfbot [gww/m\u00b2/yr]", "YlOrBr",    0,    100),
]

nvar = len(variables)

# =============================================================
# [1] Side-by-side 맵 비교 (Robinson 투영)
# =============================================================
print("\n[1] Side-by-side map (Robinson projection)...")
fig, axes = plt.subplots(nvar, 2, figsize=(18, 4.5 * nvar),
                         subplot_kw={"projection": PROJ})

for i, (var, label, cmap, vmin_fix, vmax_fix) in enumerate(variables):
    ref_grid = ref_to_grid(df_ref, var)
    nemo_grid = nemo_to_grid(df_nemo, var)

    # 통계 (NaN 제외)
    ref_finite = ref_grid[np.isfinite(ref_grid)]
    nemo_finite = nemo_grid[np.isfinite(nemo_grid)]

    if vmin_fix is not None:
        vmin, vmax = vmin_fix, vmax_fix
    else:
        all_finite = np.concatenate([ref_finite, nemo_finite])
        vmin, vmax = np.nanpercentile(all_finite, [2, 98])

    # Reference
    ax = axes[i, 0]
    pc = plot_on_robinson(ax, ref_grid, cmap, vmin, vmax)
    ax.set_title(f"Reference: {label}", fontsize=10)
    ax.text(0.02, 0.02,
            f"n={len(ref_finite)}, mean={np.nanmean(ref_finite):.1f}, "
            f"min={np.nanmin(ref_finite):.1f}, max={np.nanmax(ref_finite):.1f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")

    # NEMO
    ax = axes[i, 1]
    pc = plot_on_robinson(ax, nemo_grid, cmap, vmin, vmax)
    ax.set_title(f"NEMO-PISCES: {label}", fontsize=10)
    ax.text(0.02, 0.02,
            f"n={len(nemo_finite)}, mean={np.nanmean(nemo_finite):.1f}, "
            f"min={np.nanmin(nemo_finite):.1f}, max={np.nanmax(nemo_finite):.1f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")

fig.suptitle("FEISTY Input: Reference vs NEMO-PISCES (Robinson)", fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(f"{BASE_DIR}/compare_maps.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {BASE_DIR}/compare_maps.png")

# =============================================================
# [2] 히스토그램 비교
# =============================================================
print("\n[2] Histogram...")
fig, axes = plt.subplots(2, 4, figsize=(20, 8))

for ax, (var, label, _, _, _) in zip(axes.ravel(), variables):
    rv = df_ref[var].dropna().values
    nv = df_nemo[var].dropna().values
    lo, hi = np.nanpercentile(np.concatenate([rv, nv]), [1, 99])
    bins = np.linspace(lo, hi, 60)

    ax.hist(rv, bins=bins, alpha=0.5, density=True, color="steelblue",
            label=f"Ref (n={len(rv)})")
    ax.hist(nv, bins=bins, alpha=0.5, density=True, color="tomato",
            label=f"NEMO (n={len(nv)})")
    ax.set_title(label, fontsize=10)
    ax.legend(fontsize=7)
    ax.set_ylabel("density")
    ax.tick_params(axis="both", which="major", length=6, width=0.8)
    ax.tick_params(axis="both", which="minor", length=3, width=0.5)

fig.suptitle("Distribution: Reference vs NEMO-PISCES", fontsize=14)
plt.tight_layout()
plt.savefig(f"{BASE_DIR}/compare_hist.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {BASE_DIR}/compare_hist.png")

# =============================================================
# [3] 변수별 개별 맵 (Robinson 투영)
# =============================================================
print("\n[3] 변수별 개별 맵...")
for var, label, cmap, vmin_fix, vmax_fix in variables:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5),
                                    subplot_kw={"projection": PROJ})

    ref_grid = ref_to_grid(df_ref, var)
    nemo_grid = nemo_to_grid(df_nemo, var)

    ref_f = ref_grid[np.isfinite(ref_grid)]
    nemo_f = nemo_grid[np.isfinite(nemo_grid)]

    if vmin_fix is not None:
        vmin, vmax = vmin_fix, vmax_fix
    else:
        all_f = np.concatenate([ref_f, nemo_f])
        vmin, vmax = np.nanpercentile(all_f, [2, 98])

    pc1 = plot_on_robinson(ax1, ref_grid, cmap, vmin, vmax)
    ax1.set_title(f"Reference: {label}", fontsize=12)
    plt.colorbar(pc1, ax=ax1, fraction=0.046, pad=0.04, orientation="horizontal")

    pc2 = plot_on_robinson(ax2, nemo_grid, cmap, vmin, vmax)
    ax2.set_title(f"NEMO-PISCES: {label}", fontsize=12)
    plt.colorbar(pc2, ax=ax2, fraction=0.046, pad=0.04, orientation="horizontal")

    plt.tight_layout()
    plt.savefig(f"{BASE_DIR}/compare_{var}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  {var} → compare_{var}.png")

# =============================================================
# [4] 요약 통계
# =============================================================
print(f"\n{'='*110}")
print(f"{'변수':>8s} | {'Ref mean':>10s} {'Ref std':>10s} {'Ref min':>10s} {'Ref max':>10s} | "
      f"{'NEMO mean':>10s} {'NEMO std':>10s} {'NEMO min':>10s} {'NEMO max':>10s} | {'비율':>6s}")
print("-" * 110)
for var, label, _, _, _ in variables:
    r = df_ref[var].dropna()
    n = df_nemo[var].dropna()
    ratio = n.mean() / r.mean() if r.mean() != 0 else float("inf")
    print(f"{var:>8s} | {r.mean():10.2f} {r.std():10.2f} {r.min():10.2f} {r.max():10.2f} | "
          f"{n.mean():10.2f} {n.std():10.2f} {n.min():10.2f} {n.max():10.2f} | {ratio:5.2f}x")
print(f"{'='*110}")

print(f"\n완료!")
