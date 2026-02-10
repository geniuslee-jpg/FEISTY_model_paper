#!/usr/bin/env python3
"""
FEISTY Input 비교 검증 플롯
- 레퍼런스 (Input_global.csv) vs NEMO-PISCES (Input_NEMO_FEISTY.csv)
- 양쪽 모두 1°×1° (360×180) 정규 격자로 binning 후 pcolormesh
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import binned_statistic_2d

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

def scatter_to_grid(df, var):
    """산점 데이터를 1°×1° 격자로 binning (평균)."""
    vals = df[var].values
    mask = np.isfinite(vals)
    result = binned_statistic_2d(
        df["lon"].values[mask], df["lat"].values[mask], vals[mask],
        statistic="mean", bins=[lon_edges, lat_edges]
    )
    return result.statistic.T  # (nlat, nlon)

print("  Regridding to 1°×1°...")

# =============================================================
# 축 스타일
# =============================================================
def add_land(ax):
    ax.set_facecolor("#d0d0d0")
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude", fontsize=9)
    ax.set_ylabel("Latitude", fontsize=9)
    ax.tick_params(axis="both", which="major", length=6, width=0.8, labelsize=8)
    ax.tick_params(axis="both", which="minor", length=3, width=0.5)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(60))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(30))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(30))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(15))

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
# [1] Side-by-side 맵 비교 (pcolormesh)
# =============================================================
print("\n[1] Side-by-side map (1° gridded)...")
fig, axes = plt.subplots(nvar, 2, figsize=(20, 4 * nvar))

for i, (var, label, cmap, vmin_fix, vmax_fix) in enumerate(variables):
    ref_grid = scatter_to_grid(df_ref, var)
    nemo_grid = scatter_to_grid(df_nemo, var)

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
    add_land(ax)
    pc = ax.pcolormesh(lon_edges, lat_edges, ref_grid,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", rasterized=True, zorder=2)
    ax.set_title(f"Reference: {label}", fontsize=10)
    ax.text(0.02, 0.02,
            f"n={len(ref_finite)}, mean={np.nanmean(ref_finite):.1f}, "
            f"min={np.nanmin(ref_finite):.1f}, max={np.nanmax(ref_finite):.1f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=3)
    plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)

    # NEMO
    ax = axes[i, 1]
    add_land(ax)
    pc = ax.pcolormesh(lon_edges, lat_edges, nemo_grid,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", rasterized=True, zorder=2)
    ax.set_title(f"NEMO-PISCES: {label}", fontsize=10)
    ax.text(0.02, 0.02,
            f"n={len(nemo_finite)}, mean={np.nanmean(nemo_finite):.1f}, "
            f"min={np.nanmin(nemo_finite):.1f}, max={np.nanmax(nemo_finite):.1f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=3)
    plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)

fig.suptitle("FEISTY Input: Reference vs NEMO-PISCES (1\u00b0 grid)", fontsize=14, y=1.01)
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
# [3] 변수별 개별 맵
# =============================================================
print("\n[3] 변수별 개별 맵...")
for var, label, cmap, vmin_fix, vmax_fix in variables:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5))

    ref_grid = scatter_to_grid(df_ref, var)
    nemo_grid = scatter_to_grid(df_nemo, var)

    ref_f = ref_grid[np.isfinite(ref_grid)]
    nemo_f = nemo_grid[np.isfinite(nemo_grid)]

    if vmin_fix is not None:
        vmin, vmax = vmin_fix, vmax_fix
    else:
        all_f = np.concatenate([ref_f, nemo_f])
        vmin, vmax = np.nanpercentile(all_f, [2, 98])

    add_land(ax1)
    pc1 = ax1.pcolormesh(lon_edges, lat_edges, ref_grid,
                         cmap=cmap, vmin=vmin, vmax=vmax,
                         shading="flat", rasterized=True, zorder=2)
    ax1.set_title(f"Reference: {label}", fontsize=12)
    plt.colorbar(pc1, ax=ax1, fraction=0.046, pad=0.04)

    add_land(ax2)
    pc2 = ax2.pcolormesh(lon_edges, lat_edges, nemo_grid,
                         cmap=cmap, vmin=vmin, vmax=vmax,
                         shading="flat", rasterized=True, zorder=2)
    ax2.set_title(f"NEMO-PISCES: {label}", fontsize=12)
    plt.colorbar(pc2, ax=ax2, fraction=0.046, pad=0.04)

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
