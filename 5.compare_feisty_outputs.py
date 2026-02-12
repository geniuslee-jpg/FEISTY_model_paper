#!/usr/bin/env python3
"""
FEISTY Output 비교 검증 플롯
- 레퍼런스 (Global_fish_biomass.RData) vs NEMO-PISCES (Global_fish_biomass.RData)
- Mollweide 투영 + 1°×1° regridding (3.compare_feisty_inputs.py 스타일)
- 필요 패키지: pyreadr, cartopy, scipy
"""

import numpy as np
import pandas as pd
import pyreadr
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
# 파일 경로 (수정 필요)
# =============================================================
ref_file  = "/path/to/reference/Global_fish_biomass.RData"       # ← 레퍼런스 경로 입력
nemo_file = "/data01/labdisk/sungjin/NEMO-FEISTY/04.run/data/Global_fish_biomass.RData"
out_dir   = "/data01/labdisk/sungjin/NEMO-FEISTY/05.Figure"

# =============================================================
# 데이터 로드
# =============================================================
print("Loading data...")
df_ref  = pyreadr.read_r(ref_file)["out"]
df_nemo = pyreadr.read_r(nemo_file)["out"]

# lon 보정: 0~360 → -180~180
df_ref["lon"]  = np.where(df_ref["lon"] > 180, df_ref["lon"] - 360, df_ref["lon"])
df_nemo["lon"] = np.where(df_nemo["lon"] > 180, df_nemo["lon"] - 360, df_nemo["lon"])

# Total biomass 계산
biomass_cols = ["totB_smpel", "totB_mesopel", "totB_largepel", "totB_midwpred", "totB_dem"]
df_ref["totB"]  = df_ref[biomass_cols].sum(axis=1)
df_nemo["totB"] = df_nemo[biomass_cols].sum(axis=1)

print(f"  Reference: {len(df_ref)} points")
print(f"  NEMO:      {len(df_nemo)} points")

# =============================================================
# 1° × 1° 정규 격자 정의
# =============================================================
lon_edges = np.arange(-180, 181, 1)
lat_edges = np.arange(-90, 91, 1)
lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2
LON, LAT = np.meshgrid(lon_centers, lat_centers)

def ref_to_grid(df, var):
    """Reference (~1° 정규격자) → 1°×1° binning."""
    vals = df[var].values
    mask = np.isfinite(vals)
    result = binned_statistic_2d(
        df["lon"].values[mask], df["lat"].values[mask], vals[mask],
        statistic="mean", bins=[lon_edges, lat_edges]
    )
    return result.statistic.T

def nemo_to_grid(df, var):
    """NEMO ORCA2 (~2° 비정규격자) → 1°×1° nearest 보간 + 육지 마스킹."""
    vals = df[var].values
    mask = np.isfinite(vals)
    points = np.column_stack([df["lon"].values[mask], df["lat"].values[mask]])
    grid = griddata(points, vals[mask], (LON, LAT), method="nearest")
    tree = cKDTree(points)
    dist, _ = tree.query(np.column_stack([LON.ravel(), LAT.ravel()]))
    grid[dist.reshape(LON.shape) > 2.5] = np.nan
    return grid

print("  Regridding to 1°×1°...")

# =============================================================
# 투영 설정
# =============================================================
PROJ = ccrs.Mollweide(central_longitude=0)
DATA_CRS = ccrs.PlateCarree()

def add_map_features(ax):
    ax.set_global()
    ax.coastlines(linewidth=0.4, color="k")
    ax.add_feature(cfeature.LAND, facecolor="#d0d0d0", edgecolor="none", zorder=2)
    ax.gridlines(draw_labels=False, linewidth=0.2, color="grey",
                 alpha=0.5, linestyle="--")

def plot_on_mollweide(ax, grid, cmap, vmin, vmax):
    data_c, lon_c = add_cyclic_point(grid, coord=lon_centers)
    pc = ax.pcolormesh(lon_c, lat_centers, data_c,
                       transform=DATA_CRS,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="auto", rasterized=True)
    add_map_features(ax)
    return pc

# =============================================================
# 비교 변수 설정 (log10 스케일, 고정 컬러바 범위)
# =============================================================
variables = [
    ("totB_smpel",    "Small pelagic [g/m\u00b2]",     "YlOrRd", -3, 2),
    ("totB_mesopel",  "Mesopelagic [g/m\u00b2]",       "YlOrRd", -3, 2),
    ("totB_largepel", "Large pelagic [g/m\u00b2]",     "YlOrRd", -3, 2),
    ("totB_midwpred", "Midwater predator [g/m\u00b2]", "YlOrRd", -3, 2),
    ("totB_dem",      "Demersal [g/m\u00b2]",          "YlOrRd", -3, 2),
    ("totB",          "Total biomass [g/m\u00b2]",     "viridis", -2, 2),
]

nvar = len(variables)

# =============================================================
# [1] Side-by-side 맵 비교 (Mollweide 투영, log10 스케일)
# =============================================================
print("\n[1] Side-by-side map (Mollweide, log10 scale)...")
fig, axes = plt.subplots(nvar, 2, figsize=(18, 4.5 * nvar),
                         subplot_kw={"projection": PROJ})

for i, (var, label, cmap, vmin_fix, vmax_fix) in enumerate(variables):
    ref_grid  = np.log10(np.maximum(ref_to_grid(df_ref, var), 1e-4))
    nemo_grid = np.log10(np.maximum(nemo_to_grid(df_nemo, var), 1e-4))

    vmin, vmax = vmin_fix, vmax_fix

    # Reference
    ax = axes[i, 0]
    pc = plot_on_mollweide(ax, ref_grid, cmap, vmin, vmax)
    ax.set_title(f"Reference: {label}", fontsize=10)
    ref_finite = ref_grid[np.isfinite(ref_grid)]
    ax.text(0.02, 0.02,
            f"mean={np.nanmean(ref_finite):.2f}, "
            f"min={np.nanmin(ref_finite):.2f}, max={np.nanmax(ref_finite):.2f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label("log10(g/m\u00b2)", fontsize=8)

    # NEMO
    ax = axes[i, 1]
    pc = plot_on_mollweide(ax, nemo_grid, cmap, vmin, vmax)
    ax.set_title(f"NEMO-PISCES: {label}", fontsize=10)
    nemo_finite = nemo_grid[np.isfinite(nemo_grid)]
    ax.text(0.02, 0.02,
            f"mean={np.nanmean(nemo_finite):.2f}, "
            f"min={np.nanmin(nemo_finite):.2f}, max={np.nanmax(nemo_finite):.2f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label("log10(g/m\u00b2)", fontsize=8)

fig.suptitle("FEISTY Output: Reference vs NEMO-PISCES (Mollweide)", fontsize=14, y=1.01)
plt.tight_layout()
path1 = f"{out_dir}/compare_output_maps.png"
plt.savefig(path1, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path1}")

# =============================================================
# [2] 히스토그램 비교 (log10 스케일)
# =============================================================
print("\n[2] Histogram...")
fig, axes = plt.subplots(2, 3, figsize=(18, 9))

for ax, (var, label, *_) in zip(axes.ravel(), variables):
    rv = df_ref[var].dropna().values
    nv = df_nemo[var].dropna().values
    rv = np.log10(rv[rv > 0])
    nv = np.log10(nv[nv > 0])

    lo, hi = np.nanpercentile(np.concatenate([rv, nv]), [1, 99])
    bins = np.linspace(lo, hi, 60)

    ax.hist(rv, bins=bins, alpha=0.5, density=True, color="steelblue",
            label=f"Ref (n={len(rv)})")
    ax.hist(nv, bins=bins, alpha=0.5, density=True, color="tomato",
            label=f"NEMO (n={len(nv)})")
    ax.set_title(label, fontsize=10)
    ax.set_xlabel("log10(g/m\u00b2)")
    ax.set_ylabel("density")
    ax.legend(fontsize=7)
    ax.tick_params(axis="both", which="major", length=6, width=0.8)

fig.suptitle("Biomass Distribution: Reference vs NEMO-PISCES", fontsize=14)
plt.tight_layout()
path2 = f"{out_dir}/compare_output_hist.png"
plt.savefig(path2, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path2}")

# =============================================================
# [3] 변수별 개별 맵 + Bias (Mollweide 투영, log10 스케일)
# =============================================================
print("\n[3] 변수별 개별 맵 + Bias...")
for var, label, cmap, vmin_fix, vmax_fix in variables:
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 5),
                                         subplot_kw={"projection": PROJ})

    ref_grid  = np.log10(np.maximum(ref_to_grid(df_ref, var), 1e-4))
    nemo_grid = np.log10(np.maximum(nemo_to_grid(df_nemo, var), 1e-4))

    vmin, vmax = vmin_fix, vmax_fix

    # Reference
    pc1 = plot_on_mollweide(ax1, ref_grid, cmap, vmin, vmax)
    ax1.set_title(f"Reference: {label}", fontsize=12)
    cb1 = plt.colorbar(pc1, ax=ax1, fraction=0.046, pad=0.04, orientation="horizontal")
    cb1.set_label("log10(g/m\u00b2)", fontsize=8)

    # NEMO
    pc2 = plot_on_mollweide(ax2, nemo_grid, cmap, vmin, vmax)
    ax2.set_title(f"NEMO-PISCES: {label}", fontsize=12)
    cb2 = plt.colorbar(pc2, ax=ax2, fraction=0.046, pad=0.04, orientation="horizontal")
    cb2.set_label("log10(g/m\u00b2)", fontsize=8)

    # Bias (NEMO - Reference)
    bias_grid = nemo_grid - ref_grid
    mask = np.isfinite(ref_grid) & np.isfinite(nemo_grid)
    bias_grid[~mask] = np.nan
    bias_f = bias_grid[np.isfinite(bias_grid)]
    blim = np.nanpercentile(np.abs(bias_f), 98)

    pc3 = plot_on_mollweide(ax3, bias_grid, "RdBu_r", -blim, blim)
    ax3.set_title(f"Bias (NEMO - Ref): {label}", fontsize=12)
    ax3.text(0.02, 0.02,
             f"mean={np.nanmean(bias_f):.2f}, "
             f"std={np.nanstd(bias_f):.2f}",
             transform=ax3.transAxes, fontsize=7, va="bottom",
             bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    cb3 = plt.colorbar(pc3, ax=ax3, fraction=0.046, pad=0.04, orientation="horizontal")
    cb3.set_label("\u0394 log10(g/m\u00b2)", fontsize=8)

    plt.tight_layout()
    plt.savefig(f"{out_dir}/compare_output_{var}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  {var} → compare_output_{var}.png")

# =============================================================
# [4] 기능그룹별 비율 비교 (bar chart)
# =============================================================
print("\n[4] 기능그룹 비율 비교...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

group_labels = ["Small pel.", "Mesopelagic", "Large pel.", "Midwater pred.", "Demersal"]
colors = ["#FFD700", "#FFA500", "#FF4500", "#8B0000", "#4682B4"]

ref_means  = [df_ref[c].mean()  for c in biomass_cols]
nemo_means = [df_nemo[c].mean() for c in biomass_cols]

ax = axes[0]
bars = ax.barh(group_labels, ref_means, color=colors, edgecolor="black", linewidth=0.5)
ax.set_xlabel("Mean biomass [g/m\u00b2]")
ax.set_title("Reference")
for bar, val in zip(bars, ref_means):
    ax.text(bar.get_width() + 0.001, bar.get_y() + bar.get_height()/2,
            f"{val:.3f}", va="center", fontsize=8)

ax = axes[1]
bars = ax.barh(group_labels, nemo_means, color=colors, edgecolor="black", linewidth=0.5)
ax.set_xlabel("Mean biomass [g/m\u00b2]")
ax.set_title("NEMO-PISCES")
for bar, val in zip(bars, nemo_means):
    ax.text(bar.get_width() + 0.001, bar.get_y() + bar.get_height()/2,
            f"{val:.3f}", va="center", fontsize=8)

xmax = max(max(ref_means), max(nemo_means)) * 1.3
axes[0].set_xlim(0, xmax)
axes[1].set_xlim(0, xmax)

fig.suptitle("Functional Group Mean Biomass Comparison", fontsize=14)
plt.tight_layout()
path4 = f"{out_dir}/compare_output_groups.png"
plt.savefig(path4, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path4}")

# =============================================================
# [5] 요약 통계
# =============================================================
print(f"\n{'='*130}")
print(f"{'변수':>15s} | {'Ref mean':>10s} {'Ref std':>10s} {'Ref min':>10s} {'Ref max':>10s} | "
      f"{'NEMO mean':>10s} {'NEMO std':>10s} {'NEMO min':>10s} {'NEMO max':>10s} | {'ratio':>6s}")
print("-" * 130)
for var, label, *_ in variables:
    r = df_ref[var].dropna()
    n = df_nemo[var].dropna()
    ratio = n.mean() / r.mean() if r.mean() != 0 else float("inf")
    print(f"{var:>15s} | {r.mean():10.4f} {r.std():10.4f} {r.min():10.4f} {r.max():10.4f} | "
          f"{n.mean():10.4f} {n.std():10.4f} {n.min():10.4f} {n.max():10.4f} | {ratio:5.2f}x")
print(f"{'='*130}")

print(f"\n완료!")
print(f"  맵 전체:      {path1}")
print(f"  히스토그램:    {path2}")
print(f"  그룹 비교:    {path4}")
print(f"  개별 맵:      {out_dir}/compare_output_<var>.png")
