#!/usr/bin/env python3
"""
FEISTY 출력 비교 플롯
Reference (COBALT) vs NEMO-PISCES 어류 바이오매스 비교
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.interpolate import griddata
from scipy.spatial import cKDTree

# =============================================================
# 파일 경로
# =============================================================
BASE_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY/04.run"
OUT_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY/05.figures"

import os
os.makedirs(OUT_DIR, exist_ok=True)

ref_file = f"{BASE_DIR}/Ref_fish_biomass.csv"
nemo_file = f"{BASE_DIR}/NEMO_fish_biomass.csv"

print("Loading data...")
df_ref = pd.read_csv(ref_file)
df_nemo = pd.read_csv(nemo_file)

# Reference lon: 0~360 → -180~180
df_ref["lon"] = np.where(df_ref["lon"] > 180, df_ref["lon"] - 360, df_ref["lon"])

print(f"  Reference: {len(df_ref)} points")
print(f"  NEMO:      {len(df_nemo)} points")

# =============================================================
# 1° × 1° 정규 격자
# =============================================================
lon_edges = np.arange(-180, 181, 1)
lat_edges = np.arange(-90, 91, 1)
lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2
LON, LAT = np.meshgrid(lon_centers, lat_centers)

from scipy.stats import binned_statistic_2d

def ref_to_grid(df, var):
    vals = df[var].values
    mask = np.isfinite(vals)
    result = binned_statistic_2d(
        df["lon"].values[mask], df["lat"].values[mask], vals[mask],
        statistic="mean", bins=[lon_edges, lat_edges]
    )
    return result.statistic.T

def nemo_to_grid(df, var):
    vals = df[var].values
    mask = np.isfinite(vals)
    points = np.column_stack([df["lon"].values[mask], df["lat"].values[mask]])
    grid = griddata(points, vals[mask], (LON, LAT), method="nearest")
    tree = cKDTree(points)
    dist, _ = tree.query(np.column_stack([LON.ravel(), LAT.ravel()]))
    grid[dist.reshape(LON.shape) > 2.5] = np.nan
    return grid

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
# 변수 정의
# =============================================================
fish_groups = [
    ("totB_smpel",    "Small pelagic",    "YlOrRd"),
    ("totB_mesopel",  "Mesopelagic",      "YlGnBu"),
    ("totB_largepel", "Large pelagic",    "OrRd"),
    ("totB_midwpred", "Midwater predator","PuBuGn"),
    ("totB_dem",      "Demersal",         "YlOrBr"),
]

ngroup = len(fish_groups)

# Total biomass 계산
df_ref["totB_all"] = df_ref[[g[0] for g in fish_groups]].sum(axis=1)
df_nemo["totB_all"] = df_nemo[[g[0] for g in fish_groups]].sum(axis=1)

all_groups = fish_groups + [("totB_all", "Total biomass", "viridis")]

# =============================================================
# [1] Side-by-side 맵: 개별 그룹 + 전체
# =============================================================
print("\n[1] Side-by-side maps (6 groups)...")
fig, axes = plt.subplots(len(all_groups), 2, figsize=(20, 4 * len(all_groups)))

for i, (var, label, cmap) in enumerate(all_groups):
    ref_grid = ref_to_grid(df_ref, var)
    nemo_grid = nemo_to_grid(df_nemo, var)

    ref_f = ref_grid[np.isfinite(ref_grid)]
    nemo_f = nemo_grid[np.isfinite(nemo_grid)]
    all_f = np.concatenate([ref_f, nemo_f])
    vmin, vmax = 0, np.nanpercentile(all_f, 98)

    # Reference
    ax = axes[i, 0]
    add_land(ax)
    pc = ax.pcolormesh(lon_edges, lat_edges, ref_grid,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", rasterized=True, zorder=2)
    ax.set_title(f"Reference: {label}", fontsize=10)
    ax.text(0.02, 0.02,
            f"mean={np.nanmean(ref_f):.2f}, max={np.nanmax(ref_f):.1f}",
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
            f"mean={np.nanmean(nemo_f):.2f}, max={np.nanmax(nemo_f):.1f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=3)
    plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04)

fig.suptitle("FEISTY Fish Biomass: Reference vs NEMO-PISCES (1\u00b0 grid)",
             fontsize=14, y=1.01)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/feisty_compare_maps.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {OUT_DIR}/feisty_compare_maps.png")

# =============================================================
# [2] 변수별 개별 맵
# =============================================================
print("\n[2] 변수별 개별 맵...")
for var, label, cmap in all_groups:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 5))

    ref_grid = ref_to_grid(df_ref, var)
    nemo_grid = nemo_to_grid(df_nemo, var)

    ref_f = ref_grid[np.isfinite(ref_grid)]
    nemo_f = nemo_grid[np.isfinite(nemo_grid)]
    all_f = np.concatenate([ref_f, nemo_f])
    vmin, vmax = 0, np.nanpercentile(all_f, 98)

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
    plt.savefig(f"{OUT_DIR}/feisty_{var}.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  {var} → feisty_{var}.png")

# =============================================================
# [3] 히스토그램 비교
# =============================================================
print("\n[3] Histogram...")
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

for ax, (var, label, _) in zip(axes.ravel(), all_groups):
    rv = df_ref[var].values
    nv = df_nemo[var].values
    hi = np.nanpercentile(np.concatenate([rv, nv]), 99)
    bins = np.linspace(0, hi, 60)

    ax.hist(rv, bins=bins, alpha=0.5, density=True, color="steelblue",
            label=f"Ref (mean={np.mean(rv):.2f})")
    ax.hist(nv, bins=bins, alpha=0.5, density=True, color="tomato",
            label=f"NEMO (mean={np.mean(nv):.2f})")
    ax.set_title(label, fontsize=11)
    ax.legend(fontsize=8)
    ax.set_ylabel("density")
    ax.set_xlabel("Biomass [g/m\u00b2]")

fig.suptitle("Fish Biomass Distribution: Reference vs NEMO-PISCES", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/feisty_compare_hist.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {OUT_DIR}/feisty_compare_hist.png")

# =============================================================
# [4] 바 차트: 그룹별 평균 바이오매스 비교
# =============================================================
print("\n[4] Bar chart...")
fig, ax = plt.subplots(figsize=(10, 6))

labels = [g[1] for g in fish_groups]
ref_means = [df_ref[g[0]].mean() for g in fish_groups]
nemo_means = [df_nemo[g[0]].mean() for g in fish_groups]

x = np.arange(len(labels))
w = 0.35

bars1 = ax.bar(x - w/2, ref_means, w, label="Reference (COBALT)", color="steelblue", alpha=0.8)
bars2 = ax.bar(x + w/2, nemo_means, w, label="NEMO-PISCES", color="tomato", alpha=0.8)

# 비율 표시
for i, (r, n) in enumerate(zip(ref_means, nemo_means)):
    ratio = n / r if r > 0 else 0
    ax.text(i, max(r, n) + 0.1, f"{ratio:.2f}x", ha="center", fontsize=9, fontweight="bold")

ax.set_ylabel("Mean biomass [g/m\u00b2]", fontsize=11)
ax.set_title("Mean Fish Biomass by Functional Group", fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=10)
ax.legend(fontsize=10)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig(f"{OUT_DIR}/feisty_compare_bar.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {OUT_DIR}/feisty_compare_bar.png")

# =============================================================
# [5] 위도별 평균 바이오매스 (zonal mean)
# =============================================================
print("\n[5] Zonal mean...")
fig, axes = plt.subplots(2, 3, figsize=(18, 10))

for ax, (var, label, _) in zip(axes.ravel(), all_groups):
    ref_grid = ref_to_grid(df_ref, var)
    nemo_grid = nemo_to_grid(df_nemo, var)

    ref_zonal = np.nanmean(ref_grid, axis=1)
    nemo_zonal = np.nanmean(nemo_grid, axis=1)

    ax.plot(lat_centers, ref_zonal, color="steelblue", linewidth=1.5, label="Reference")
    ax.plot(lat_centers, nemo_zonal, color="tomato", linewidth=1.5, label="NEMO")
    ax.set_title(label, fontsize=11)
    ax.set_xlabel("Latitude", fontsize=9)
    ax.set_ylabel("Biomass [g/m\u00b2]", fontsize=9)
    ax.legend(fontsize=8)
    ax.set_xlim(-90, 90)
    ax.grid(alpha=0.3)

fig.suptitle("Zonal Mean Fish Biomass: Reference vs NEMO-PISCES", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/feisty_compare_zonal.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {OUT_DIR}/feisty_compare_zonal.png")

# =============================================================
# [6] NEMO 단독: 그룹별 비율 파이 차트
# =============================================================
print("\n[6] Pie chart (NEMO composition)...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

labels_short = ["Small pel.", "Mesopel.", "Large pel.", "Midw. pred.", "Demersal"]
colors = ["#e74c3c", "#3498db", "#e67e22", "#2ecc71", "#9b59b6"]

ref_means = [df_ref[g[0]].mean() for g in fish_groups]
nemo_means = [df_nemo[g[0]].mean() for g in fish_groups]

ax1.pie(ref_means, labels=labels_short, autopct="%1.1f%%", colors=colors, startangle=90)
ax1.set_title(f"Reference (total={sum(ref_means):.2f})", fontsize=11)

ax2.pie(nemo_means, labels=labels_short, autopct="%1.1f%%", colors=colors, startangle=90)
ax2.set_title(f"NEMO-PISCES (total={sum(nemo_means):.2f})", fontsize=11)

fig.suptitle("Fish Biomass Composition", fontsize=14)
plt.tight_layout()
plt.savefig(f"{OUT_DIR}/feisty_compare_pie.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {OUT_DIR}/feisty_compare_pie.png")

# =============================================================
# [7] 요약 통계
# =============================================================
print(f"\n{'='*90}")
print(f"{'그룹':>18s} | {'Ref mean':>10s} {'Ref max':>10s} | {'NEMO mean':>10s} {'NEMO max':>10s} | {'비율':>6s}")
print("-" * 90)
for var, label, _ in all_groups:
    r = df_ref[var]
    n = df_nemo[var]
    ratio = n.mean() / r.mean() if r.mean() > 0 else 0
    print(f"{label:>18s} | {r.mean():10.3f} {r.max():10.3f} | {n.mean():10.3f} {n.max():10.3f} | {ratio:5.2f}x")
print(f"{'='*90}")

print(f"\n완료! 그림 {OUT_DIR}/ 에 저장")
