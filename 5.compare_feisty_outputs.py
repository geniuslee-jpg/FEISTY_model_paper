#!/usr/bin/env python3
"""
FEISTY Output 비교 검증 플롯
- 레퍼런스 (Global_fish_biomass.RData) vs NEMO-PISCES (Global_fish_biomass.RData)
- RData 파일을 읽기 위해 pyreadr 필요: pip install pyreadr
"""

import numpy as np
import pandas as pd
import pyreadr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
# 비교 변수 설정
# =============================================================
variables = [
    ("totB_smpel",    "Small pelagic [g/m2]",     "YlOrRd"),
    ("totB_mesopel",  "Mesopelagic [g/m2]",       "YlOrRd"),
    ("totB_largepel", "Large pelagic [g/m2]",     "YlOrRd"),
    ("totB_midwpred", "Midwater predator [g/m2]", "YlOrRd"),
    ("totB_dem",      "Demersal [g/m2]",          "YlOrRd"),
    ("totB",          "Total biomass [g/m2]",     "viridis"),
]

nvar = len(variables)

# =============================================================
# 플롯 1: 나란히 비교 맵 (log scale)
# =============================================================
print("\n[1] Side-by-side map 생성 (log scale)...")
fig, axes = plt.subplots(nvar, 2, figsize=(20, 4 * nvar))

for i, (var, label, cmap) in enumerate(variables):
    ref_vals  = df_ref[var].values
    nemo_vals = df_nemo[var].values

    # log10 변환 (0 이하는 작은 값으로 대체)
    ref_log  = np.log10(np.maximum(ref_vals, 1e-4))
    nemo_log = np.log10(np.maximum(nemo_vals, 1e-4))

    # 공통 색상 범위
    all_log = np.concatenate([ref_log[np.isfinite(ref_log)],
                              nemo_log[np.isfinite(nemo_log)]])
    vmin, vmax = np.nanpercentile(all_log, [2, 98])

    # Reference
    ax_ref = axes[i, 0]
    sc = ax_ref.scatter(df_ref["lon"], df_ref["lat"], c=ref_log, s=0.3,
                        cmap=cmap, vmin=vmin, vmax=vmax,
                        edgecolors="none", rasterized=True)
    ax_ref.set_xlim(-180, 180)
    ax_ref.set_ylim(-90, 90)
    ax_ref.set_title(f"Reference: {label}", fontsize=10)
    cb = plt.colorbar(sc, ax=ax_ref, fraction=0.046, pad=0.04)
    cb.set_label("log10(g/m2)", fontsize=7)

    # NEMO
    ax_nemo = axes[i, 1]
    sc2 = ax_nemo.scatter(df_nemo["lon"], df_nemo["lat"], c=nemo_log, s=0.5,
                          cmap=cmap, vmin=vmin, vmax=vmax,
                          edgecolors="none", rasterized=True)
    ax_nemo.set_xlim(-180, 180)
    ax_nemo.set_ylim(-90, 90)
    ax_nemo.set_title(f"NEMO-PISCES: {label}", fontsize=10)
    cb2 = plt.colorbar(sc2, ax=ax_nemo, fraction=0.046, pad=0.04)
    cb2.set_label("log10(g/m2)", fontsize=7)

    # 통계 표시 (원래 스케일)
    ref_finite  = ref_vals[np.isfinite(ref_vals) & (ref_vals > 0)]
    nemo_finite = nemo_vals[np.isfinite(nemo_vals) & (nemo_vals > 0)]
    ax_ref.text(0.02, 0.02,
                f"n={len(ref_finite)}, mean={np.nanmean(ref_finite):.3f}, "
                f"median={np.nanmedian(ref_finite):.3f}, max={np.nanmax(ref_finite):.2f}",
                transform=ax_ref.transAxes, fontsize=7, va="bottom",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax_nemo.text(0.02, 0.02,
                 f"n={len(nemo_finite)}, mean={np.nanmean(nemo_finite):.3f}, "
                 f"median={np.nanmedian(nemo_finite):.3f}, max={np.nanmax(nemo_finite):.2f}",
                 transform=ax_nemo.transAxes, fontsize=7, va="bottom",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

fig.suptitle("FEISTY Output Comparison: Reference vs NEMO-PISCES", fontsize=14, y=1.01)
plt.tight_layout()

path1 = f"{out_dir}/compare_output_maps.png"
plt.savefig(path1, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path1}")

# =============================================================
# 플롯 2: 히스토그램 비교 (log scale)
# =============================================================
print("\n[2] Histogram 비교 생성...")
fig, axes = plt.subplots(2, 3, figsize=(18, 9))

for ax, (var, label, _) in zip(axes.ravel(), variables):
    ref_vals  = df_ref[var].dropna().values
    nemo_vals = df_nemo[var].dropna().values

    # log10 (양수만)
    ref_pos  = ref_vals[ref_vals > 0]
    nemo_pos = nemo_vals[nemo_vals > 0]
    ref_log  = np.log10(ref_pos)
    nemo_log = np.log10(nemo_pos)

    all_log = np.concatenate([ref_log, nemo_log])
    lo, hi = np.nanpercentile(all_log, [1, 99])
    bins = np.linspace(lo, hi, 60)

    ax.hist(ref_log, bins=bins, alpha=0.5, label=f"Ref (n={len(ref_pos)})",
            density=True, color="steelblue")
    ax.hist(nemo_log, bins=bins, alpha=0.5, label=f"NEMO (n={len(nemo_pos)})",
            density=True, color="tomato")
    ax.set_title(label, fontsize=10)
    ax.set_xlabel("log10(g/m2)")
    ax.set_ylabel("density")
    ax.legend(fontsize=7)

fig.suptitle("Biomass Distribution: Reference vs NEMO-PISCES", fontsize=14)
plt.tight_layout()

path2 = f"{out_dir}/compare_output_hist.png"
plt.savefig(path2, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path2}")

# =============================================================
# 플롯 3: 기능그룹별 비율 비교 (파이 차트 스타일 bar)
# =============================================================
print("\n[3] 기능그룹 비율 비교...")
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

group_labels = ["Small pel.", "Mesopelagic", "Large pel.", "Midwater pred.", "Demersal"]
colors = ["#FFD700", "#FFA500", "#FF4500", "#8B0000", "#4682B4"]

ref_means  = [df_ref[c].mean()  for c in biomass_cols]
nemo_means = [df_nemo[c].mean() for c in biomass_cols]

# Reference
ax = axes[0]
bars = ax.barh(group_labels, ref_means, color=colors, edgecolor="black", linewidth=0.5)
ax.set_xlabel("Mean biomass [g/m2]")
ax.set_title("Reference")
for bar, val in zip(bars, ref_means):
    ax.text(bar.get_width() + 0.001, bar.get_y() + bar.get_height()/2,
            f"{val:.3f}", va="center", fontsize=8)

# NEMO
ax = axes[1]
bars = ax.barh(group_labels, nemo_means, color=colors, edgecolor="black", linewidth=0.5)
ax.set_xlabel("Mean biomass [g/m2]")
ax.set_title("NEMO-PISCES")
for bar, val in zip(bars, nemo_means):
    ax.text(bar.get_width() + 0.001, bar.get_y() + bar.get_height()/2,
            f"{val:.3f}", va="center", fontsize=8)

# 동일한 x축 범위
xmax = max(max(ref_means), max(nemo_means)) * 1.3
axes[0].set_xlim(0, xmax)
axes[1].set_xlim(0, xmax)

fig.suptitle("Functional Group Mean Biomass Comparison", fontsize=14)
plt.tight_layout()

path3 = f"{out_dir}/compare_output_groups.png"
plt.savefig(path3, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path3}")

# =============================================================
# 요약 통계 비교 테이블
# =============================================================
print("\n[4] 요약 통계 비교")
all_vars = variables
print(f"\n{'변수':>15s} | {'Ref mean':>10s} {'Ref med':>10s} {'Ref max':>10s} | {'NEMO mean':>10s} {'NEMO med':>10s} {'NEMO max':>10s} | {'ratio':>6s}")
print("-" * 100)
for var, label, _ in all_vars:
    r = df_ref[var].dropna()
    n = df_nemo[var].dropna()
    ratio = n.mean() / r.mean() if r.mean() != 0 else float("inf")
    print(f"{var:>15s} | {r.mean():10.4f} {r.median():10.4f} {r.max():10.4f} | "
          f"{n.mean():10.4f} {n.median():10.4f} {n.max():10.4f} | {ratio:6.2f}x")

print(f"\n완료!")
print(f"  맵 비교:      {path1}")
print(f"  히스토그램:    {path2}")
print(f"  그룹별 비교:  {path3}")
