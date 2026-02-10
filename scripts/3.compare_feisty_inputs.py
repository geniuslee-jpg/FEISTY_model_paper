#!/usr/bin/env python3
"""
FEISTY Input 비교 검증 플롯
- 레퍼런스 (Input_global.csv) vs NEMO-PISCES (Input_NEMO_FEISTY.csv)
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================
# 파일 경로
# =============================================================
ref_file = "/home/user/FEISTY_model_paper/quarto/data/Input_global.csv"
nemo_file = "/data01/labdisk/sungjin/NEMO-FEISTY/03.feisty_input/Input_NEMO_FEISTY.csv"
out_dir = "/data01/labdisk/sungjin/NEMO-FEISTY/03.feisty_input"

print("Loading data...")
df_ref = pd.read_csv(ref_file)
df_nemo = pd.read_csv(nemo_file)

# 레퍼런스 lon: 0~360 → -180~180 변환
df_ref["lon"] = np.where(df_ref["lon"] > 180, df_ref["lon"] - 360, df_ref["lon"])

# NA 제거
df_ref = df_ref.dropna(subset=["Tp", "depth"])

print(f"  Reference: {len(df_ref)} points")
print(f"  NEMO:      {len(df_nemo)} points")

# =============================================================
# 비교 변수 설정
# =============================================================
variables = [
    ("Tp",     "Tp [C]",              "RdYlBu_r", None,       None),
    ("Tm",     "Tm [C]",              "RdYlBu_r", None,       None),
    ("Tb",     "Tb [C]",              "RdYlBu_r", None,       None),
    ("depth",  "Depth [m]",           "viridis_r", None,       None),
    ("photic", "Euphotic depth [m]",  "YlGn_r",    None,       None),
    ("szprod", "szprod [gww/m2/yr]",  "YlOrRd",    None,       None),
    ("lzprod", "lzprod [gww/m2/yr]",  "YlOrRd",    None,       None),
    ("dfbot",  "dfbot [gww/m2/yr]",   "YlOrBr",    None,       None),
]

nvar = len(variables)

# =============================================================
# 플롯 1: 나란히 비교 (Reference vs NEMO)
# =============================================================
print("\n[1] Side-by-side map 생성...")
fig, axes = plt.subplots(nvar, 2, figsize=(20, 4 * nvar))

for i, (var, label, cmap, vmin_fix, vmax_fix) in enumerate(variables):
    # 레퍼런스 데이터
    ref_vals = df_ref[var].values
    ref_finite = ref_vals[np.isfinite(ref_vals)]

    # NEMO 데이터
    nemo_vals = df_nemo[var].values
    nemo_finite = nemo_vals[np.isfinite(nemo_vals)]

    # 공통 색상 범위 (2~98 percentile)
    all_finite = np.concatenate([ref_finite, nemo_finite])
    if vmin_fix is not None:
        vmin, vmax = vmin_fix, vmax_fix
    else:
        vmin, vmax = np.nanpercentile(all_finite, [2, 98])

    # Reference
    ax_ref = axes[i, 0]
    sc = ax_ref.scatter(df_ref["lon"], df_ref["lat"], c=ref_vals, s=0.3,
                        cmap=cmap, vmin=vmin, vmax=vmax,
                        edgecolors="none", rasterized=True)
    ax_ref.set_xlim(-180, 180)
    ax_ref.set_ylim(-90, 90)
    ax_ref.set_title(f"Reference: {label}", fontsize=10)
    plt.colorbar(sc, ax=ax_ref, fraction=0.046, pad=0.04)

    # NEMO
    ax_nemo = axes[i, 1]
    sc2 = ax_nemo.scatter(df_nemo["lon"], df_nemo["lat"], c=nemo_vals, s=0.5,
                          cmap=cmap, vmin=vmin, vmax=vmax,
                          edgecolors="none", rasterized=True)
    ax_nemo.set_xlim(-180, 180)
    ax_nemo.set_ylim(-90, 90)
    ax_nemo.set_title(f"NEMO-PISCES: {label}", fontsize=10)
    plt.colorbar(sc2, ax=ax_nemo, fraction=0.046, pad=0.04)

    # 통계 표시
    ax_ref.text(0.02, 0.02,
                f"n={len(ref_finite)}, mean={np.nanmean(ref_finite):.1f}, "
                f"min={np.nanmin(ref_finite):.1f}, max={np.nanmax(ref_finite):.1f}",
                transform=ax_ref.transAxes, fontsize=7, va="bottom",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    ax_nemo.text(0.02, 0.02,
                 f"n={len(nemo_finite)}, mean={np.nanmean(nemo_finite):.1f}, "
                 f"min={np.nanmin(nemo_finite):.1f}, max={np.nanmax(nemo_finite):.1f}",
                 transform=ax_nemo.transAxes, fontsize=7, va="bottom",
                 bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))

fig.suptitle("FEISTY Input Comparison: Reference vs NEMO-PISCES", fontsize=14, y=1.01)
plt.tight_layout()

path1 = f"{out_dir}/compare_ref_vs_nemo_maps.png"
plt.savefig(path1, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path1}")

# =============================================================
# 플롯 2: 히스토그램 비교
# =============================================================
print("\n[2] Histogram 비교 생성...")
fig, axes = plt.subplots(2, 4, figsize=(20, 8))

for ax, (var, label, _, _, _) in zip(axes.ravel(), variables):
    ref_vals = df_ref[var].dropna().values
    nemo_vals = df_nemo[var].dropna().values

    # 공통 범위
    all_vals = np.concatenate([ref_vals, nemo_vals])
    lo, hi = np.nanpercentile(all_vals, [1, 99])
    bins = np.linspace(lo, hi, 60)

    ax.hist(ref_vals, bins=bins, alpha=0.5, label=f"Ref (n={len(ref_vals)})",
            density=True, color="steelblue")
    ax.hist(nemo_vals, bins=bins, alpha=0.5, label=f"NEMO (n={len(nemo_vals)})",
            density=True, color="tomato")
    ax.set_title(label, fontsize=10)
    ax.legend(fontsize=7)
    ax.set_ylabel("density")

fig.suptitle("Distribution Comparison: Reference vs NEMO-PISCES", fontsize=14)
plt.tight_layout()

path2 = f"{out_dir}/compare_ref_vs_nemo_hist.png"
plt.savefig(path2, dpi=150, bbox_inches="tight")
plt.close()
print(f"  저장: {path2}")

# =============================================================
# 플롯 3: 요약 통계 비교 테이블
# =============================================================
print("\n[3] 요약 통계 비교")
print(f"\n{'변수':>8s} | {'Ref mean':>10s} {'Ref std':>10s} {'Ref min':>10s} {'Ref max':>10s} | {'NEMO mean':>10s} {'NEMO std':>10s} {'NEMO min':>10s} {'NEMO max':>10s} | {'ratio':>6s}")
print("-" * 110)
for var, label, _, _, _ in variables:
    r = df_ref[var].dropna()
    n = df_nemo[var].dropna()
    ratio = n.mean() / r.mean() if r.mean() != 0 else float("inf")
    print(f"{var:>8s} | {r.mean():10.2f} {r.std():10.2f} {r.min():10.2f} {r.max():10.2f} | {n.mean():10.2f} {n.std():10.2f} {n.min():10.2f} {n.max():10.2f} | {ratio:6.2f}x")

print(f"\n완료!")
print(f"  맵 비교: {path1}")
print(f"  히스토그램: {path2}")
