#!/usr/bin/env python3
"""
COBALT 원본 인풋 vs EXP02 인풋 동물플랑크톤 파라미터 비교
- 변수: szprod, lzprod, dfbot, photic, depth, Tb, Tm, Tp
- Mollweide 투영 + 1°×1° 그리드
- Side-by-side 맵 + Bias + 히스토그램 + Zonal mean
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
cobalt_file = "/data01/labdisk/sungjin/NEMO-FEISTY/00.org_COBALR/Input_global.csv"
exp02_file  = "/data01/labdisk/sungjin/NEMO-FEISTY/06.ratio_application_input/data/EXP_02/Input_NEMO_FEISTY_corrected.csv"
out_dir     = "/data01/labdisk/sungjin/NEMO-FEISTY/06.ratio_application_input/figures"

# =============================================================
# 데이터 로드
# =============================================================
print("Loading data...")
df_cobalt = pd.read_csv(cobalt_file)
df_exp02  = pd.read_csv(exp02_file)

# lon 보정: 0~360 → -180~180
for df in [df_cobalt, df_exp02]:
    df["lon"] = np.where(df["lon"] > 180, df["lon"] - 360, df["lon"])

print(f"  COBALT (original): {len(df_cobalt)} points, columns: {list(df_cobalt.columns)}")
print(f"  EXP02 (corrected): {len(df_exp02)} points, columns: {list(df_exp02.columns)}")

# =============================================================
# 1° × 1° 정규 격자
# =============================================================
lon_edges = np.arange(-180, 181, 1)
lat_edges = np.arange(-90, 91, 1)
lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2
LON, LAT = np.meshgrid(lon_centers, lat_centers)


def to_grid_binning(df, var):
    """~1° 정규격자 → 1°×1° binning (COBALT 용)."""
    vals = df[var].values
    mask = np.isfinite(vals)
    if mask.sum() == 0:
        return np.full(LON.shape, np.nan)
    result = binned_statistic_2d(
        df["lon"].values[mask], df["lat"].values[mask], vals[mask],
        statistic="mean", bins=[lon_edges, lat_edges]
    )
    return result.statistic.T


def to_grid_interp(df, var):
    """NEMO ORCA2 (~2° 비정규격자) → 1°×1° nearest 보간."""
    vals = df[var].values
    mask = np.isfinite(vals)
    if mask.sum() == 0:
        return np.full(LON.shape, np.nan)
    points = np.column_stack([df["lon"].values[mask], df["lat"].values[mask]])
    grid = griddata(points, vals[mask], (LON, LAT), method="nearest")
    tree = cKDTree(points)
    dist, _ = tree.query(np.column_stack([LON.ravel(), LAT.ravel()]))
    grid[dist.reshape(LON.shape) > 2.5] = np.nan
    return grid


# EXP02가 ORCA2 비정규격자인지 확인 후 적절한 함수 결정
# COBALT는 1° 정규격자, EXP02도 동일 1° 격자로 변환되었을 수 있음
# 격자 간격으로 판단
def detect_grid_type(df):
    """격자 타입 추정: regular (~1°) vs irregular (ORCA2)."""
    lon_unique = np.sort(df["lon"].unique())
    if len(lon_unique) > 10:
        diffs = np.diff(lon_unique[:20])
        mean_diff = np.mean(diffs)
        if 0.8 < mean_diff < 1.2:
            return "regular"
    return "irregular"


cobalt_grid_type = detect_grid_type(df_cobalt)
exp02_grid_type  = detect_grid_type(df_exp02)
print(f"  COBALT grid type: {cobalt_grid_type}")
print(f"  EXP02  grid type: {exp02_grid_type}")

grid_cobalt = to_grid_binning if cobalt_grid_type == "regular" else to_grid_interp
grid_exp02  = to_grid_binning if exp02_grid_type == "regular" else to_grid_interp

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
# 비교 변수 설정
# =============================================================
# 동물플랑크톤 관련 변수 (log10 스케일)
zoo_variables = [
    ("szprod", "Small Zoo Production\n[g WW m⁻² yr⁻¹]", "YlOrRd", True),
    ("lzprod", "Large Zoo Production\n[g WW m⁻² yr⁻¹]", "YlOrRd", True),
    ("dfbot",  "Detrital Flux to Bottom\n[g WW m⁻² yr⁻¹]", "YlGnBu", True),
]

# 물리 변수 (선형 스케일)
phys_variables = [
    ("Tp",     "Pelagic Temp (0-100m)\n[°C]", "RdYlBu_r", False),
    ("Tm",     "Mid-water Temp (500-1500m)\n[°C]", "RdYlBu_r", False),
    ("Tb",     "Bottom Temp\n[°C]", "RdYlBu_r", False),
    ("photic", "Euphotic Depth\n[m]", "viridis", False),
    ("depth",  "Seafloor Depth\n[m]", "cividis", False),
]

# 두 데이터셋에서 공통으로 존재하는 변수만 선택
common_cols = set(df_cobalt.columns) & set(df_exp02.columns)
zoo_variables  = [(v, l, c, lg) for v, l, c, lg in zoo_variables  if v in common_cols]
phys_variables = [(v, l, c, lg) for v, l, c, lg in phys_variables if v in common_cols]

all_variables = zoo_variables + phys_variables
print(f"\n  비교할 변수: {[v[0] for v in all_variables]}")

# =============================================================
# [1] Side-by-side 맵 비교 + Bias (3열: COBALT | EXP02 | Bias)
# =============================================================
print("\n[1] Side-by-side maps with bias (3-column)...")
nvar = len(all_variables)
fig, axes = plt.subplots(nvar, 3, figsize=(24, 4.5 * nvar),
                         subplot_kw={"projection": PROJ})
if nvar == 1:
    axes = axes.reshape(1, -1)

for i, (var, label, cmap, use_log) in enumerate(all_variables):
    print(f"  Processing {var}...")

    cobalt_grid = grid_cobalt(df_cobalt, var)
    exp02_grid  = grid_exp02(df_exp02, var)

    if use_log:
        with np.errstate(divide="ignore", invalid="ignore"):
            cobalt_plot = np.log10(np.maximum(cobalt_grid, 1e-4))
            exp02_plot  = np.log10(np.maximum(exp02_grid, 1e-4))
        unit_label = "log₁₀(g WW m⁻² yr⁻¹)"
        # 자동 범위 설정
        all_vals = np.concatenate([cobalt_plot[np.isfinite(cobalt_plot)],
                                   exp02_plot[np.isfinite(exp02_plot)]])
        vmin = np.nanpercentile(all_vals, 2)
        vmax = np.nanpercentile(all_vals, 98)
    else:
        cobalt_plot = cobalt_grid
        exp02_plot  = exp02_grid
        unit_label = label.split("\n")[-1].strip("[]") if "\n" in label else ""
        all_vals = np.concatenate([cobalt_plot[np.isfinite(cobalt_plot)],
                                   exp02_plot[np.isfinite(exp02_plot)]])
        vmin = np.nanpercentile(all_vals, 2)
        vmax = np.nanpercentile(all_vals, 98)

    # COBALT
    ax = axes[i, 0]
    pc = plot_on_mollweide(ax, cobalt_plot, cmap, vmin, vmax)
    ax.set_title(f"COBALT: {label}", fontsize=10, fontweight="bold")
    cob_f = cobalt_plot[np.isfinite(cobalt_plot)]
    ax.text(0.02, 0.02,
            f"mean={np.nanmean(cob_f):.2f}, "
            f"min={np.nanmin(cob_f):.2f}, max={np.nanmax(cob_f):.2f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label(unit_label, fontsize=8)

    # EXP02
    ax = axes[i, 1]
    pc = plot_on_mollweide(ax, exp02_plot, cmap, vmin, vmax)
    ax.set_title(f"EXP02 (NEMO-PISCES corrected): {label}", fontsize=10, fontweight="bold")
    exp_f = exp02_plot[np.isfinite(exp02_plot)]
    ax.text(0.02, 0.02,
            f"mean={np.nanmean(exp_f):.2f}, "
            f"min={np.nanmin(exp_f):.2f}, max={np.nanmax(exp_f):.2f}",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label(unit_label, fontsize=8)

    # Bias (EXP02 - COBALT)
    ax = axes[i, 2]
    bias_grid = exp02_plot - cobalt_plot
    mask = np.isfinite(cobalt_plot) & np.isfinite(exp02_plot)
    bias_grid[~mask] = np.nan
    bias_f = bias_grid[np.isfinite(bias_grid)]
    if len(bias_f) > 0:
        blim = max(np.nanpercentile(np.abs(bias_f), 98), 0.01)
    else:
        blim = 1.0
    pc = plot_on_mollweide(ax, bias_grid, "RdBu_r", -blim, blim)
    ax.set_title(f"Bias (EXP02 − COBALT): {var}", fontsize=10, fontweight="bold")
    if len(bias_f) > 0:
        ax.text(0.02, 0.02,
                f"mean={np.nanmean(bias_f):.2f}, "
                f"std={np.nanstd(bias_f):.2f}",
                transform=ax.transAxes, fontsize=7, va="bottom",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    if use_log:
        cb_label = "Δ log₁₀"
    else:
        cb_label = f"Δ {unit_label}"
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label(cb_label, fontsize=8)

fig.suptitle("FEISTY Input Comparison: COBALT (original) vs EXP02 (NEMO-PISCES corrected)",
             fontsize=16, y=1.01, fontweight="bold")
plt.tight_layout()
path1 = f"{out_dir}/fig_input_cobalt_vs_exp02_maps.png"
plt.savefig(path1, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {path1}")

# =============================================================
# [2] 동물플랑크톤 변수만 확대 비교 (3열: COBALT | EXP02 | Bias)
# =============================================================
if zoo_variables:
    print("\n[2] Zooplankton-focused maps...")
    nzoo = len(zoo_variables)
    fig, axes = plt.subplots(nzoo, 3, figsize=(24, 5 * nzoo),
                             subplot_kw={"projection": PROJ})
    if nzoo == 1:
        axes = axes.reshape(1, -1)

    for i, (var, label, cmap, use_log) in enumerate(zoo_variables):
        cobalt_grid = grid_cobalt(df_cobalt, var)
        exp02_grid  = grid_exp02(df_exp02, var)

        with np.errstate(divide="ignore", invalid="ignore"):
            cobalt_plot = np.log10(np.maximum(cobalt_grid, 1e-4))
            exp02_plot  = np.log10(np.maximum(exp02_grid, 1e-4))

        all_vals = np.concatenate([cobalt_plot[np.isfinite(cobalt_plot)],
                                   exp02_plot[np.isfinite(exp02_plot)]])
        vmin = np.nanpercentile(all_vals, 1)
        vmax = np.nanpercentile(all_vals, 99)

        # COBALT
        ax = axes[i, 0]
        pc = plot_on_mollweide(ax, cobalt_plot, cmap, vmin, vmax)
        ax.set_title(f"COBALT: {label}", fontsize=12, fontweight="bold")
        cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
        cb.set_label("log₁₀(g WW m⁻² yr⁻¹)", fontsize=9)

        # EXP02
        ax = axes[i, 1]
        pc = plot_on_mollweide(ax, exp02_plot, cmap, vmin, vmax)
        ax.set_title(f"EXP02: {label}", fontsize=12, fontweight="bold")
        cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
        cb.set_label("log₁₀(g WW m⁻² yr⁻¹)", fontsize=9)

        # Bias
        ax = axes[i, 2]
        bias_grid = exp02_plot - cobalt_plot
        mask = np.isfinite(cobalt_plot) & np.isfinite(exp02_plot)
        bias_grid[~mask] = np.nan
        bias_f = bias_grid[np.isfinite(bias_grid)]
        blim = max(np.nanpercentile(np.abs(bias_f), 98), 0.01) if len(bias_f) > 0 else 1.0
        pc = plot_on_mollweide(ax, bias_grid, "RdBu_r", -blim, blim)
        ax.set_title(f"Bias (EXP02 − COBALT): {var}", fontsize=12, fontweight="bold")
        if len(bias_f) > 0:
            ax.text(0.02, 0.02,
                    f"mean bias={np.nanmean(bias_f):.3f}, "
                    f"RMSE={np.sqrt(np.nanmean(bias_f**2)):.3f}",
                    transform=ax.transAxes, fontsize=8, va="bottom",
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
        cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
        cb.set_label("Δ log₁₀", fontsize=9)

    fig.suptitle("Zooplankton Input: COBALT vs EXP02",
                 fontsize=16, y=1.01, fontweight="bold")
    plt.tight_layout()
    path2 = f"{out_dir}/fig_input_zoo_cobalt_vs_exp02.png"
    plt.savefig(path2, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {path2}")

# =============================================================
# [3] 히스토그램 비교
# =============================================================
print("\n[3] Histograms...")
ncols = min(4, len(all_variables))
nrows = int(np.ceil(len(all_variables) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
axes = np.atleast_2d(axes)

for idx, (var, label, cmap, use_log) in enumerate(all_variables):
    ax = axes.ravel()[idx]
    cv = df_cobalt[var].dropna().values
    ev = df_exp02[var].dropna().values

    if use_log:
        cv = np.log10(cv[cv > 0])
        ev = np.log10(ev[ev > 0])
        xlabel = f"log₁₀({var})"
    else:
        xlabel = var

    lo = np.nanpercentile(np.concatenate([cv, ev]), 1)
    hi = np.nanpercentile(np.concatenate([cv, ev]), 99)
    bins = np.linspace(lo, hi, 50)

    ax.hist(cv, bins=bins, alpha=0.55, density=True, color="steelblue",
            label=f"COBALT (n={len(cv)})", edgecolor="white", linewidth=0.3)
    ax.hist(ev, bins=bins, alpha=0.55, density=True, color="tomato",
            label=f"EXP02 (n={len(ev)})", edgecolor="white", linewidth=0.3)

    short_label = label.split("\n")[0] if "\n" in label else label
    ax.set_title(short_label, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel("Density", fontsize=9)
    ax.legend(fontsize=7, loc="upper right")

# 빈 서브플롯 숨기기
for idx in range(len(all_variables), nrows * ncols):
    axes.ravel()[idx].set_visible(False)

fig.suptitle("Input Distribution: COBALT vs EXP02", fontsize=14, fontweight="bold")
plt.tight_layout()
path3 = f"{out_dir}/fig_input_hist_cobalt_vs_exp02.png"
plt.savefig(path3, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {path3}")

# =============================================================
# [4] Zonal Mean 비교
# =============================================================
print("\n[4] Zonal mean comparison...")
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
axes = axes.ravel()

for idx, (var, label, cmap, use_log) in enumerate(all_variables):
    if idx >= len(axes):
        break
    ax = axes[idx]

    cobalt_grid = grid_cobalt(df_cobalt, var)
    exp02_grid  = grid_exp02(df_exp02, var)

    if use_log:
        with np.errstate(divide="ignore", invalid="ignore"):
            cobalt_grid = np.log10(np.maximum(cobalt_grid, 1e-4))
            exp02_grid  = np.log10(np.maximum(exp02_grid, 1e-4))

    cobalt_zonal = np.nanmean(cobalt_grid, axis=1)
    exp02_zonal  = np.nanmean(exp02_grid, axis=1)

    ax.plot(cobalt_zonal, lat_centers, color="steelblue", linewidth=1.5,
            label="COBALT", alpha=0.9)
    ax.plot(exp02_zonal, lat_centers, color="tomato", linewidth=1.5,
            label="EXP02", alpha=0.9)
    ax.fill_betweenx(lat_centers, cobalt_zonal, exp02_zonal,
                     alpha=0.15, color="grey")

    short_label = label.split("\n")[0] if "\n" in label else label
    ax.set_title(short_label, fontsize=10)
    if use_log:
        ax.set_xlabel(f"log₁₀({var})", fontsize=9)
    else:
        ax.set_xlabel(var, fontsize=9)
    ax.set_ylabel("Latitude", fontsize=9)
    ax.set_ylim(-90, 90)
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

for idx in range(len(all_variables), len(axes)):
    axes[idx].set_visible(False)

fig.suptitle("Zonal Mean: COBALT vs EXP02", fontsize=14, fontweight="bold")
plt.tight_layout()
path4 = f"{out_dir}/fig_input_zonal_cobalt_vs_exp02.png"
plt.savefig(path4, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved: {path4}")

# =============================================================
# [5] 요약 통계 테이블
# =============================================================
print(f"\n{'='*120}")
print(f"{'Variable':>10s} | {'COBALT mean':>12s} {'COBALT std':>12s} {'COBALT min':>12s} {'COBALT max':>12s} | "
      f"{'EXP02 mean':>12s} {'EXP02 std':>12s} {'EXP02 min':>12s} {'EXP02 max':>12s} | {'ratio':>6s}")
print("-" * 120)
for var, label, *_ in all_variables:
    c = df_cobalt[var].dropna()
    e = df_exp02[var].dropna()
    ratio = e.mean() / c.mean() if c.mean() != 0 else float("inf")
    print(f"{var:>10s} | {c.mean():12.4f} {c.std():12.4f} {c.min():12.4f} {c.max():12.4f} | "
          f"{e.mean():12.4f} {e.std():12.4f} {e.min():12.4f} {e.max():12.4f} | {ratio:5.2f}x")
print(f"{'='*120}")

print(f"\nDone!")
print(f"  [1] All variables map:  {path1}")
if zoo_variables:
    print(f"  [2] Zoo variables map:  {path2}")
print(f"  [3] Histograms:         {path3}")
print(f"  [4] Zonal means:        {path4}")
