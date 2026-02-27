#!/usr/bin/env python3
"""
동물플랑크톤 인풋 비교: COBALT vs EXP02
  - 2행: szprod (소형), lzprod (중대형)
  - 3열: COBALT | EXP02 | EXP02 - COBALT
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

# ---- 파일 경로 ----
cobalt_file = "/data01/labdisk/sungjin/NEMO-FEISTY/00.org_COBALR/Input_global.csv"
exp02_file  = "/data01/labdisk/sungjin/NEMO-FEISTY/06.ratio_application_input/data/EXP_02/Input_NEMO_FEISTY_corrected.csv"
out_dir     = "/data01/labdisk/sungjin/NEMO-FEISTY/06.ratio_application_input/figures"

# ---- 데이터 로드 ----
print("Loading data...")
df_cobalt = pd.read_csv(cobalt_file)
df_exp02  = pd.read_csv(exp02_file)

for df in [df_cobalt, df_exp02]:
    df["lon"] = np.where(df["lon"] > 180, df["lon"] - 360, df["lon"])

# ---- 1x1 정규격자 ----
lon_edges   = np.arange(-180, 181, 1)
lat_edges   = np.arange(-90, 91, 1)
lon_centers = (lon_edges[:-1] + lon_edges[1:]) / 2
lat_centers = (lat_edges[:-1] + lat_edges[1:]) / 2
LON, LAT    = np.meshgrid(lon_centers, lat_centers)


def to_grid_binning(df, var):
    vals = df[var].values
    mask = np.isfinite(vals)
    if mask.sum() == 0:
        return np.full(LON.shape, np.nan)
    result = binned_statistic_2d(
        df["lon"].values[mask], df["lat"].values[mask], vals[mask],
        statistic="mean", bins=[lon_edges, lat_edges])
    return result.statistic.T


def to_grid_interp(df, var):
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


def detect_grid_type(df):
    lon_unique = np.sort(df["lon"].unique())
    if len(lon_unique) > 10:
        diffs = np.diff(lon_unique[:20])
        if 0.8 < np.mean(diffs) < 1.2:
            return "regular"
    return "irregular"


grid_cobalt = to_grid_binning if detect_grid_type(df_cobalt) == "regular" else to_grid_interp
grid_exp02  = to_grid_binning if detect_grid_type(df_exp02)  == "regular" else to_grid_interp

# ---- 투영 설정 ----
PROJ     = ccrs.Mollweide(central_longitude=0)
DATA_CRS = ccrs.PlateCarree()


def add_map_features(ax):
    ax.set_global()
    ax.coastlines(linewidth=0.4, color="k")
    ax.add_feature(cfeature.LAND, facecolor="#d0d0d0", edgecolor="none", zorder=2)
    ax.gridlines(draw_labels=False, linewidth=0.2, color="grey",
                 alpha=0.5, linestyle="--")


def plot_map(ax, grid, cmap, vmin, vmax):
    data_c, lon_c = add_cyclic_point(grid, coord=lon_centers)
    pc = ax.pcolormesh(lon_c, lat_centers, data_c, transform=DATA_CRS,
                       cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="auto", rasterized=True)
    add_map_features(ax)
    return pc


# ---- 변수 정의: 2행 ----
variables = [
    ("szprod", "Small Zooplankton Production", "YlOrRd"),
    ("lzprod", "Large Zooplankton Production", "YlOrRd"),
]

# ---- 그림 그리기: 2행 x 3열 ----
print("Plotting...")
fig, axes = plt.subplots(2, 3, figsize=(24, 10), subplot_kw={"projection": PROJ})

for i, (var, label, cmap) in enumerate(variables):
    print(f"  {var}...")
    cg = grid_cobalt(df_cobalt, var)
    eg = grid_exp02(df_exp02, var)

    # log10 변환
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.log10(np.maximum(cg, 1e-4))
        ep = np.log10(np.maximum(eg, 1e-4))

    # 공통 색상 범위
    all_vals = np.concatenate([cp[np.isfinite(cp)], ep[np.isfinite(ep)]])
    vmin = np.nanpercentile(all_vals, 2)
    vmax = np.nanpercentile(all_vals, 98)
    unit = "log$_{10}$(g WW m$^{-2}$ yr$^{-1}$)"

    # (i, 0) COBALT
    ax = axes[i, 0]
    pc = plot_map(ax, cp, cmap, vmin, vmax)
    ax.set_title(f"COBALT: {label}", fontsize=11, fontweight="bold")
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label(unit, fontsize=9)

    # (i, 1) EXP02
    ax = axes[i, 1]
    pc = plot_map(ax, ep, cmap, vmin, vmax)
    ax.set_title(f"EXP02: {label}", fontsize=11, fontweight="bold")
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label(unit, fontsize=9)

    # (i, 2) Bias (EXP02 - COBALT)
    ax = axes[i, 2]
    bias = ep - cp
    mask = np.isfinite(cp) & np.isfinite(ep)
    bias[~mask] = np.nan
    bf = bias[np.isfinite(bias)]
    blim = max(np.nanpercentile(np.abs(bf), 98), 0.01) if len(bf) > 0 else 1.0
    pc = plot_map(ax, bias, "RdBu_r", -blim, blim)
    ax.set_title(f"EXP02 $-$ COBALT: {label}", fontsize=11, fontweight="bold")
    if len(bf) > 0:
        ax.text(0.02, 0.02,
                f"mean={np.nanmean(bf):.3f}, RMSE={np.sqrt(np.nanmean(bf**2)):.3f}",
                transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.8), zorder=5)
    cb = plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.04, orientation="horizontal")
    cb.set_label("$\\Delta$ log$_{10}$", fontsize=9)

fig.suptitle("Zooplankton Input: COBALT (original) vs EXP02 (NEMO-PISCES corrected)",
             fontsize=14, fontweight="bold", y=1.02)
plt.tight_layout()
outpath = f"{out_dir}/fig_zoo_input_cobalt_vs_exp02.png"
plt.savefig(outpath, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved: {outpath}")
