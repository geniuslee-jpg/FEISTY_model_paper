#!/usr/bin/env python3
"""
NEMO-PISCES → FEISTY Input CSV 변환
디버깅 출력 포함 버전
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

# 단위 변환: mol C/m2/s → g wet weight/m2/yr
# = sec_per_yr(365*86400) × g_per_mol(12) × ww_per_C(9)
CONV = 365 * 86400 * 12.0 * 9.0  # = 3,405,888,000

# Martin curve exponent
MARTIN_B = 0.858

print("=" * 60)
print("NEMO-PISCES → FEISTY Input")
print(f"grid: {grid_file}")
print(f"diad: {diad_file}")
print("=" * 60)

# =============================================================
# 파일 로드
# =============================================================
print("\n[0] 파일 로드...")
ds_grid = xr.open_dataset(grid_file, decode_times=False)
ds_diad = xr.open_dataset(diad_file, decode_times=False)

# 시간 차원 자동 감지
time_dim_grid = [d for d in ds_grid["thetao"].dims if "time" in d][0]
time_dim_diad = [d for d in ds_diad["Heup"].dims if "time" in d][0]

print(f"  grid time dim: {time_dim_grid}, steps: {len(ds_grid[time_dim_grid])}")
print(f"  diad time dim: {time_dim_diad}, steps: {len(ds_diad[time_dim_diad])}")
print(f"  grid vars: {list(ds_grid.data_vars)}")
print(f"  diad vars: {list(ds_diad.data_vars)}")

# 좌표
lon2d = ds_grid["nav_lon"].values
lat2d = ds_grid["nav_lat"].values
deptht = ds_grid["deptht"].values
ny, nx = lon2d.shape

print(f"  격자: ({ny}, {nx}), 총 {ny*nx}")
print(f"  deptht ({len(deptht)} levels): {np.round(deptht, 1)}")

# =============================================================
# 레벨 두께 계산
# =============================================================
print("\n[0b] 레벨 두께...")
if "deptht_bounds" in ds_grid:
    dz_bnds = ds_grid["deptht_bounds"].values  # (31, 2)
    layer_dz = dz_bnds[:, 1] - dz_bnds[:, 0]
    level_bottom = dz_bnds[:, 1]
    print("  deptht_bounds 사용")
else:
    # fallback: 중간점 방식
    layer_dz = np.zeros(len(deptht))
    layer_dz[0] = (deptht[0] + deptht[1]) / 2.0
    for k in range(1, len(deptht) - 1):
        layer_dz[k] = (deptht[k+1] - deptht[k-1]) / 2.0
    layer_dz[-1] = layer_dz[-2]
    level_bottom = deptht + layer_dz / 2.0
    print("  deptht_bounds 없음 -> 중간점 방식")

print(f"  layer_dz: {np.round(layer_dz, 1)}")
print(f"  level_bottom: {np.round(level_bottom, 1)}")

dz_da = xr.DataArray(layer_dz, dims=["deptht"],
                      coords={"deptht": deptht})

# =============================================================
# 1. Bathymetry 재구성 (thetao 3D mask)
# =============================================================
print("\n[1] Bathymetry 재구성...")

theta0 = ds_grid["thetao"].isel({time_dim_grid: 0}).values  # (31, y, x)
theta0 = np.where(np.abs(theta0) > 1e10, np.nan, theta0)
valid_mask = ~np.isnan(theta0)  # (31, y, x)

print(f"  theta0 shape: {theta0.shape}")
print(f"  valid per level: {[int(valid_mask[k].sum()) for k in range(len(deptht))]}")

# 깊은 레벨부터 올라오며 가장 깊은 유효 레벨의 하한 = bathymetry
depth_bot = np.full((ny, nx), np.nan)
for k in range(len(deptht) - 1, -1, -1):
    update = np.isnan(depth_bot) & valid_mask[k]
    depth_bot[update] = level_bottom[k]

ocean_mask = ~np.isnan(depth_bot)
n_ocean = ocean_mask.sum()

print(f"  바다 격자: {n_ocean}")
print(f"  depth range: {np.nanmin(depth_bot):.1f} ~ {np.nanmax(depth_bot):.1f} m")

if n_ocean == 0:
    print("  *** ERROR: 바다 격자 0개! thetao 값 확인 필요 ***")
    print(f"  theta0 min/max: {np.nanmin(theta0)}, {np.nanmax(theta0)}")
    print(f"  theta0 NaN count: {np.isnan(theta0).sum()} / {theta0.size}")
    raise ValueError("No ocean grid cells found")

# =============================================================
# 2. Tb: sbt (sea bottom temperature)
# =============================================================
print("\n[2] Tb (sbt)...")

Tb = ds_grid["sbt"].mean(dim=time_dim_grid).values
Tb = np.where(np.abs(Tb) > 1e10, np.nan, Tb)

n_valid_Tb = np.count_nonzero(~np.isnan(Tb))
print(f"  shape: {Tb.shape}, valid: {n_valid_Tb}")
print(f"  range: {np.nanmin(Tb):.2f} ~ {np.nanmax(Tb):.2f} C")

# =============================================================
# 3. Tp: thetao 0-100m 깊이 가중 평균
# =============================================================
print("\n[3] Tp (thetao 0-100m)...")

idx_tp = np.where(deptht <= 100)[0]
print(f"  levels: {np.round(deptht[idx_tp], 1)}")
print(f"  weights: {np.round(layer_dz[idx_tp], 1)}")

wt_tp = xr.DataArray(layer_dz[idx_tp], dims=["deptht"])
temp_tp = ds_grid["thetao"].isel(deptht=idx_tp)
Tp = (temp_tp * wt_tp).sum(dim="deptht") / wt_tp.sum()
Tp = Tp.mean(dim=time_dim_grid).values
Tp = np.where(np.abs(Tp) > 1e10, np.nan, Tp)

n_valid_Tp = np.count_nonzero(~np.isnan(Tp))
print(f"  shape: {Tp.shape}, valid: {n_valid_Tp}")
print(f"  range: {np.nanmin(Tp):.2f} ~ {np.nanmax(Tp):.2f} C")

# =============================================================
# 4. Tm: thetao 500-1500m 깊이 가중 평균
# =============================================================
print("\n[4] Tm (thetao 500-1500m)...")

idx_tm = np.where((deptht >= 500) & (deptht <= 1500))[0]
print(f"  levels: {np.round(deptht[idx_tm], 1)}")

if len(idx_tm) > 0:
    wt_tm = xr.DataArray(layer_dz[idx_tm], dims=["deptht"])
    temp_tm = ds_grid["thetao"].isel(deptht=idx_tm)
    Tm = (temp_tm * wt_tm).sum(dim="deptht") / wt_tm.sum()
    Tm = Tm.mean(dim=time_dim_grid).values
    Tm = np.where(np.abs(Tm) > 1e10, np.nan, Tm)
else:
    print("  WARNING: 500-1500m 없음 -> 최하층 사용")
    Tm = ds_grid["thetao"].isel(deptht=-1).mean(dim=time_dim_grid).values
    Tm = np.where(np.abs(Tm) > 1e10, np.nan, Tm)

# 얕은 해역 (500m 미만): Tm NaN -> Tb로 대체
n_shallow = int((np.isnan(Tm) & ocean_mask).sum())
Tm = np.where(np.isnan(Tm) & ocean_mask, Tb, Tm)
print(f"  얕은 해역 Tb 대체: {n_shallow}개")
print(f"  range: {np.nanmin(Tm):.2f} ~ {np.nanmax(Tm):.2f} C")

ds_grid.close()

# =============================================================
# 5. Heup -> photic (유광층 깊이)
# =============================================================
print("\n[5] photic (Heup)...")

photic = ds_diad["Heup"].mean(dim=time_dim_diad).values
photic = np.where(np.abs(photic) > 1e10, np.nan, photic)

print(f"  shape: {photic.shape}")
print(f"  range: {np.nanmin(photic):.1f} ~ {np.nanmax(photic):.1f} m")

# =============================================================
# 6. EPC100 -> dfbot (Martin curve)
# =============================================================
print("\n[6] dfbot (EPC100 -> Martin curve)...")

epc100 = ds_diad["EPC100"].mean(dim=time_dim_diad).values  # mol C/m2/s
epc100 = np.where(np.abs(epc100) > 1e10, np.nan, epc100)

# 단위 변환: mol C/m2/s -> g ww/m2/yr at 100m
epc100_gww = epc100 * CONV
print(f"  EPC100 (mol C/m2/s): {np.nanmin(epc100):.2e} ~ {np.nanmax(epc100):.2e}")
print(f"  EPC100 (g ww/m2/yr at 100m): {np.nanmin(epc100_gww):.2f} ~ {np.nanmax(epc100_gww):.2f}")

# Martin curve: F(z_bot) = F(100m) x (z_bot/100)^(-b)
depth_martin = np.maximum(depth_bot, 100.0)
dfbot = epc100_gww * np.power(depth_martin / 100.0, -MARTIN_B)
print(f"  dfbot (at seafloor): {np.nanmin(dfbot):.4f} ~ {np.nanmax(dfbot):.2f}")

# =============================================================
# 7. GRAZ1 -> szprod, GRAZ2 -> lzprod (수직적분 + 단위 변환)
# =============================================================
print("\n[7] szprod/lzprod (GRAZ1/GRAZ2)...")

graz1_mean = ds_diad["GRAZ1"].mean(dim=time_dim_diad)  # (deptht, y, x) mol/m3/s
graz2_mean = ds_diad["GRAZ2"].mean(dim=time_dim_diad)

print(f"  GRAZ1 after time mean: {graz1_mean.shape} (should be 3D: deptht,y,x)")
print(f"  GRAZ2 after time mean: {graz2_mean.shape}")

# 수직적분: sum( GRAZ(k) x dz(k) ) -> mol C/m2/s
graz1_int = (graz1_mean * dz_da).sum(dim="deptht").values  # (y, x)
graz2_int = (graz2_mean * dz_da).sum(dim="deptht").values

graz1_int = np.where(np.abs(graz1_int) > 1e10, np.nan, graz1_int)
graz2_int = np.where(np.abs(graz2_int) > 1e10, np.nan, graz2_int)

print(f"  after vertical integration: {graz1_int.shape} (should be 2D: y,x)")

# mol C/m2/s -> g ww/m2/yr
szprod = graz1_int * CONV
lzprod = graz2_int * CONV

print(f"  szprod: {np.nanmin(szprod):.2f} ~ {np.nanmax(szprod):.2f} g ww/m2/yr")
print(f"  lzprod: {np.nanmin(lzprod):.2f} ~ {np.nanmax(lzprod):.2f} g ww/m2/yr")

ds_diad.close()

# =============================================================
# 8. DataFrame 생성
# =============================================================
print("\n[8] DataFrame 생성...")

# shape 검증
arrays = {
    "lon": lon2d, "lat": lat2d, "depth": depth_bot,
    "Tp": Tp, "Tm": Tm, "Tb": Tb,
    "photic": photic, "dfbot": dfbot,
    "szprod": szprod, "lzprod": lzprod,
}
print("  Shape 검증:")
for name, arr in arrays.items():
    print(f"    {name:8s}: {arr.shape}")
    assert arr.shape == (ny, nx), f"{name} shape mismatch! expected ({ny},{nx})"

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
        print(f"  {col}: 음수 {n_neg}개 -> 0")
        df[col] = df[col].clip(lower=0)

# NaN 잔여 확인
print("\n  NaN 개수:")
for col in df.columns:
    n_nan = df[col].isna().sum()
    if n_nan > 0:
        print(f"    {col}: {n_nan}")

# =============================================================
# 9. 저장
# =============================================================
out_file = f"{OUT_DIR}/Input_NEMO_FEISTY.csv"
df.to_csv(out_file, index=False)

print(f"\n[9] 저장: {out_file}")
print(f"  격자 수: {len(df)}")

# =============================================================
# 10. 요약 통계
# =============================================================
print(f"\n{'='*60}")
print("데이터 요약")
print(f"{'='*60}")
print(df.describe().round(3).to_string())

# =============================================================
# 11. 검증 플롯
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

    fig.suptitle("NEMO-PISCES -> FEISTY Input (10yr mean)", fontsize=13)
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
