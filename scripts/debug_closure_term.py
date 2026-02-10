#!/usr/bin/env python3
"""
Closure term (mzrat × ZOO²) 디버깅 스크립트
ZOO/ZOO2 값 확인 + closure term 추정 + 온도별 보정 + Reference 비교
"""
import xarray as xr
import numpy as np

BASE_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY"
grid_file = f"{BASE_DIR}/02.subset/grid_T_10yr.nc"
zoo_file = f"{BASE_DIR}/01.org/ORCA2_1y_00010101_00501231_ptrc_T.nc"
domcfg_file = f"{BASE_DIR}/01.org/ORCA_R2_zps_domcfg.nc"
ref_file = f"{BASE_DIR}/03.feisty_input/Input_global.csv"

# namelist_pisces_ref 값
MZRAT = 0.02     # micro (/day/(mmol/m3))
MZRAT2 = 0.01    # meso  (/day/(mmol/m3))
CONV_CLOSURE = 1e-3 * 365 * 12.0 * 9.0  # = 39.42

# 온도별 보정 비율 (Stock et al. 2014, 2017)
# ratio = loss_to_fish / max_production
T_KNOTS = np.array([0.0, 10.0, 20.0, 30.0])
R_KNOTS = np.array([0.78, 0.85, 0.89, 0.91])

print("=" * 60)
print("Closure term 디버깅 (온도 보정 포함)")
print("=" * 60)

# ─────────────────────────────────────────────
# [1] ZOO/ZOO2 원시값 확인
# ─────────────────────────────────────────────
print("\n[1] ZOO/ZOO2 원시값...")
ds = xr.open_dataset(zoo_file, decode_times=False)
time_dim = [d for d in ds["ZOO"].dims if "time" in d][0]
print(f"  시간 스텝: {len(ds[time_dim])}")

zoo_mean = ds["ZOO"].mean(dim=time_dim).values
zoo2_mean = ds["ZOO2"].mean(dim=time_dim).values
zoo_mean = np.where(np.abs(zoo_mean) > 1e10, np.nan, zoo_mean)
zoo2_mean = np.where(np.abs(zoo2_mean) > 1e10, np.nan, zoo2_mean)

print(f"\n  ZOO (micro, mmol/m3):")
print(f"    shape: {zoo_mean.shape}")
print(f"    min:   {np.nanmin(zoo_mean):.6f}")
print(f"    max:   {np.nanmax(zoo_mean):.6f}")
print(f"    mean:  {np.nanmean(zoo_mean):.6f}")

print(f"\n  ZOO2 (meso, mmol/m3):")
print(f"    min:   {np.nanmin(zoo2_mean):.6f}")
print(f"    max:   {np.nanmax(zoo2_mean):.6f}")
print(f"    mean:  {np.nanmean(zoo2_mean):.6f}")

print(f"\n  표층 (레벨 0):")
print(f"    ZOO  mean: {np.nanmean(zoo_mean[0]):.6f} mmol/m3")
print(f"    ZOO2 mean: {np.nanmean(zoo2_mean[0]):.6f} mmol/m3")

print(f"\n  레벨별 ZOO 평균 (상위 10):")
for k in range(min(10, zoo_mean.shape[0])):
    z_val = np.nanmean(zoo_mean[k])
    z2_val = np.nanmean(zoo2_mean[k])
    print(f"    lev {k:2d}: ZOO={z_val:.6f}  ZOO2={z2_val:.6f}")

ds.close()

# ─────────────────────────────────────────────
# [2] Tp 계산 (0-100m 깊이가중 평균, 보정비율에 필요)
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[2] Tp (0-100m 평균 수온) 계산...")

ds_grid = xr.open_dataset(grid_file, decode_times=False)
time_dim_grid = [d for d in ds_grid["thetao"].dims if "time" in d][0]
deptht = ds_grid["deptht"].values
nz_grid = len(deptht)

ds_dom = xr.open_dataset(domcfg_file)
bottom_level = ds_dom["bottom_level"].values
e3t_0 = ds_dom["e3t_0"].values
ds_dom.close()

nz = e3t_0.shape[0]
ny, nx = bottom_level.shape
ocean_mask = bottom_level > 0
k3d = np.arange(nz)[:, None, None]
tmask = (k3d < bottom_level[None, :, :]).astype(np.float64)

thetao_mean = ds_grid["thetao"].mean(dim=time_dim_grid).values
thetao_mean = np.where(np.abs(thetao_mean) > 1e10, np.nan, thetao_mean)
thetao_masked = np.where(tmask == 1, thetao_mean, np.nan)

idx_tp = np.where(deptht <= 100)[0]
e3t_tp = e3t_0[idx_tp] * tmask[idx_tp]
theta_tp = thetao_masked[idx_tp]
weight_sum = np.nansum(e3t_tp, axis=0)
weighted_theta = np.nansum(theta_tp * e3t_tp, axis=0)
Tp = np.where(weight_sum > 0, weighted_theta / weight_sum, np.nan)
Tp[~ocean_mask] = np.nan

ds_grid.close()

print(f"  Tp range: {np.nanmin(Tp):.2f} ~ {np.nanmax(Tp):.2f} C")
print(f"  Tp mean:  {np.nanmean(Tp):.2f} C")

# ─────────────────────────────────────────────
# [3] 전구 closure term 계산 (수직적분)
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[3] 전구 closure term (수직적분)...")

closure_zoo = MZRAT * zoo_mean**2
closure_zoo2 = MZRAT2 * zoo2_mean**2

szprod = np.nansum(closure_zoo * e3t_0 * tmask, axis=0) * CONV_CLOSURE
lzprod_raw = np.nansum(closure_zoo2 * e3t_0 * tmask, axis=0) * CONV_CLOSURE

print(f"\n  szprod (gww/m2/yr):")
print(f"    min:  {np.nanmin(szprod[ocean_mask]):.4f}")
print(f"    max:  {np.nanmax(szprod[ocean_mask]):.4f}")
print(f"    mean: {np.nanmean(szprod[ocean_mask]):.4f}")

print(f"\n  lzprod 보정 전 (gww/m2/yr):")
print(f"    min:  {np.nanmin(lzprod_raw[ocean_mask]):.4f}")
print(f"    max:  {np.nanmax(lzprod_raw[ocean_mask]):.4f}")
print(f"    mean: {np.nanmean(lzprod_raw[ocean_mask]):.4f}")

# ─────────────────────────────────────────────
# [4] 온도별 보정 적용
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[4] lzprod 온도별 보정...")
print(f"  ratio = loss/maxprod → lzprod = closure / ratio(Tp)")
print(f"  0C→0.78, 10C→0.85, 20C→0.89, 30C→0.91")

ratio_map = np.interp(np.clip(Tp, 0, 30), T_KNOTS, R_KNOTS)
ratio_map = np.where(ocean_mask, ratio_map, np.nan)

print(f"\n  ratio_map:")
print(f"    min:  {np.nanmin(ratio_map):.4f}")
print(f"    max:  {np.nanmax(ratio_map):.4f}")
print(f"    mean: {np.nanmean(ratio_map):.4f}")

lzprod = lzprod_raw / ratio_map
szprod[~ocean_mask] = np.nan
lzprod[~ocean_mask] = np.nan

print(f"\n  lzprod 보정 후 (gww/m2/yr):")
print(f"    min:  {np.nanmin(lzprod):.4f}")
print(f"    max:  {np.nanmax(lzprod):.4f}")
print(f"    mean: {np.nanmean(lzprod):.4f}")

# ─────────────────────────────────────────────
# [5] Reference와 비교
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[5] Reference 비교")

try:
    import pandas as pd
    ref = pd.read_csv(ref_file)

    ref_sz = ref["szprod"].mean()
    ref_lz = ref["lzprod"].mean()
    nemo_sz = np.nanmean(szprod)
    nemo_lz_raw = np.nanmean(lzprod_raw[ocean_mask])
    nemo_lz = np.nanmean(lzprod)

    print(f"\n  {'':16s} {'Ref':>10s} {'NEMO':>10s} {'비율':>8s}")
    print(f"  {'-'*48}")
    print(f"  {'szprod':16s} {ref_sz:10.2f} {nemo_sz:10.2f} {nemo_sz/ref_sz:7.2f}x")
    print(f"  {'lzprod(보정전)':16s} {ref_lz:10.2f} {nemo_lz_raw:10.2f} {nemo_lz_raw/ref_lz:7.2f}x")
    print(f"  {'lzprod(보정후)':16s} {ref_lz:10.2f} {nemo_lz:10.2f} {nemo_lz/ref_lz:7.2f}x")

    print(f"\n  참고 (이전 GRAZ 방식):")
    print(f"  {'szprod':16s} {ref_sz:10.2f} {'1524.04':>10s} {'36.67x':>8s}")
    print(f"  {'lzprod':16s} {ref_lz:10.2f} {'340.55':>10s} {'9.73x':>8s}")

except Exception as e:
    print(f"  Reference 파일 읽기 실패: {e}")

print("\n" + "=" * 60)
print("디버깅 완료!")
print("=" * 60)
