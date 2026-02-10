#!/usr/bin/env python3
"""
Closure term (mzrat × ZOO²) 디버깅 스크립트
ZOO/ZOO2 값 확인 + closure term 추정 + Reference 비교
"""
import xarray as xr
import numpy as np

BASE_DIR = "/data01/labdisk/sungjin/NEMO-FEISTY"
zoo_file = f"{BASE_DIR}/01.org/ORCA2_1y_00010101_00501231_ptrc_T.nc"
domcfg_file = f"{BASE_DIR}/01.org/ORCA_R2_zps_domcfg.nc"
ref_file = f"{BASE_DIR}/03.feisty_input/Input_global.csv"

# namelist_pisces_ref 값
MZRAT = 0.02     # micro (/day/(mmol/m3))
MZRAT2 = 0.01    # meso  (/day/(mmol/m3))
CONV_CLOSURE = 1e-3 * 365 * 12.0 * 9.0  # = 39.42

print("=" * 60)
print("Closure term 디버깅")
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

# 표층 (레벨 0) 값
print(f"\n  표층 (레벨 0):")
print(f"    ZOO  mean: {np.nanmean(zoo_mean[0]):.6f} mmol/m3")
print(f"    ZOO2 mean: {np.nanmean(zoo2_mean[0]):.6f} mmol/m3")

# 레벨별 평균 (상위 10레벨)
print(f"\n  레벨별 ZOO 평균 (상위 10):")
for k in range(min(10, zoo_mean.shape[0])):
    z_val = np.nanmean(zoo_mean[k])
    z2_val = np.nanmean(zoo2_mean[k])
    print(f"    lev {k:2d}: ZOO={z_val:.6f}  ZOO2={z2_val:.6f}")

ds.close()

# ─────────────────────────────────────────────
# [2] Closure term 계산 (단일 격자점 예시)
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[2] Closure term 단일 격자점 예시")

# 표층 대표값
zoo_ex = np.nanmean(zoo_mean[0])
zoo2_ex = np.nanmean(zoo2_mean[0])

closure_ex = MZRAT * zoo_ex**2
closure2_ex = MZRAT2 * zoo2_ex**2

print(f"  ZOO = {zoo_ex:.4f} mmol/m3")
print(f"  mzrat × ZOO² = {MZRAT} × {zoo_ex:.4f}² = {closure_ex:.6f} mmol/m3/day")
print(f"")
print(f"  ZOO2 = {zoo2_ex:.4f} mmol/m3")
print(f"  mzrat2 × ZOO2² = {MZRAT2} × {zoo2_ex:.4f}² = {closure2_ex:.6f} mmol/m3/day")

# 100m 적분 가정
depth_approx = 100.0  # m
szprod_ex = closure_ex * depth_approx * CONV_CLOSURE
lzprod_ex = closure2_ex * depth_approx * CONV_CLOSURE
print(f"\n  100m 적분 가정:")
print(f"    szprod ≈ {szprod_ex:.2f} gww/m2/yr")
print(f"    lzprod ≈ {lzprod_ex:.2f} gww/m2/yr")

# ─────────────────────────────────────────────
# [3] 전구 closure term 계산 (domcfg 사용)
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[3] 전구 closure term (수직적분)...")

ds_dom = xr.open_dataset(domcfg_file)
bottom_level = ds_dom["bottom_level"].values
e3t_0 = ds_dom["e3t_0"].values
ds_dom.close()

nz = e3t_0.shape[0]
ocean_mask = bottom_level > 0
k3d = np.arange(nz)[:, None, None]
tmask = (k3d < bottom_level[None, :, :]).astype(np.float64)

# closure term
closure_zoo = MZRAT * zoo_mean**2
closure_zoo2 = MZRAT2 * zoo2_mean**2

# 수직적분
szprod = np.nansum(closure_zoo * e3t_0 * tmask, axis=0) * CONV_CLOSURE
lzprod = np.nansum(closure_zoo2 * e3t_0 * tmask, axis=0) * CONV_CLOSURE
szprod[~ocean_mask] = np.nan
lzprod[~ocean_mask] = np.nan

print(f"\n  szprod (gww/m2/yr):")
print(f"    min:  {np.nanmin(szprod):.4f}")
print(f"    max:  {np.nanmax(szprod):.4f}")
print(f"    mean: {np.nanmean(szprod):.4f}")

print(f"\n  lzprod (gww/m2/yr):")
print(f"    min:  {np.nanmin(lzprod):.4f}")
print(f"    max:  {np.nanmax(lzprod):.4f}")
print(f"    mean: {np.nanmean(lzprod):.4f}")

# ─────────────────────────────────────────────
# [4] Reference와 비교
# ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("[4] Reference 비교")

try:
    import pandas as pd
    ref = pd.read_csv(ref_file)
    print(f"\n  {'':12s} {'Ref mean':>10s} {'NEMO mean':>10s} {'비율':>8s}")
    print(f"  {'-'*42}")

    ref_sz = ref["szprod"].mean()
    ref_lz = ref["lzprod"].mean()
    nemo_sz = np.nanmean(szprod)
    nemo_lz = np.nanmean(lzprod)

    print(f"  {'szprod':12s} {ref_sz:10.2f} {nemo_sz:10.2f} {nemo_sz/ref_sz:7.2f}x")
    print(f"  {'lzprod':12s} {ref_lz:10.2f} {nemo_lz:10.2f} {nemo_lz/ref_lz:7.2f}x")

    # 이전 GRAZ 방식과 비교
    print(f"\n  참고 (이전 GRAZ 방식):")
    print(f"  {'szprod':12s} {ref_sz:10.2f} {'1524.04':>10s} {'36.67x':>8s}")
    print(f"  {'lzprod':12s} {ref_lz:10.2f} {'340.55':>10s} {'9.73x':>8s}")

except Exception as e:
    print(f"  Reference 파일 읽기 실패: {e}")

print("\n" + "=" * 60)
print("디버깅 완료!")
print("=" * 60)
