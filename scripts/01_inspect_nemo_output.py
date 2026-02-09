"""
01_inspect_nemo_output.py
=========================
NEMO-FEISTY NetCDF 출력 파일의 변수, 차원, 속성을 탐색하는 스크립트.
각 파일에 어떤 변수가 들어있는지 확인한 후 02_extract 스크립트를 수정할 때 참고.

사용법:
    python 01_inspect_nemo_output.py [DATA_DIR]

기본 DATA_DIR: /data01/labdisk/sungjin/NEMO-FEISTY
"""

import sys
import os
import glob
import netCDF4 as nc

def inspect_file(filepath):
    """단일 NetCDF 파일의 변수 정보를 출력."""
    fname = os.path.basename(filepath)
    fsize_mb = os.path.getsize(filepath) / 1e6
    print(f"\n{'='*80}")
    print(f"FILE: {fname}  ({fsize_mb:.1f} MB)")
    print(f"{'='*80}")

    ds = nc.Dataset(filepath, 'r')

    # 차원 정보
    print(f"\n  Dimensions:")
    for dname, dim in ds.dimensions.items():
        print(f"    {dname:20s} = {len(dim):6d}  (unlimited: {dim.isunlimited()})")

    # 변수 정보
    print(f"\n  Variables:")
    for vname, var in ds.variables.items():
        dims_str = str(var.dimensions)
        shape_str = str(var.shape)
        long_name = getattr(var, 'long_name', '')
        units = getattr(var, 'units', '')
        print(f"    {vname:25s}  shape={shape_str:30s}  dims={dims_str}")
        if long_name:
            print(f"    {'':25s}  long_name: {long_name}")
        if units:
            print(f"    {'':25s}  units: {units}")

    # 시간 범위 확인
    if 'time_counter' in ds.variables:
        time_var = ds.variables['time_counter']
        t_units = getattr(time_var, 'units', 'unknown')
        t_cal = getattr(time_var, 'calendar', 'unknown')
        print(f"\n  Time info:")
        print(f"    units: {t_units}")
        print(f"    calendar: {t_cal}")
        print(f"    range: {time_var[0]:.1f} to {time_var[-1]:.1f}")
        print(f"    n_timesteps: {len(time_var)}")

    ds.close()


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else "/data01/labdisk/sungjin/NEMO-FEISTY"

    if not os.path.isdir(data_dir):
        print(f"ERROR: Directory not found: {data_dir}")
        sys.exit(1)

    nc_files = sorted(glob.glob(os.path.join(data_dir, "*.nc")))

    if not nc_files:
        print(f"ERROR: No .nc files found in {data_dir}")
        sys.exit(1)

    print(f"Found {len(nc_files)} NetCDF files in {data_dir}")
    print(f"Files: {[os.path.basename(f) for f in nc_files]}")

    for fpath in nc_files:
        try:
            inspect_file(fpath)
        except Exception as e:
            print(f"\n  ERROR reading {os.path.basename(fpath)}: {e}")

    # 요약: FEISTY 입력에 필요한 변수 매핑 가이드
    print(f"\n{'='*80}")
    print("FEISTY INPUT VARIABLE MAPPING GUIDE")
    print(f"{'='*80}")
    print("""
FEISTY Input_global.csv에 필요한 변수:
  lon, lat       <- grid_T 또는 grid_T_static의 nav_lon, nav_lat
  Tp (0-100m)    <- grid_T의 votemper/thetao, 0-100m 깊이 평균
  Tm (500-1500m) <- grid_T의 votemper/thetao, 500-1500m 깊이 평균
  Tb (bottom)    <- grid_T의 votemper/thetao, 최저층 값
  depth          <- grid_T_static의 bathymetry/mbathy/gdept_0
  photic         <- diad_T의 ZEU 또는 Chl로부터 계산
  lzprod         <- diad_T의 대형 동물플랑크톤 생산 관련 변수
  szprod         <- diad_T의 소형 동물플랑크톤 생산 관련 변수
  dfbot          <- diad_T의 해저 detrital flux 관련 변수

위 출력 결과를 보고 실제 변수명을 확인한 후,
02_extract_feisty_input.py 스크립트의 변수명을 수정하세요.
""")


if __name__ == "__main__":
    main()
