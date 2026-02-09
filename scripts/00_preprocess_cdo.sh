#!/bin/bash
# =============================================================
# 00_preprocess_cdo.sh
# NEMO-FEISTY 출력에서 필요한 변수만 추출 + 10년 시간 자르기
#
# 사용법:
#   bash 00_preprocess_cdo.sh
#
# NOTE: grid_T_static의 depth는 time_counter=0 문제로 값이 비어있음.
#       → Python 스크립트에서 thetao 3D mask로 bathymetry 재구성.
# =============================================================

set -e

BASE_DIR="/data01/labdisk/sungjin/NEMO-FEISTY"
ORG_DIR="${BASE_DIR}/01.org"
OUT_DIR="${BASE_DIR}/02.subset"

mkdir -p "${OUT_DIR}"

echo "=============================================="
echo " NEMO-FEISTY CDO 전처리"
echo " 원본: ${ORG_DIR}"
echo " 출력: ${OUT_DIR}"
echo "=============================================="

# static (참고용 복사, depth는 Python에서 thetao로 재구성)
echo ""
echo "[1/4] grid_T_static 복사..."
cp "${ORG_DIR}/ORCA2_1y_00010101_00501231_grid_T_static.nc" \
   "${OUT_DIR}/static_depth.nc"

# grid_T 5d → thetao + sbt, 10년 (73 steps/yr × 10 = 730)
echo ""
echo "[2/4] grid_T 5d → thetao + sbt, 10년..."
cdo seltimestep,1/730 -selvar,thetao,sbt \
    "${ORG_DIR}/ORCA2_5d_00010101_00501231_grid_T.nc" \
    "${OUT_DIR}/grid_T_10yr.nc"

# diad_T 1y → 10년
echo ""
echo "[3/4] diad_T 1y → Heup + EPC100 + GRAZ1 + GRAZ2, 10년..."
cdo seltimestep,1/10 -selvar,Heup,EPC100,GRAZ1,GRAZ2 \
    "${ORG_DIR}/ORCA2_1y_00010101_00501231_diad_T.nc" \
    "${OUT_DIR}/diad_T_1y_10yr.nc"

# diad_T 1m → 10년 (120개월)
echo ""
echo "[4/4] diad_T 1m → Heup + EPC100 + GRAZ1 + GRAZ2, 10년..."
cdo seltimestep,1/120 -selvar,Heup,EPC100,GRAZ1,GRAZ2 \
    "${ORG_DIR}/ORCA2_1m_00010101_00501231_diad_T.nc" \
    "${OUT_DIR}/diad_T_1m_10yr.nc"

echo ""
echo "=============================================="
echo " 완료!"
echo "=============================================="
ls -lh "${OUT_DIR}"
echo ""
echo " 다음: python 02_extract_feisty_input.py ${OUT_DIR}"
