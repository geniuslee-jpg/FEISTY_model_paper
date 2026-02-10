#!/bin/bash

set -e

BASE_DIR="/data01/labdisk/sungjin/NEMO-FEISTY"
ORG_DIR="${BASE_DIR}/01.org"
OUT_DIR="${BASE_DIR}/02.subset"
DOMCFG="${ORG_DIR}/ORCA_R2_zps_domcfg.nc"
TARGET="r360x180"   # 1도 x 1도

mkdir -p "${OUT_DIR}"

# =============================================================
# 1. 변수 + 기간 서브셋 (ORCA2 원본 격자)
# =============================================================
echo "=== [1] Variable & time subsetting ==="

# grid_T (5d): thetao + sbt → 10년분 (730 steps = 73/yr × 10yr)
cdo seltimestep,1/730 -selvar,thetao,sbt \
    "${ORG_DIR}/ORCA2_5d_00010101_00501231_grid_T.nc" \
    "${OUT_DIR}/grid_T_10yr.nc"

# diad_T (1y): Heup + EPC100 → 10년분
cdo seltimestep,1/10 -selvar,Heup,EPC100 \
    "${ORG_DIR}/ORCA2_1y_00010101_00501231_diad_T.nc" \
    "${OUT_DIR}/diad_T_1y_10yr.nc"

# diad_T (1m): Heup + EPC100 → 10년분 (120 steps = 12/yr × 10yr)
cdo seltimestep,1/120 -selvar,Heup,EPC100 \
    "${ORG_DIR}/ORCA2_1m_00010101_00501231_diad_T.nc" \
    "${OUT_DIR}/diad_T_1m_10yr.nc"

# ptrc_T (1y): ZOO + ZOO2 → closure term (mzrat × ZOO²) 계산용
cdo selvar,ZOO,ZOO2 \
    "${ORG_DIR}/ORCA2_1y_00010101_00501231_ptrc_T.nc" \
    "${OUT_DIR}/ptrc_T_1y_zoo.nc"

echo ""
echo "=== [2] Regrid to ${TARGET} ==="

# grid_T → 1도 (bilinear)
cdo remapbil,${TARGET} \
    "${OUT_DIR}/grid_T_10yr.nc" \
    "${OUT_DIR}/grid_T_10yr_1deg.nc"
echo "  grid_T done"

# diad_T 1m → 1도
cdo remapbil,${TARGET} \
    "${OUT_DIR}/diad_T_1m_10yr.nc" \
    "${OUT_DIR}/diad_T_1m_10yr_1deg.nc"
echo "  diad_T_1m done"

# diad_T 1y → 1도
cdo remapbil,${TARGET} \
    "${OUT_DIR}/diad_T_1y_10yr.nc" \
    "${OUT_DIR}/diad_T_1y_10yr_1deg.nc"
echo "  diad_T_1y done"

# =============================================================
# 3. domcfg regrid (bottom_level, e3t_0)
# =============================================================
echo ""
echo "=== [3] domcfg regrid ==="

# domcfg에 격자 정보 없을 수 있으므로 grid_T에서 복사
GRID_DESC="${OUT_DIR}/_orca2_grid.txt"
cdo griddes "${OUT_DIR}/grid_T_10yr.nc" > "${GRID_DESC}"

# bottom_level → nearest neighbor (정수 보존)
cdo remapnn,${TARGET} -setgrid,${GRID_DESC} \
    -selvar,bottom_level "${DOMCFG}" \
    "${OUT_DIR}/bottom_level_1deg.nc"
echo "  bottom_level done (remapnn)"

# e3t_0 → bilinear
cdo remapbil,${TARGET} -setgrid,${GRID_DESC} \
    -selvar,e3t_0 "${DOMCFG}" \
    "${OUT_DIR}/e3t_0_1deg.nc"
echo "  e3t_0 done (remapbil)"

echo ""
echo "=== Done ==="
ls -lh "${OUT_DIR}"
