#!/bin/bash
# =============================================================
# 00_preprocess_cdo.sh
# NEMO-FEISTY 출력에서 필요한 변수만 추출 + 10년 시간 자르기
#
# 사용법:
#   cd /data01/labdisk/sungjin/NEMO-FEISTY
#   bash 00_preprocess_cdo.sh
#
# 출력: ./subset/ 디렉토리에 전처리된 파일 생성
#
# NOTE: grid_T_static의 depth는 time_counter=0 문제로 값이 비어있음.
#       → Python 스크립트에서 thetao 3D mask로 bathymetry 재구성.
# =============================================================

set -e

# --- 설정 ---
INDIR="."
OUTDIR="./subset"
mkdir -p ${OUTDIR}

echo "=============================================="
echo " NEMO-FEISTY CDO 전처리"
echo " 입력: ${INDIR}"
echo " 출력: ${OUTDIR}"
echo "=============================================="

# =============================================================
# 1. grid_T 5d → thetao, sbt / 10년 (timestep 1-730)
#    noleap calendar, 5일 간격: 73 steps/yr × 10yr = 730 steps
#    thetao: 3D 온도 (Tp, Tm 계산용 + bathymetry 재구성용)
#    sbt: 해저 온도 (Tb)
# =============================================================
echo ""
echo "[1/3] grid_T 5d → thetao + sbt, 10년 추출..."
IN="${INDIR}/ORCA2_5d_00010101_00501231_grid_T.nc"

# 변수 추출 후 시간 자르기 (파이프)
TMP1="${OUTDIR}/_tmp_gridT_vars.nc"
OUT="${OUTDIR}/grid_T_10yr.nc"

cdo selvar,thetao,sbt ${IN} ${TMP1}
echo "  변수 추출 완료 (thetao, sbt)"

cdo seltimestep,1/730 ${TMP1} ${OUT}
rm -f ${TMP1}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 2. diad_T 1y → Heup, EPC100, GRAZ1, GRAZ2 / 10년 (timestep 1-10)
#    Heup: 유광층 깊이
#    EPC100: 100m export → Martin curve로 해저 flux 계산
#    GRAZ1/GRAZ2: micro/meso-zoo grazing → szprod/lzprod
# =============================================================
echo ""
echo "[2/3] diad_T 1y → Heup + EPC100 + GRAZ1 + GRAZ2, 10년..."
IN="${INDIR}/ORCA2_1y_00010101_00501231_diad_T.nc"
OUT="${OUTDIR}/diad_T_1y_10yr.nc"

cdo seltimestep,1/10 -selvar,Heup,EPC100,GRAZ1,GRAZ2 ${IN} ${OUT}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 3. diad_T 1m → 같은 변수 / 10년 (timestep 1-120)
#    (월별 데이터 백업)
# =============================================================
echo ""
echo "[3/3] diad_T 1m → Heup + EPC100 + GRAZ1 + GRAZ2, 10년..."
IN="${INDIR}/ORCA2_1m_00010101_00501231_diad_T.nc"
OUT="${OUTDIR}/diad_T_1m_10yr.nc"

cdo seltimestep,1/120 -selvar,Heup,EPC100,GRAZ1,GRAZ2 ${IN} ${OUT}
echo "  → $(du -h ${OUT} | cut -f1)  ${OUT}"

# =============================================================
# 요약
# =============================================================
echo ""
echo "=============================================="
echo " 전처리 완료!"
echo "=============================================="
echo ""
echo " 생성된 파일:"
ls -lh ${OUTDIR}/*.nc
echo ""
echo " NOTE: bathymetry는 static 파일 문제로 thetao에서 재구성됩니다."
echo ""
echo " 다음 단계:"
echo "   python 02_extract_feisty_input.py ${OUTDIR}"
